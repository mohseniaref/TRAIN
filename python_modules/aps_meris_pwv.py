#!/usr/bin/env python3
"""Process Envisat MERIS / Sentinel-3 OLCI water-vapour products for TRAIN.

Python replacement for aps_meris_SAR.m and aps_meris_InSAR.m.

Reads per-date GeoTIFF files containing MERIS IWV (integrated water vapour),
applies cloud/quality masking, fills gaps by Gaussian interpolation and
reprojects to the InSAR geometry grid, then saves the LOS delay in the same
.xyz format expected by the rest of the TRAIN pipeline.

MERIS IWV → LOS slant wet delay conversion
    SWD [m] = conversion * IWV [kg/m²] / cos(inc_angle [rad])

where the default conversion factor (Pi = 6.2e-4 m / (kg m⁻²)) follows
Bekaert et al. (2015) and is calibrated from sounding data in TRAIN.

Usage:
    python aps_meris_pwv.py \\
        --meris-files meris_files.txt \\
        --dem  dem.rsc \\           # or --lat-range ... --lon-range ...
        --inc-angle 23.0 \\
        --conversion 6.2e-4 \\
        --outdir ./meris_delays \\
        --res  0.001                # output resolution in degrees

Sentinel-3 OLCI mode (--sensor olci):
    Same interface; OLCI band 19 carries the IWV column in [kg/m²].

Output per date:
    <outdir>/<YYYYMMDD>_SWD.xyz   (lon lat SWD[m])
"""

import argparse
import os
import sys
import warnings

import numpy as np


##############################################################################
# Helpers
##############################################################################

def _import_rasterio():
    try:
        import rasterio
        from rasterio.warp import reproject, Resampling
        return rasterio, reproject, Resampling
    except ImportError:
        sys.exit('[ERROR] rasterio is required: pip install rasterio')


def _import_scipy():
    try:
        from scipy.ndimage import gaussian_filter
        from scipy.interpolate import RegularGridInterpolator
        return gaussian_filter, RegularGridInterpolator
    except ImportError:
        sys.exit('[ERROR] scipy is required: pip install scipy')


##############################################################################
# MERIS / OLCI cloud masking
##############################################################################

def _meris_cloud_mask(band_flags):
    """Return boolean mask (True = valid) from MERIS quality/flag band (band index 33).

    Bit definitions follow the MERIS Level-2 flag coding:
        bit 1  (value   2) : cloud
        bit 22 (value 4194304) : low pressure
        bit 23 (value 8388608) : p-confidence

    Args:
        band_flags : 2-D uint32 array

    Returns:
        valid : bool array, True where pixel is cloud-free and confident
    """
    cloud   = (band_flags & (1 << 1))  != 0
    low_p   = (band_flags & (1 << 22)) != 0
    low_pco = (band_flags & (1 << 23)) != 0
    invalid = cloud | low_p | low_pco
    return ~invalid


def _olci_cloud_mask(quality_flags):
    """Return boolean valid mask from Sentinel-3 OLCI quality flag.

    Bit 1 (INVALID) and bit 10 (CLOUD) are masked out.
    """
    invalid = (quality_flags & (1 << 1))  != 0
    cloud   = (quality_flags & (1 << 10)) != 0
    return ~(invalid | cloud)


##############################################################################
# Gap filling
##############################################################################

def fill_gaps(data, sigma_px=5):
    """Fill NaN gaps by Gaussian-weighted interpolation.

    1. Replace NaNs with the Gaussian-blurred field (computed from valid pixels).
    2. Iterate twice to reduce boundary artifacts.
    """
    gaussian_filter, _ = _import_scipy()

    filled = data.copy()
    for _ in range(2):
        nan_mask = np.isnan(filled)
        if not nan_mask.any():
            break
        tmp = np.where(nan_mask, 0.0, filled)
        weight = np.where(nan_mask, 0.0, 1.0)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            smooth_val = gaussian_filter(tmp, sigma=sigma_px)
            smooth_wgt = gaussian_filter(weight, sigma=sigma_px)
        smooth_wgt = np.where(smooth_wgt < 1e-6, np.nan, smooth_wgt)
        filled[nan_mask] = (smooth_val / smooth_wgt)[nan_mask]
    return filled


##############################################################################
# Resampling
##############################################################################

def resample_to_grid(src_data, src_lats, src_lons, tgt_lats, tgt_lons):
    """Bilinear resample src_data from (src_lats, src_lons) to (tgt_lats, tgt_lons).

    src_lats / src_lons must be 1-D and monotonically increasing.
    tgt_lats / tgt_lons are 2-D output grids.
    """
    _, RegularGridInterpolator = _import_scipy()

    interp = RegularGridInterpolator(
        (src_lats, src_lons),
        src_data,
        method='linear',
        bounds_error=False,
        fill_value=np.nan,
    )
    pts = np.column_stack([tgt_lats.ravel(), tgt_lons.ravel()])
    out = interp(pts).reshape(tgt_lats.shape)
    return out.astype(np.float32)


##############################################################################
# IWV → LOS SWD conversion
##############################################################################

def iwv_to_swd(iwv, conversion, inc_rad):
    """Convert IWV [kg/m²] to slant wet delay [m].

    SWD = conversion * IWV / cos(inc_angle)
    """
    cos_inc = np.cos(inc_rad)
    cos_inc[cos_inc < 0.1] = np.nan     # guard against near-90° incidence
    return conversion * iwv / cos_inc


##############################################################################
# Core per-date processing
##############################################################################

def process_meris_date(meris_tif, tgt_lats, tgt_lons, inc_angle_deg,
                       conversion=6.2e-4, sensor='meris', iwv_band=1,
                       flag_band=None):
    """Read one MERIS/OLCI GeoTIFF, mask, fill gaps and resample to InSAR grid.

    Parameters
    ----------
    meris_tif     : str, path to GeoTIFF
    tgt_lats      : 2-D array, target latitude grid [°]
    tgt_lons      : 2-D array, target longitude grid [°]
    inc_angle_deg : float or 2-D array, incidence angle [°]
    conversion    : float, IWV→SWD conversion factor [m/(kg m⁻²)]
    sensor        : 'meris' or 'olci'
    iwv_band      : int (1-based), band index containing IWV
    flag_band     : int (1-based) or None; band containing quality/flag data

    Returns
    -------
    swd : 2-D float32 array on the InSAR grid [m], NaN where no valid data
    """
    rasterio, _, Resampling = _import_rasterio()

    with rasterio.open(meris_tif) as ds:
        # ---- read IWV band ----
        iwv_raw = ds.read(iwv_band).astype(np.float32)
        profile = ds.profile
        transform = ds.transform

        # pixel centres
        rows, cols = np.meshgrid(np.arange(ds.height), np.arange(ds.width), indexing='ij')
        src_lons_2d, src_lats_2d = rasterio.transform.xy(transform, rows, cols)
        src_lats_2d = np.array(src_lats_2d, dtype=np.float64)
        src_lons_2d = np.array(src_lons_2d, dtype=np.float64)

        # ---- build cloud / quality mask ----
        if flag_band is not None:
            flags = ds.read(flag_band).astype(np.uint32)
        elif ds.count >= 33 and sensor == 'meris':
            # MERIS standard product: band 33 is the flag band
            flags = ds.read(33).astype(np.uint32)
        else:
            flags = None

        if flags is not None:
            if sensor == 'olci':
                valid = _olci_cloud_mask(flags)
            else:
                valid = _meris_cloud_mask(flags)
        else:
            # no flag band: mask only fill values (0 or NaN)
            valid = (iwv_raw > 0) & np.isfinite(iwv_raw)

        iwv = iwv_raw.astype(np.float64)
        iwv[~valid] = np.nan

    # ---- scale factor for MERIS standard product (IWV in g/cm² → kg/m²) ----
    # MERIS L2 IWV is stored as uint16 with scale 0.001 g/cm² ; here we assume
    # the caller passes a pre-scaled float array.  Adjust if raw uint16:
    # iwv *= 0.1   # g/cm² → kg/m²  (check sensor product manual)

    # ---- extend cloud mask edges (dilate) ----
    try:
        from scipy.ndimage import binary_dilation
        struct = np.ones((3, 3), dtype=bool)
        nan_mask = np.isnan(iwv)
        dilated = binary_dilation(nan_mask, structure=struct)
        iwv[dilated] = np.nan
    except ImportError:
        pass

    # ---- fill gaps ----
    iwv_filled = fill_gaps(iwv)

    # ---- get 1-D lat/lon vectors for regular-grid interpolation ----
    # Use one coordinate per row/col (do not sort independently of the data)
    src_lat_1d = np.round(src_lats_2d[:, 0], 5)
    src_lon_1d = np.round(src_lons_2d[0, :], 5)

    if src_lat_1d[0] > src_lat_1d[-1]:
        src_lat_1d = src_lat_1d[::-1]
        iwv_filled = iwv_filled[::-1, :]
    if src_lon_1d[0] > src_lon_1d[-1]:
        src_lon_1d = src_lon_1d[::-1]
        iwv_filled = iwv_filled[:, ::-1]
    # ---- resample to InSAR grid ----
    iwv_insar = resample_to_grid(iwv_filled, src_lat_1d, src_lon_1d,
                                  tgt_lats, tgt_lons)

    # ---- IWV → LOS SWD ----
    inc_rad = np.deg2rad(np.full_like(tgt_lats, inc_angle_deg, dtype=np.float32)
                         if np.isscalar(inc_angle_deg)
                         else inc_angle_deg)
    swd = iwv_to_swd(iwv_insar, conversion, inc_rad).astype(np.float32)
    return swd


##############################################################################
# Output
##############################################################################

def save_xyz(outfile, lats_2d, lons_2d, data_2d):
    """Save a delay map as whitespace-separated lon lat value lines."""
    mask = np.isfinite(data_2d)
    lons_f = lons_2d[mask].ravel()
    lats_f = lats_2d[mask].ravel()
    vals_f = data_2d[mask].ravel()
    with open(outfile, 'w') as fh:
        for lo, la, v in zip(lons_f, lats_f, vals_f):
            fh.write(f'{lo:.6f} {la:.6f} {v:.6f}\n')
    print(f'[INFO] saved {outfile}')


##############################################################################
# Build target lat/lon grid from .rsc or bounding-box
##############################################################################

def build_target_grid(rsc_file=None,
                      lat_range=None, lon_range=None, res=0.001):
    """Return (lat_2d, lon_2d) on the output grid."""
    if rsc_file and os.path.isfile(rsc_file):
        meta = _read_rsc(rsc_file)
        lon0 = float(meta['X_FIRST'])
        lat0 = float(meta['Y_FIRST'])
        dlon = float(meta['X_STEP'])
        dlat = float(meta['Y_STEP'])
        nrow = int(meta['LENGTH'])
        ncol = int(meta['WIDTH'])
        lats = lat0 + dlat * np.arange(nrow)
        lons = lon0 + dlon * np.arange(ncol)
    elif lat_range and lon_range:
        lats = np.arange(lat_range[0], lat_range[1], res)
        lons = np.arange(lon_range[0], lon_range[1], res)
    else:
        raise ValueError('Provide either --dem (with .rsc) or --lat-range / --lon-range.')
    return np.meshgrid(lats, lons, indexing='ij')   # (nrow, ncol)


def _read_rsc(rsc_file):
    meta = {}
    with open(rsc_file) as fh:
        for line in fh:
            kv = line.strip().split()
            if len(kv) >= 2:
                meta[kv[0].upper()] = kv[1]
    return meta


##############################################################################
# CLI
##############################################################################

def create_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--meris-files', required=True,
                   help='Text file listing full paths to per-date MERIS/OLCI GeoTIFFs '
                        '(one file per line; first line may be "files").')
    p.add_argument('--dem', default=None,
                   help='.rsc file describing the output grid (WIDTH, LENGTH, X_FIRST, etc.).')
    p.add_argument('--lat-range', nargs=2, type=float, metavar=('LAT_MIN', 'LAT_MAX'),
                   help='Lat range [°] for output grid (used when --dem is not given).')
    p.add_argument('--lon-range', nargs=2, type=float, metavar=('LON_MIN', 'LON_MAX'),
                   help='Lon range [°] for output grid.')
    p.add_argument('--res', type=float, default=0.001,
                   help='Output grid resolution in degrees (default 0.001 ≈ 100 m).')
    p.add_argument('--inc-angle', type=float, default=23.0,
                   help='Mean incidence angle in degrees (default 23.0).')
    p.add_argument('--conversion', type=float, default=6.2e-4,
                   help='IWV→SWD conversion factor [m/(kg m⁻²)] (default 6.2e-4).')
    p.add_argument('--sensor', choices=['meris', 'olci'], default='meris',
                   help='Sensor type: meris (Envisat) or olci (Sentinel-3). default: meris')
    p.add_argument('--iwv-band', type=int, default=1,
                   help='1-based band index for IWV in the GeoTIFF (default 1).')
    p.add_argument('--flag-band', type=int, default=None,
                   help='1-based band index for quality flags. Auto-detected for MERIS (band 33).')
    p.add_argument('--outdir', default='./meris_delays',
                   help='Output directory for *_SWD.xyz files.')
    return p


def main(argv=None):
    p = create_parser()
    inps = p.parse_args(argv)

    # ---- read file list ----
    with open(inps.meris_files) as fh:
        lines = [l.strip() for l in fh if l.strip() and l.strip().lower() != 'files']
    print(f'[INFO] {len(lines)} MERIS/OLCI files to process.')

    # ---- build target grid ----
    rsc = inps.dem if inps.dem and inps.dem.endswith('.rsc') else (
          inps.dem + '.rsc' if inps.dem and not inps.dem.endswith('.rsc') else None)
    lats_2d, lons_2d = build_target_grid(
        rsc_file=rsc,
        lat_range=inps.lat_range,
        lon_range=inps.lon_range,
        res=inps.res,
    )

    os.makedirs(inps.outdir, exist_ok=True)

    # ---- process each date ----
    for tif_path in lines:
        tif_path = tif_path.strip()
        if not os.path.isfile(tif_path):
            print(f'[WARN] file not found, skipping: {tif_path}')
            continue

        # extract YYYYMMDD from the parent folder name (TRAIN convention) or filename
        import re
        date_match = re.search(r'(\d{8})', os.path.basename(os.path.dirname(tif_path)))
        if not date_match:
            date_match = re.search(r'(\d{8})', os.path.basename(tif_path))
        date_str = date_match.group(1) if date_match else 'unknown'

        outfile = os.path.join(inps.outdir, f'{date_str}_SWD.xyz')
        if os.path.isfile(outfile):
            print(f'[INFO] already exists, skipping: {outfile}')
            continue

        print(f'[INFO] processing {date_str}: {tif_path}')
        try:
            swd = process_meris_date(
                tif_path,
                lats_2d, lons_2d,
                inc_angle_deg=inps.inc_angle,
                conversion=inps.conversion,
                sensor=inps.sensor,
                iwv_band=inps.iwv_band,
                flag_band=inps.flag_band,
            )
            save_xyz(outfile, lats_2d, lons_2d, swd)
        except Exception as exc:
            print(f'[ERROR] failed for {date_str}: {exc}')
            continue

    print('[INFO] MERIS/OLCI processing complete.')


if __name__ == '__main__':
    main()
