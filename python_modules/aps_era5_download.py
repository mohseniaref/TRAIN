#!/usr/bin/env python3
"""Download ERA5 / ERA5T grib files for TRAIN using PyAPS3 (MintPy).

Replaces aps_era5_ECMWF_Python.m and the download portion of aps_era5_files.m.
ERA5T (near-real-time, last ~3 months) is served automatically by the CDS API
when final ERA5 is not yet available; no separate model flag is needed.

Usage (called from MATLAB via system() or standalone):
    python aps_era5_download.py --dates 20200101 20200215 \
        --hour 06 --snwe 30 50 -10 40 --outdir /path/to/ERA5

Requirements:
    pyaps3  (MintPy branch with CDS support)
    cdsapi  ~/.cdsapirc with valid UID:key
"""

import argparse
import os
import re
import sys

import numpy as np


##############################################################################
WEATHER_MODEL_HOURS = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
                       10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]


def create_parser():
    p = argparse.ArgumentParser(
        description='Download ERA5/ERA5T grib files for TRAIN via PyAPS3',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument('--dates', nargs='+', required=True,
                   help='SAR acquisition dates in YYYYMMDD format')
    p.add_argument('--hour', default=None,
                   help='UTC hour (2-digit, zero-padded). Auto-detected from --utc-sec if omitted.')
    p.add_argument('--utc-sec', dest='utc_sec', type=float, default=None,
                   help='Satellite pass time in UTC seconds (used to pick closest ERA5 hour).')
    p.add_argument('--snwe', nargs=4, type=float, default=None,
                   metavar=('S', 'N', 'W', 'E'),
                   help='Bounding box in degrees. Rounded to nearest 10° for ERA5.')
    p.add_argument('--outdir', default='.', help='Output directory for grib files.')
    p.add_argument('--check-only', action='store_true',
                   help='Print required file list without downloading.')
    p.add_argument('--era5t', action='store_true',
                   help='Force ERA5T (near-real-time) dataset. Normally auto-detected by CDS.')
    return p


##############################################################################
def closest_hour(utc_sec):
    """Return zero-padded hour string closest to utc_sec seconds."""
    h = int(min(WEATHER_MODEL_HOURS, key=lambda x: abs(x - utc_sec / 3600.0)))
    return f'{h:02d}'


def snwe_to_str(snwe):
    """Convert (S, N, W, E) floats to the TRAIN/MintPy file-name fragment."""
    s, n, w, e = [int(v) for v in snwe]
    area = ''
    area += f'_S{abs(s)}' if s < 0 else f'_N{abs(s)}'
    area += f'_S{abs(n)}' if n < 0 else f'_N{abs(n)}'
    area += f'_W{abs(w)}' if w < 0 else f'_E{abs(w)}'
    area += f'_W{abs(e)}' if e < 0 else f'_E{abs(e)}'
    return area


def round_snwe(snwe, step=10):
    """Round bounding box outward to multiples of `step` degrees (ERA5 convention)."""
    s, n, w, e = snwe
    s = int(np.floor(s / step) * step)
    n = int(np.ceil(n / step) * step)
    w = int(np.floor(w / step) * step)
    e = int(np.ceil(e / step) * step)
    return (s, n, w, e)


def build_grib_filenames(dates, hour, outdir, snwe=None):
    """Return list of (date, expected_grib_path) pairs."""
    area = snwe_to_str(snwe) if snwe else ''
    pairs = []
    for d in dates:
        fname = f'ERA5{area}_{d}_{hour}.grb'
        pairs.append((d, os.path.join(outdir, fname)))
    return pairs


def check_existing(grib_files):
    """Return (existing, missing) lists; flag files < 90 % of median size as corrupt."""
    existing = [f for f in grib_files if os.path.isfile(f)]
    if not existing:
        return [], grib_files

    sizes = [os.path.getsize(f) for f in existing]
    median_sz = float(np.median(sizes))
    corrupt = [f for f, sz in zip(existing, sizes) if sz < 0.9 * median_sz]
    for f in corrupt:
        print(f'[WARN] removing corrupt file ({os.path.getsize(f)} B): {f}')
        os.remove(f)
        existing.remove(f)

    missing = [f for f in grib_files if f not in existing]
    return existing, missing


##############################################################################
def download_era5(dates_missing, hour, outdir, snwe, era5t=False):
    """Call PyAPS3 ECMWFdload for the missing dates."""
    try:
        import pyaps3 as pa
    except ImportError:
        sys.exit('[ERROR] pyaps3 is not installed. Install the MintPy PyAPS3 package.')

    print(f'[INFO] PyAPS3 version: {pa.__version__}')

    snwe_int = tuple(int(v) for v in snwe) if snwe else None
    flist = build_grib_filenames(dates_missing, hour, outdir, snwe=snwe_int)
    flist_paths = [p for _, p in flist]

    os.makedirs(outdir, exist_ok=True)

    kwargs = dict(model='ERA5', snwe=snwe_int, flist=flist_paths)
    if era5t:
        # Pass the era5t flag through to PyAPS if the installed version supports it
        kwargs['era5t'] = True

    print(f'[INFO] Downloading {len(dates_missing)} ERA5 file(s) to {outdir} ...')
    for attempt in range(1, 4):
        try:
            pa.ECMWFdload(dates_missing, hour, outdir, **kwargs)
            break
        except Exception as exc:
            if attempt < 3:
                print(f'[WARN] Attempt {attempt} failed: {exc}. Retrying ...')
            else:
                print(f'[ERROR] Download failed after 3 attempts: {exc}')
                print('[INFO] Continuing with whatever was downloaded.')


##############################################################################
def main(argv=None):
    p = create_parser()
    inps = p.parse_args(argv)

    # ---- Resolve UTC hour ----
    if inps.hour is None:
        if inps.utc_sec is not None:
            inps.hour = closest_hour(inps.utc_sec)
            print(f'[INFO] Closest ERA5 hour to {inps.utc_sec:.0f} s UTC: {inps.hour}:00')
        else:
            p.error('Provide --hour or --utc-sec.')

    # ---- Round bounding box ----
    snwe_rounded = round_snwe(inps.snwe) if inps.snwe else None
    if snwe_rounded:
        print(f'[INFO] SNWE (rounded to 10°): {snwe_rounded}')

    # ---- Build expected file list ----
    file_pairs = build_grib_filenames(inps.dates, inps.hour, inps.outdir, snwe=snwe_rounded)
    all_files = [p for _, p in file_pairs]

    if inps.check_only:
        print('[INFO] Required ERA5 files:')
        for f in all_files:
            print(f'  {os.path.basename(f)}')
        return

    # ---- Download only missing / non-corrupt files ----
    existing, missing = check_existing(all_files)
    print(f'[INFO] Already present : {len(existing)}')
    print(f'[INFO] To download     : {len(missing)}')

    if missing:
        dates_missing = [re.findall(r'\d{8}', os.path.basename(f))[0] for f in missing]
        download_era5(dates_missing, inps.hour, inps.outdir, snwe_rounded, era5t=inps.era5t)

    # ---- Final check ----
    existing_after, still_missing = check_existing(all_files)
    if still_missing:
        print(f'[WARN] {len(still_missing)} file(s) still missing after download:')
        for f in still_missing:
            print(f'  {f}')
    else:
        print('[INFO] All ERA5 files present.')


if __name__ == '__main__':
    main()
