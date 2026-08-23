# TRAIN
**Toolbox for Reducing Atmospheric InSAR Noise**

TRAIN provides state-of-the-art tropospheric correction methods for InSAR time-series analysis.  
It is processor-independent and compatible with **StaMPS**, **MintPy**, **ISCE**, and ROI\_PAC-derived workflows.

```
We welcome community contributions and request users to contribute back to the repo.
```

---

## Methods

| Category | Method | Script |
|---|---|---|
| **Empirical** | Linear phase–elevation | `aps_linear.m` |
| **Empirical** | Spatially-variable power-law | `aps_powerlaw.m` |
| **NWM** | ERA5 / ERA5T (via PyAPS3) | `aps_weather_model('era5',…)` |
| **NWM** | ERA-Interim *(discontinued Aug 2019 — legacy only)* | `aps_weather_model('era',…)` |
| **NWM** | MERRA / MERRA2 | `aps_weather_model('merra',…)` |
| **NWM** | NARR | `aps_weather_model('narr',…)` |
| **NWM** | WRF (regional) | `aps_wrf.m` |
| **Spectrometer** | Envisat MERIS (via Python) | `aps_meris.m` |
| **Spectrometer** | MODIS | `aps_modis.m` |
| **Empirical** | GACOS | `aps_weather_model('gacos',…)` |

---

## Python dependencies

ERA5 download and MERIS processing now use Python scripts (`python_modules/`).

```bash
pip install pyaps3 cdsapi rasterio scipy numpy
```

Configure your **CDS API** key for ERA5 (see `python_modules/README`).  
ERA5T (near-real-time, last ~3 months) is served automatically — no extra flag needed.

---

## Installation

1. Add `matlab/` to your MATLAB path:
   ```matlab
   addpath /path/to/TRAIN/matlab
   ```
2. Install Python dependencies (see above).
3. Copy `APS_CONFIG.sh` (or `.tcsh`) and set `$STAMPS` and `$TRAIN_DIR`.

---

## Quick start

```matlab
% ERA5 correction (Steps 0–3: list files, download, SAR delays, InSAR delays)
aps_weather_model('era5', 0, 3)

% ERA5T (near-real-time) — same call, ERA5T served automatically by CDS
aps_weather_model('era5t', 0, 3)

% Empirical power-law correction
aps_powerlaw(1, 5)

% MERIS spectrometer correction
aps_meris(1, 3)

% GACOS correction
aps_weather_model('gacos', 0, 3)
```

---

## Citation

Please cite the primary TRAIN paper:

> Bekaert, D.P.S., Walters, R.J., Wright, T.J., Hooper, A.J., and Parker, D.J. (2015),
> *Statistical comparison of InSAR tropospheric correction techniques*,
> Remote Sensing of Environment, doi:[10.1016/j.rse.2015.08.035](https://doi.org/10.1016/j.rse.2015.08.035)

For the **power-law** method also cite:

> Bekaert, D.P.S., Hooper, A.J., and Wright, T.J. (2015),
> *A spatially-variable power-law tropospheric correction technique for InSAR data*,
> JGR Solid Earth, doi:[10.1029/2014JB011558](https://doi.org/10.1029/2014JB011558)

For the **ERA5 / PyAPS3 download** scripts (`python_modules/aps_era5_download.py`) also cite:

> Mohseni Aref, M. et al. (2026), contributions via MintPy/PyAPS3 integration,
> see [MintPy](https://github.com/insarlab/MintPy) and PyAPS3.

For **MERIS/OLCI** Python processing (`python_modules/aps_meris_pwv.py`) also cite:

> Mohseni Aref, M. (2026), Python MERIS/OLCI atmospheric correction for TRAIN.

---

## Zenodo

A citable Zenodo archive of this release is available at:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

*(Replace `XXXXXXX` with the DOI assigned after uploading a new release to [zenodo.org](https://zenodo.org).)*

**To create a new Zenodo release** with these Python updates:
1. Tag this commit: `git tag -a v2.0 -m "Python ERA5/MERIS update via PyAPS3"`
2. Push the tag: `git push origin v2.0`
3. On zenodo.org, link the GitHub repo and publish the release.
4. Add yourself (M. Mohseni Aref) as a contributor/author in the Zenodo metadata.
5. Replace `XXXXXXX` above with the assigned DOI.

---

## Acknowledgements

Thanks to Richard J. Walters, Hannes Bathke, Simran Sangha, Tim J. Wright, Andy J. Hooper,
Doug J. Parker, Zhenhong Li, the Leeds InSAR group, COMET members, and the community
for their feedback and contributions.

Python ERA5 download and MERIS processing contributed by Mohammad Mohseni Aref (2026),
integrated via [MintPy](https://github.com/insarlab/MintPy) PyAPS3.

---

## License

TRAIN is distributed under a **GNU GPL** licence.

