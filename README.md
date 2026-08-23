# TRAIN
**Toolbox for Reducing Atmospheric InSAR Noise — v3.1-alpha**

TRAIN provides state-of-the-art tropospheric correction methods for InSAR time-series analysis.  
It is processor-independent and compatible with **StaMPS**, **ISCE**, and ROI\_PAC-derived workflows.

This repository (v3.1-alpha, M. Mohseni Aref) extends **Version 3beta** by D. Bekaert (University of Leeds / JPL),
distributed under a **GNU GPL** licence. The extension adds Python-based ERA5/ERA5T download and
MERIS/OLCI processing; the core MATLAB toolbox and licence remain unchanged.

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

For the **Python ERA5 download and MERIS processing** scripts also cite this repository (v3.1-alpha, doi in Zenodo section below).

---

## Zenodo

A citable archive of this toolbox (v3.1-alpha) is available on Zenodo:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

*(Update the DOI badge once a new release is published on [zenodo.org](https://zenodo.org).)*

---

## License

TRAIN is distributed under a **GNU GPL** licence.

