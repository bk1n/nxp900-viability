# nxp900-viability

Processing of NXP900 viability data from GDSC, Oncolines, Temps, and Brognard.
Runs QC, visualisation, normalisation, and dose-response curve fitting and saves response metrics/parameters to `out/viability.csv` for downstream analysis.

Install packages, pull data, then run: `Rscript drm.R`

## Packages

Packages are managed with `renv`

```R
renv::restore() # to install
```

## Pull data

Create a `.env` file in root with data path and run `fetch_data.sh`:

```bash
DATA_PATH="/full/path/to/data"
```
