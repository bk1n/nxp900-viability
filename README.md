# nxp900-viability

An end-to-end dose-response modelling pipeline to pre-process, QC, visualise, and model data from large-scale viability screens (GDSC and Oncolines) against NXP900.

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

## Run

```R
Rscript drm.R
```

Response metrics/parameters are saved to `out/viability.csv` for downstream analysis.
