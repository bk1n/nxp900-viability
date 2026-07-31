# nxp900-viability

Processing of NXP900 viability data from GDSC, Oncolines, Temps, and Brognard.

Produces fully processed viability dataset in `out/viability.csv` for downstream analysis.

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