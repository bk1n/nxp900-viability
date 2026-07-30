# Fit dose-response curves per screen, combine into viability_raw.csv.
#
# Each entry of SCREENS describes one screen and how to fit it:
#   path         - directory holding raw_processing output and where DRC artefacts land.
#                  Expects a drc_ready.csv in there (written by raw_processing/<screen>.R);
#                  writes <RESPONSE>_models.qs2 and drc_metrics.csv into the same dir.
#   responses    - subset of c("IC50","GI50","GR50") to fit (depends on what the screen supports)
#   top_dose, bottom_dose - dose range (uM) used for AOC / DSS integration
#   combine_name - column-prefix used by combine_drcfits; screens sharing one are row-bound
#   grouped      - TRUE for screens with multiple dose ranges (handled per group)
#
# See preprocessing/viability/README.md for the output naming conventions.
#
# Toggles below are character vectors of screen keys. All default to nothing
# because per-screen fits take hours — list a screen explicitly to do work.

devtools::load_all(".")

library(parallel)

DATA_PATH <- "data/"
OUT_PATH <- "out/"

# LOAD METADATA
metadata <- dplyr::read_csv(file.path(DATA_PATH, "Model.csv"))
cell_line_name_translator <- read.csv("data/viability/oncolines/oncolines_cell_line_translator.csv", header = TRUE)

# PREPROCESSING ----
process_gdsc(
    main_csv   = file.path(DATA_PATH, "gdsc/drug_data/raw_data_GDSC_007_the.university.of.edinburgh.asier.unciti-broceta_21Feb19_2131.csv"),
    growth_csv = file.path(DATA_PATH, "gdsc/drug_data/day1_data_GDSC_007_21Feb19_2131.csv"),
    out_path   = file.path(OUTPATH, "preprocessing", "gdsc.csv"),
    info       = metadata
)
process_oncolines(
    xlsx_path    = file.path(DATA_PATH, "oncolines/v1/Raw data_21OL897.xlsx"),
    out_path     = file.path(DATA_PATH, "preprocessing", "oncolines.csv"),
    clt          = cell_line_name_translator,
    info         = metadata,
    dataset_name = "oncolines",
    skip         = 4
)
process_oncolines(
    xlsx_path    = file.path(DATA_PATH, "oncolines/v2/Raw data_21OL897.xlsx"),
    out_path     = file.path(DATA_PATH, "preprocessing", "oncolines_v2.csv"),
    clt          = cell_line_name_translator,
    info         = metadata,
    dataset_name = "oncolines_v2",
    skip         = 4
)
process_temps(
    main_csv = file.path(DATA_PATH, "temps/nxp900_viability.csv"),
    out_path = file.path(OUTPATH, "preprocessing", "temps.csv"),
    info = metadata
)
process_brognard(
    main_csv = file.path(DATA_PATH, "brognard/nxp900_viability.csv"),
    out_path = file.path(OUTPATH, "preprocessing", "brognard.csv"),
    info = metadata
)

# DOSE-RESPONSE MODELLING ----
SCREENS <- list(
    gdsc = list(
        path         = file.path(OUT_PATH, "gdsc"),
        combine_name = "GDSC",
        responses    = c("IC50", "GI50", "GR50"),
        top_dose     = 10,
        bottom_dose  = 0.01
    ),
    oncolines = list(
        path         = file.path(OUT_PATH, "oncolines"),
        combine_name = "ONCO",
        responses    = c("IC50", "GI50", "GR50"),
        top_dose     = 31.6,
        bottom_dose  = 0.00316
    ),
    oncolines_v2 = list(
        path         = file.path(OUT_PATH, "oncolines_v2"),
        combine_name = "ONCO",
        responses    = c("IC50", "GI50", "GR50"),
        top_dose     = 31.6,
        bottom_dose  = 0.00316
    ),
    carolin = list(
        path         = file.path(OUT_PATH, "carolin"),
        combine_name = "CAROLIN",
        responses    = "GI50",
        # carolin has two dose ranges (top/bottom dose vary by cell line) —
        # one fit per cell line, AOC / DSS computed per-group from data
        grouped      = TRUE
    ),
    brognard = list(
        path         = file.path(OUT_PATH, "brognard"),
        combine_name = "BROGNARD",
        responses    = "IC50",
        top_dose     = 31.6,
        bottom_dose  = 0.00316
    )
)

# To do work, list screen keys here. REFIT/PLOT must be subsets of RUN.
RUN <- c("gdsc", "oncolines", "oncolines_v2", "carolin", "brognard") # e.g. c("gdsc", "oncolines") — assemble + write drc_metrics.csv
REFIT <- c("gdsc", "oncolines", "oncolines_v2", "carolin", "brognard") # subset of RUN — refit DRCs from scratch (else load existing rds)
PLOT <- c("gdsc", "oncolines", "oncolines_v2", "carolin", "brognard") # subset of RUN — regenerate diagnostic DRC plot grid
COMBINE <- TRUE # rebuild outputs/preprocessing/viability/viability_raw.csv

stopifnot(all(RUN %in% names(SCREENS)))
stopifnot(all(REFIT %in% RUN))
stopifnot(all(PLOT %in% RUN))

for (screen in RUN) {
    process_screen(
        SCREENS[[screen]],
        refit    = screen %in% REFIT,
        plot_drc = screen %in% PLOT
    )
}

if (COMBINE) {
    combine_drcfits(
        screens  = SCREENS,
        out_path = file.path(OUT_BASE, "viability_raw.csv")
    )
}

# POSTPROCESSING ----
postprocess_viability(
    in_path  = file.path(OUT_BASE, "viability_raw.csv"),
    out_path = file.path(OUT_BASE, "viability.csv")
)
