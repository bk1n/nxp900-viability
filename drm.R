# =============================================================================
# NXP900 viability: raw plate data -> dose-response metrics -> viability.csv
# =============================================================================
#
# Driver script for the whole viability pipeline. Run it top-to-bottom with
# `Rscript drm.R` (or source it interactively) from the repo root. All the real
# work lives in R/ and is loaded by devtools::load_all(".").
#
# ---- Pipeline -------------------------------------------------------------
#
#   1. PREPROCESSING   process_<screen>()  R/preprocessing.R
#        Reads each vendor's raw plate export, runs the screen's QC, maps cell
#        line names to DepMapID via Model.csv, and computes the response
#        transforms (R/response_transforms.R). Writes, into <OUT_PATH>/<screen>/:
#          drc_ready.csv   long per-dose data, one row per cell line x conc
#          qc_*.csv        per-plate QC tables (only for screens with controls)
#        Every drc_ready.csv is validated by check_drm_ready() before being
#        written, so a schema violation fails here rather than during fitting.
#
#   2. DOSE-RESPONSE   process_screen()    R/drm.R
#        Fits a 3-parameter log-logistic curve (drc::LL.3u, upper fixed at 100)
#        per cell line per response, then derives IC50/GI50/GR50, AOC, DSS1-3,
#        RMSE, AIC and flat/increasing flags. Writes, into <OUT_PATH>/<screen>/:
#          <RESPONSE>_models.qs2   serialised fits (reloaded when refit = FALSE)
#          drc_metrics.csv         long: DepMapID, dataset, halfmax, type, value
#          <RESPONSE>_drc.png      diagnostic curve grid (only when plot_drc)
#
#   3. COMBINE         combine_drcfits()   R/drm.R
#        Row-binds the per-screen drc_metrics.csv into <OUT_PATH>/viability_raw.csv
#        and errors on any duplicated (DepMapID, dataset, halfmax, type) key.
#
#   4. POSTPROCESSING  postprocess_viability()  R/postprocessing.R
#        Filters to the metrics of interest, censors bad fits to the top dose,
#        and writes the analysis-ready <OUT_PATH>/viability.csv.
#
# ---- Screens --------------------------------------------------------------
#
# Each entry of SCREENS describes one screen and how to fit it:
#   path         - directory holding the step-1 output and where the DRC
#                  artefacts land. Must already contain drc_ready.csv.
#   responses    - subset of c("IC50", "GI50", "GR50") to fit. A response is
#                  only available if drc_ready.csv has the matching
#                  TRT_INTENSITY_<RESPONSE> column, which depends on whether
#                  the screen supplied day-0 counts (GI50/GR50 need them).
#   top_dose,    - dose range (uM) used as the integration limits for AOC and
#   bottom_dose    DSS. Ignored when grouped = TRUE. Keep these matched to the
#                  screen's actual measured dose range, otherwise AOC/DSS
#                  integrate extrapolated curve outside the data. Datasets also
#                  listed in postprocessing.R::VIABILITY_TOP_DOSE must agree
#                  with their top_dose here.
#   combine_name - dataset label written into drc_metrics.csv. Screens sharing
#                  a label (oncolines + oncolines_v2 -> "ONCO") are row-bound
#                  into one dataset, so their cell lines must not overlap.
#   grouped      - TRUE when cell lines within the screen were dosed over
#                  different ranges. Fitting is still one curve per cell line,
#                  but AOC/DSS integration limits are taken per dose-range
#                  group from the data instead of from top_dose/bottom_dose.
#
# ---- Toggles --------------------------------------------------------------
#
# RUN / REFIT / PLOT are character vectors of SCREENS keys, and REFIT and PLOT
# must be subsets of RUN:
#   RUN    assemble drc_metrics.csv for these screens
#   REFIT  refit the curves from scratch; omitted screens reload their
#          <RESPONSE>_models.qs2 instead. Fitting a large screen takes hours,
#          so leave a screen out of REFIT once its models are on disk.
#   PLOT   regenerate the diagnostic <RESPONSE>_drc.png grid (slow, 30x30in)
#
# A screen in RUN but not REFIT will fail if its .qs2 files do not exist yet.
# Step 3 only reads drc_metrics.csv, so every screen in RUN must have been run at
# least once before it can be combined.
# =============================================================================

devtools::load_all(".")

DATA_PATH <- "data"
OUT_PATH <- "out"

# LOAD METADATA
# Model.csv is the DepMap model table; it supplies the CELL_LINE_NAME -> DepMapID
# mapping used by every finalise step. Oncolines ships non-CCLE cell line names,
# so it needs an extra translator table on top.
# DepMap names the id column ModelID; the pipeline refers to it as DepMapID
# throughout, so rename it once here rather than at every join.
metadata <- readr::read_csv(file.path(DATA_PATH, "Model.csv")) %>%
    dplyr::rename(DepMapID = ModelID)
cell_line_name_translator <- read.csv(file.path(DATA_PATH, "oncolines/oncolines_cell_line_translator.csv"), header = TRUE)

# PREPROCESSING ----
# out_path is a *directory* per screen — each process_* writes drc_ready.csv
# into it, and these directories are the same ones SCREENS$<screen>$path points
# at below. Comment out a screen once its drc_ready.csv is up to date.
process_gdsc(
    main_csv   = file.path(DATA_PATH, "gdsc/drug_data/raw_data_GDSC_007_the.university.of.edinburgh.asier.unciti-broceta_21Feb19_2131.csv"),
    growth_csv = file.path(DATA_PATH, "gdsc/drug_data/day1_data_GDSC_007_21Feb19_2131.csv"),
    out_path   = file.path(OUT_PATH, "gdsc"),
    info       = metadata
)
process_oncolines(
    xlsx_path    = file.path(DATA_PATH, "oncolines/v1/Raw data_21OL897.xlsx"),
    out_path     = file.path(OUT_PATH, "oncolines"),
    clt          = cell_line_name_translator,
    info         = metadata,
    dataset_name = "oncolines",
    skip         = 4
)
# skip = 0, not 4: oncolines_raw_v2.xlsx has already been trimmed to the header
# row, whereas Raw data_21OL897.xlsx still carries four rows of preamble above it.
# The column layout is identical between the two, so the same reader handles both.
process_oncolines(
    xlsx_path    = file.path(DATA_PATH, "oncolines/v1/oncolines_raw_v2.xlsx"),
    out_path     = file.path(OUT_PATH, "oncolines_v2"),
    clt          = cell_line_name_translator,
    info         = metadata,
    dataset_name = "oncolines_v2",
    skip         = 0
)
process_temps(
    csv_path = file.path(DATA_PATH, "temps/nxp900_viability.csv"),
    out_path = file.path(OUT_PATH, "temps"),
    info     = metadata
)
process_brognard(
    csv_path = file.path(DATA_PATH, "brognard/nxp900_viability.csv"),
    out_path = file.path(OUT_PATH, "brognard"),
    info     = metadata
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
    temps = list(
        path         = file.path(OUT_PATH, "temps"),
        combine_name = "TEMPS",
        responses    = "GI50",
        # temps has two dose ranges (top/bottom dose vary by cell line) —
        # one fit per cell line, AOC / DSS computed per-group from data
        grouped      = TRUE
    ),
    brognard = list(
        path         = file.path(OUT_PATH, "brognard"),
        combine_name = "BROGNARD",
        responses    = "IC50",
        top_dose     = 10,
        bottom_dose  = 0.0251
    )
)

# To do work, list screen keys here. REFIT/PLOT must be subsets of RUN.
RUN <- c("gdsc", "oncolines", "oncolines_v2", "temps", "brognard")
REFIT <- c("gdsc", "oncolines", "oncolines_v2", "temps", "brognard")
PLOT <- c("gdsc", "oncolines", "oncolines_v2", "temps", "brognard")

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

combine_drcfits(
    screens  = SCREENS[RUN],
    out_path = file.path(OUT_PATH, "viability_raw.csv")
)

# POSTPROCESSING ----
postprocess_viability(
    in_path  = file.path(OUT_PATH, "viability_raw.csv"),
    out_path = file.path(OUT_PATH, "viability.csv")
)
