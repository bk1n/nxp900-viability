# Post-process the raw combined DRC table from drm.R into the canonical
# viability artefact used by downstream EDA / modelling.
#
# Both input and output are long format (DepMapID, dataset, halfmax, type, value).
# Two corrections on top of `.fix_response` (which already runs inside
# drm_dataframe per screen):
#   1. Cap good-fit absolute50 past TOP_DOSE[dataset]. `.fix_response`
#      deliberately does not cap extrapolated good fits — that policy is
#      applied here, post-fit.
#   2. Override absolute50 / aoc / aoct / DSS for flat or increasing curves
#      (fit did not converge meaningfully). The absolute50 / aoc / aoct
#      overrides are redundant with `.fix_response`, but DSS is only fixed here.
#
# The output is then restricted to VIABILITY_KEEP_TYPES (the curated metric set
# the old eda/viability/postprocess.R produced via `matches(...)`).

# Per-dataset top dose (uM) for capping extrapolated absolute50. Must stay in
# sync with the corresponding `top_dose` entries in drm.R::SCREENS. Only
# datasets listed here have their good-fit absolute50 capped; others pass
# through unchanged (e.g. TEMPS, whose top dose is per-group, not a single
# screen-wide value).
VIABILITY_TOP_DOSE <- c(GDSC = 10^1, ONCO = 10^1.5, BROGNARD = 10^1.5)

# Metric types kept in the output (mirrors the column selection the old
# eda/viability/postprocess.R applied via `matches(...)`).
VIABILITY_KEEP_TYPES <- c(
    "absolute50", "c", "aoc", "aoct",
    "DSS1", "DSS2", "DSS3",
    "aic", "b", "incr", "isflat"
)

VIABILITY_ID_COLS <- c("DepMapID", "dataset", "halfmax")

# One badfit flag per (DepMapID, dataset, halfmax) by pivoting isflat / incr.
viability_badfit_flags <- function(v) {
    v %>%
        dplyr::filter(type %in% c("isflat", "incr")) %>%
        tidyr::pivot_wider(names_from = "type", values_from = "value") %>%
        dplyr::mutate(badfit = (isflat == 1) | (incr == 1)) %>%
        dplyr::select(dplyr::all_of(c(VIABILITY_ID_COLS, "badfit")))
}

# Apply the two post-fit corrections described at the top of this file.
correct_viability_fits <- function(v, top_dose = VIABILITY_TOP_DOSE) {
    v %>%
        dplyr::left_join(viability_badfit_flags(v), by = VIABILITY_ID_COLS) %>%
        dplyr::mutate(
            # 1. cap good-fit absolute50 past dataset top_dose
            value = dplyr::case_when(
                type == "absolute50" &
                    dataset %in% names(top_dose) &
                    value > top_dose[dataset] ~ top_dose[dataset],
                .default = value
            ),
            # 2. flat / increasing override
            value = dplyr::case_when(
                badfit & type == "absolute50" &
                    dataset %in% names(top_dose) ~ top_dose[dataset],
                badfit & type %in% c("aoc", "aoct", "DSS1", "DSS2", "DSS3") ~ 0,
                .default = value
            )
        ) %>%
        dplyr::select(-badfit)
}

# in:  <in_path>  viability_raw.csv from combine_drcfits()
# out: <out_path> viability.csv
postprocess_viability <- function(in_path, out_path,
                                  top_dose = VIABILITY_TOP_DOSE,
                                  keep_types = VIABILITY_KEEP_TYPES) {
    v <- read.csv(in_path)
    stopifnot(all(c(VIABILITY_ID_COLS, "type", "value") %in% colnames(v)))

    v <- v %>%
        dplyr::filter(type %in% keep_types, !is.na(value)) %>%
        correct_viability_fits(top_dose = top_dose)

    write.csv(v, out_path, row.names = FALSE)
    cat("postprocess: wrote", nrow(v), "rows ->", out_path, "\n")
    invisible(v)
}
