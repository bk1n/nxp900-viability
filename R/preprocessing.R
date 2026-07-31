# Constants ----
## Oncolines ----
# Dose well code -> log10(concentration in M). Converted to uM downstream.
ONCOLINES_DOSES <- c(
    "d1" = -4.5, "d2" = -5.0, "d3" = -5.5, "d4" = -6.0, "d5" = -6.5,
    "d6" = -7.0, "d7" = -7.5, "d8" = -8.0, "d9" = -8.5
)

# Oncolines-style cell line names whose translated form will not match CCLE 24Q2.
# Accepted as missing — anything else that drops out is a regression.
ONCOLINES_KNOWN_EXCL <- c(
    "FAC", "Jurkat E6.1", "KPPC", "LS 174T", "MM.1R", "SB1", "UWB1.289+BRCA1", "MSK921", "HCC1588"
)

## GDSC ----
# The single GDSC drug we keep (eCF506).
GDSC_DRUG_ID <- 2437

# Z-factor metrics computed on the MAIN plate (NC1 vs each high-signal well).
# Drives the per-metric qc_zfactor_*.png filenames and the qc_metrics.csv columns.
GDSC_ZFACTORS <- c(B = "Z_B", PC1 = "Z_PC1", PC2 = "Z_PC2")

# Shared ----
# Some screens arrive already summarised to per-dose responses, with no raw control wells.
# Fill in the plate columns the fitting code expects so they match the per-well screens.
standardise_plate_cols <- function(d) {
    if (!"INTENSITY" %in% colnames(d)) {
        trt_col <- intersect(c("TRT_INTENSITY_GI50", "TRT_INTENSITY_IC50"), colnames(d))
        if (length(trt_col) == 0) {
            stop("standardise_plate_cols: no INTENSITY and no TRT_INTENSITY_GI50/IC50 to derive it from")
        }
        d$INTENSITY <- d[[trt_col[1]]]
    }
    if (!"PLATE_ID" %in% colnames(d)) {
        d$PLATE_ID <- ave(d$CELL_LINE_NAME, d$CELL_LINE_NAME, FUN = seq_along)
    }
    d$PLATE_TYPE <- "MAIN"
    if (!"WELL_TYPE" %in% colnames(d)) {
        d$WELL_TYPE <- if ("TAG" %in% colnames(d)) as.character(d$TAG) else "TRT"
    }
    d
}

# Join a pre-summarised screen to CCLE ids on the stripped cell line name, then put the
# columns in drc-ready order. Cell lines absent from CCLE are dropped by the inner join.
finalise_summarised_screen <- function(d, info) {
    ccle <- info %>% dplyr::select(DepMapID, StrippedCellLineName)

    d %>%
        inner_join(ccle, by = c("CELL_LINE_NAME" = "StrippedCellLineName")) %>%
        standardise_plate_cols() %>%
        relocate(DepMapID, CONC, INTENSITY, PLATE_ID, PLATE_TYPE, WELL_TYPE, CELL_LINE_NAME)
}

# GDSC ----

gdsc_read <- function(main_csv, growth_csv) {
    d <- data.table::fread(main_csv)
    d1 <- data.table::fread(growth_csv)

    dropped <- unique(d$CELL_LINE_NAME[!d$COSMIC_ID %in% d1$COSMIC_ID])
    if (length(dropped) > 0) {
        cat("Cell lines absent from day-1 data, dropped from MAIN plate:\n")
        print(dropped)
    }

    main_plate <- d %>%
        mutate(
            PLATE_ID = SCAN_ID,
            PLATE_TYPE = "MAIN",
            WELL_TYPE = gsub("-", "", TAG),
            DepMapID = NA_character_
        ) %>%
        dplyr::select(DepMapID, COSMIC_ID, CELL_LINE_NAME, PLATE_ID, PLATE_TYPE, WELL_TYPE, CONC, INTENSITY, DRUG_ID)

    growth_plate <- d1 %>%
        mutate(
            PLATE_ID = SCAN_ID,
            PLATE_TYPE = "GROWTH",
            WELL_TYPE = gsub("-", "", TAG),
            CONC = NA_real_,
            DepMapID = NA_character_
        ) %>%
        dplyr::select(DepMapID, COSMIC_ID, CELL_LINE_NAME, PLATE_ID, PLATE_TYPE, WELL_TYPE, CONC, INTENSITY)

    list(main = main_plate, growth = growth_plate)
}

# CV, NCR, Z-factor QC for the MAIN plate; writes 5 histograms with thresholds drawn.
gdsc_qc_main <- function(main_plate, out_path) {
    plate_meta <- main_plate %>% dplyr::distinct(COSMIC_ID, CELL_LINE_NAME, PLATE_ID)
    cv <- qc_calc_cv(main_plate, well_type = "NC1")
    ncr <- qc_calc_ncr(main_plate, nc0_well = "NC0", nc1_well = "NC1")
    z_factor <- qc_calc_z(main_plate, nc1_well = "NC1", b_well = "B", pc1_well = "PC1D1S", pc2_well = "PC2D1S")

    png(file.path(out_path, "qc_cv.png"), width = 8, height = 8, units = "in", res = 600)
    hist(cv$CV, main = "NC1 CV", breaks = 30, xlab = "CV")
    abline(v = 0.18, col = "black", lty = "dashed")
    dev.off()

    png(file.path(out_path, "qc_ncr.png"), width = 8, height = 8, units = "in", res = 600)
    hist(ncr$NCR, main = "NC0/NC1 ratio", breaks = 30, xlab = "NC0/NC1 ratio")
    abline(v = 0.8, col = "darkblue", lty = "dashed")
    abline(v = 1.2, col = "darkred", lty = "dashed")
    dev.off()

    for (label in names(GDSC_ZFACTORS)) {
        col <- GDSC_ZFACTORS[[label]]
        png(file.path(out_path, sprintf("qc_zfactor_%s.png", tolower(label))),
            width = 8, height = 8, units = "in", res = 600
        )
        hist(z_factor[[col]], main = paste("Z-factor", label), breaks = 30, xlab = "Z-factor")
        abline(v = 0.3, col = "black", lty = "dashed")
        dev.off()
    }

    qc_result <- plate_meta %>%
        left_join(cv, by = c("CELL_LINE_NAME", "PLATE_ID")) %>%
        left_join(ncr, by = c("CELL_LINE_NAME", "PLATE_ID")) %>%
        left_join(z_factor, by = c("CELL_LINE_NAME", "PLATE_ID")) %>%
        mutate(qcPASS = CV < 0.18 & (NCR >= 0.8 & NCR <= 1.2) & (Z_B > 0.3 & Z_PC1 > 0.3 & Z_PC2 > 0.3)) %>%
        qc_append_standard_columns() %>%
        relocate(COSMIC_ID, CELL_LINE_NAME, PLATE_ID)

    write.csv(qc_result, file.path(out_path, "qc_metrics.csv"), row.names = FALSE)
    qc_result
}

# CV + NC0-vs-B Z-factor QC for the day-1 GROWTH plate.
gdsc_qc_growth <- function(growth_plate, out_path) {
    growth_meta <- growth_plate %>% dplyr::distinct(COSMIC_ID, CELL_LINE_NAME, PLATE_ID)
    cv_day1 <- qc_calc_cv(growth_plate, well_type = "NC0")
    z_day1 <- qc_calc_z_single(growth_plate, high_well = "NC0", low_well = "B", out_col = "Z_B")

    qc_result <- growth_meta %>%
        left_join(cv_day1, by = c("CELL_LINE_NAME", "PLATE_ID")) %>%
        left_join(z_day1, by = c("CELL_LINE_NAME", "PLATE_ID")) %>%
        mutate(qcPASS = CV < 0.18 & Z_B > 0.3) %>%
        qc_append_standard_columns() %>%
        relocate(COSMIC_ID, CELL_LINE_NAME, PLATE_ID)

    write.csv(qc_result, file.path(out_path, "qc_metrics_growth.csv"), row.names = FALSE)
    qc_result
}

# Drop failed plates, restrict to one drug, blank-subtract per (COSMIC_ID, PLATE_ID),
# keep treatment wells only (TAG starts with "L").
gdsc_summarise_trt <- function(main_plate, drug_id, rm_scans) {
    blanks <- main_plate %>%
        dplyr::filter(WELL_TYPE == "B") %>%
        group_by(COSMIC_ID, PLATE_ID) %>%
        summarise(B_MAIN = mean(INTENSITY), .groups = "drop")

    main_plate %>%
        dplyr::filter(!PLATE_ID %in% rm_scans, DRUG_ID == drug_id) %>%
        dplyr::select(COSMIC_ID, PLATE_ID, CELL_LINE_NAME, WELL_TYPE, CONC, INTENSITY) %>%
        inner_join(blanks, by = c("COSMIC_ID", "PLATE_ID")) %>%
        mutate(INTENSITY = INTENSITY - B_MAIN) %>%
        dplyr::select(-B_MAIN) %>%
        dplyr::filter(grepl("^L", WELL_TYPE))
}

# Per-plate NC0 / NC1 baselines from MAIN (blank-subtracted).
gdsc_summarise_ctrl <- function(main_plate, rm_scans) {
    main_plate %>%
        dplyr::filter(!PLATE_ID %in% rm_scans, WELL_TYPE %in% c("NC0", "NC1", "B")) %>%
        group_by(COSMIC_ID, CELL_LINE_NAME, PLATE_ID, WELL_TYPE) %>%
        summarise(mean = mean(INTENSITY), .groups = "drop") %>%
        pivot_wider(names_from = "WELL_TYPE", values_from = "mean") %>%
        mutate(
            NC0 = NC0 - B,
            NC1 = NC1 - B
        ) %>%
        dplyr::select(-B)
}

# Per-cell-line NC0 baseline from GROWTH (blank-subtracted, averaged across plates).
gdsc_summarise_growth <- function(growth_plate, rm_scans_gi50) {
    growth_plate %>%
        dplyr::filter(!PLATE_ID %in% rm_scans_gi50, WELL_TYPE %in% c("NC0", "B")) %>%
        group_by(COSMIC_ID, CELL_LINE_NAME, PLATE_ID, WELL_TYPE) %>%
        summarise(mean = mean(INTENSITY), .groups = "drop") %>%
        pivot_wider(names_from = "WELL_TYPE", values_from = "mean") %>%
        mutate(NC0_GROWTH = NC0 - B) %>%
        dplyr::select(COSMIC_ID, NC0_GROWTH) %>%
        group_by(COSMIC_ID) %>%
        summarise(NC0_GROWTH = mean(NC0_GROWTH), .groups = "drop")
}

# Apply IC50 / GI50 / GR50 transforms; join to CCLE ids.
gdsc_finalise <- function(trt, ctrl, ctrl_day1, info) {
    trt %>%
        inner_join(ctrl, by = c("COSMIC_ID", "CELL_LINE_NAME", "PLATE_ID")) %>%
        inner_join(ctrl_day1, by = "COSMIC_ID") %>%
        mutate(
            TRT_INTENSITY_IC50 = ic50(INTENSITY, NC1),
            TRT_INTENSITY_GI50 = gi50(INTENSITY, NC1, NC0_GROWTH),
            TRT_INTENSITY_GR50 = gr(INTENSITY, NC1, NC0_GROWTH),
            PLATE_TYPE = "MAIN",
            TAG = WELL_TYPE
        ) %>%
        dplyr::inner_join(info, by = c("COSMIC_ID" = "COSMICID")) %>%
        relocate(DepMapID, CONC, INTENSITY, PLATE_ID, PLATE_TYPE, WELL_TYPE) %>%
        relocate(CELL_LINE_NAME, .after = COSMIC_ID)
}

# End-to-end: raw csvs -> drc-ready trt CSV plus QC artefacts.
process_gdsc <- function(main_csv, growth_csv, out_path, info,
                         drug_id = GDSC_DRUG_ID, dataset_name = "gdsc") {
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

    plates <- gdsc_read(main_csv, growth_csv)

    qc_main <- gdsc_qc_main(plates$main, out_path)
    qc_growth <- gdsc_qc_growth(plates$growth, out_path)

    rm_scans <- qc_main$PLATE_ID[!qc_main$qcPASS]
    rm_scans_gi50 <- qc_growth$PLATE_ID[!qc_growth$qcPASS]

    trt <- gdsc_summarise_trt(plates$main, drug_id, rm_scans)
    ctrl <- gdsc_summarise_ctrl(plates$main, rm_scans)
    ctrl_day1 <- gdsc_summarise_growth(plates$growth, rm_scans_gi50)

    trt_processed <- gdsc_finalise(trt, ctrl, ctrl_day1, info)

    check_drm_ready(trt_processed, dataset_name = dataset_name)
    write.csv(trt_processed, file.path(out_path, "drc_ready.csv"), row.names = FALSE)

    cat(dataset_name, "DONE\n")
    invisible(trt_processed)
}

# Oncolines ----
# Read the raw Oncolines xlsx and return long-format per-well plate data.
oncolines_read <- function(xlsx_path, skip = 0) {
    d <- readxl::read_xlsx(xlsx_path, skip = skip)

    names <- d$`Cell line name`
    names <- names[-1]
    names[is.na(names)] <- names[which(is.na(names)) - 1]

    df <- data.frame(
        CELL_LINE_NAME = names,
        PLATE_ID = ave(names, names, FUN = seq_along),
        NC1_main_1 = d[-1, 13, drop = TRUE],
        NC1_main_2 = d[-1, 14, drop = TRUE],
        NC0_GI50_1 = d[-1, 2, drop = TRUE],
        NC0_GI50_2 = d[-1, 3, drop = TRUE],
        d1 = d[-1, 4, drop = TRUE],
        d2 = d[-1, 5, drop = TRUE],
        d3 = d[-1, 6, drop = TRUE],
        d4 = d[-1, 7, drop = TRUE],
        d5 = d[-1, 8, drop = TRUE],
        d6 = d[-1, 9, drop = TRUE],
        d7 = d[-1, 10, drop = TRUE],
        d8 = d[-1, 11, drop = TRUE],
        d9 = d[-1, 12, drop = TRUE]
    )

    df %>%
        pivot_longer(cols = NC1_main_1:d9, names_to = "RAW_TAG", values_to = "INTENSITY") %>%
        mutate(
            WELL_TYPE = case_when(
                RAW_TAG %in% c("NC1_main_1", "NC1_main_2") ~ "NC1",
                RAW_TAG %in% c("NC0_GI50_1", "NC0_GI50_2") ~ "NC0",
                .default = RAW_TAG
            ),
            PLATE_TYPE = if_else(WELL_TYPE == "NC0", "GROWTH", "MAIN"),
            CONC = if_else(WELL_TYPE %in% names(ONCOLINES_DOSES),
                10^ONCOLINES_DOSES[WELL_TYPE] * 10^6, NA_real_
            ),
            DepMapID = NA_character_
        ) %>%
        dplyr::select(DepMapID, CELL_LINE_NAME, PLATE_ID, PLATE_TYPE, WELL_TYPE, CONC, INTENSITY)
}

# Per-plate controls + per-dose trt intensities -> trt table with IC50/GI50/GR50 transforms.
oncolines_summarise <- function(plate_df) {
    ctrl <- plate_df %>%
        dplyr::filter(WELL_TYPE %in% c("NC1", "NC0")) %>%
        group_by(CELL_LINE_NAME, PLATE_ID, WELL_TYPE) %>%
        summarise(mean = mean(INTENSITY), .groups = "drop") %>%
        pivot_wider(names_from = "WELL_TYPE", values_from = "mean") %>%
        dplyr::rename(NC0_GROWTH = NC0)

    trt <- plate_df %>%
        dplyr::filter(WELL_TYPE %in% names(ONCOLINES_DOSES)) %>%
        group_by(CELL_LINE_NAME, PLATE_ID, WELL_TYPE, CONC) %>%
        summarise(INTENSITY = mean(INTENSITY), .groups = "drop")

    trt %>%
        inner_join(ctrl, by = c("CELL_LINE_NAME", "PLATE_ID")) %>%
        mutate(
            TRT_INTENSITY_IC50 = ic50(INTENSITY, NC1),
            TRT_INTENSITY_GI50 = gi50(INTENSITY, NC1, NC0_GROWTH),
            TRT_INTENSITY_GR50 = gr(INTENSITY, NC1, NC0_GROWTH),
            PLATE_TYPE = "MAIN",
            TAG = WELL_TYPE
        )
}

# CV-based per-plate QC; writes a histogram with the 0.18 threshold marked.
oncolines_qc <- function(plate_df, out_path) {
    cv <- qc_calc_cv(plate_df, well_type = "NC1")

    png(file.path(out_path, "qc_cv.png"), width = 8, height = 8, units = "in", res = 600)
    hist(cv$CV, main = "NC1 CV", breaks = 30, xlab = "CV")
    abline(v = 0.18, col = "black", lty = "dashed")
    dev.off()

    qc_result <- cv %>%
        mutate(qcPASS = CV < 0.18) %>%
        qc_append_standard_columns()
    write.csv(qc_result, file.path(out_path, "qc_metrics.csv"), row.names = FALSE)

    qc_result
}

# Rewrite Oncolines names into DepMap form using the hand-curated translator.
make_oncolines_translator <- function(clt) {
    function(x) if_else(x %in% clt$Oncolines, clt$DepMap[match(x, clt$Oncolines)], x)
}

# Apply QC filter, validate the unmapped set against known_excl, join to CCLE ids.
oncolines_finalise <- function(trt_processed, qc_result, clt, info,
                               known_excl = ONCOLINES_KNOWN_EXCL) {
    qcr <- qc_result %>% dplyr::select(CELL_LINE_NAME, PLATE_ID, qcPASS)
    translate <- make_oncolines_translator(clt)

    qc_passed <- trt_processed %>%
        left_join(qcr, by = c("CELL_LINE_NAME", "PLATE_ID")) %>%
        dplyr::filter(qcPASS)

    # canary: any QC-passed cell line whose translated form does not appear in CCLE 24Q2
    # must already be on the accepted exclusion list. Compared in Oncolines-style names so
    # translator / CCLE release / raw xlsx going stale all surface here instead of silently
    # shrinking the output via the inner_join below.
    oncolines_seen <- unique(qc_passed$CELL_LINE_NAME)
    unmapped <- oncolines_seen[!translate(oncolines_seen) %in% info$CellLineName]
    unexpected <- setdiff(unmapped, known_excl)
    if (length(unexpected) > 0) {
        stop(
            "Unmapped Oncolines cell lines not in known_excl: ",
            paste(unexpected, collapse = ", ")
        )
    }

    qc_passed %>%
        mutate(CELL_LINE_NAME = translate(CELL_LINE_NAME)) %>%
        inner_join(info, by = c("CELL_LINE_NAME" = "CellLineName")) %>%
        relocate(DepMapID, CONC, INTENSITY, PLATE_ID, PLATE_TYPE, WELL_TYPE)
}

# End-to-end: raw xlsx -> drc-ready trt CSV plus QC artefacts.
process_oncolines <- function(xlsx_path, out_path, clt, info, dataset_name, skip = 0) {
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

    plate_df <- oncolines_read(xlsx_path, skip = skip)
    trt_processed <- oncolines_summarise(plate_df)
    qc_result <- oncolines_qc(plate_df, out_path)
    trt_processed <- oncolines_finalise(trt_processed, qc_result, clt, info)

    check_drm_ready(trt_processed, dataset_name = dataset_name)
    write.csv(trt_processed, file.path(out_path, "drc_ready.csv"), row.names = FALSE)

    cat(dataset_name, "DONE\n")
    invisible(trt_processed)
}

# Temps ----
# Read the raw Temps csv. Doses are molar, and OVCAR3 is the only name that does not
# already match its CCLE stripped form.
temps_read <- function(csv_path) {
    read.csv(csv_path) %>%
        mutate(
            CONC = CONC * 1e6, # M -> uM
            CELL_LINE_NAME = if_else(CELL_LINE_NAME == "OVCAR3", "NIHOVCAR3", CELL_LINE_NAME)
        )
}

# End-to-end: raw csv -> drc-ready trt CSV. No QC artefacts (no control wells supplied).
process_temps <- function(csv_path, out_path, info, dataset_name = "temps") {
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

    trt_processed <- temps_read(csv_path) %>%
        finalise_summarised_screen(info)

    check_drm_ready(trt_processed, dataset_name = dataset_name)
    write.csv(trt_processed, file.path(out_path, "drc_ready.csv"), row.names = FALSE)

    cat(dataset_name, "DONE\n")
    invisible(trt_processed)
}

# Brognard ----
# Read the raw Brognard csv: drop the untreated rows and the per-replicate columns,
# keeping only the replicate mean. Every cell line gets a single plate (bar a special few),
# and PLATE_ID is already supplied.
brognard_read <- function(csv_path) {
    read.csv(csv_path) %>%
        dplyr::filter(logCONC != 0) %>%
        dplyr::select(-logCONC, -c(n1, n2, n3)) %>%
        mutate(TRT_INTENSITY_IC50 = TRT_INTENSITY_IC50 * 100) # fraction -> %
}

# End-to-end: raw csv -> drc-ready trt CSV. No QC artefacts (no control wells supplied).
# MSK921 is absent from CCLE 24Q2 and is expected to drop out of the join.
process_brognard <- function(csv_path, out_path, info, dataset_name = "brognard") {
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

    trt_processed <- brognard_read(csv_path) %>%
        finalise_summarised_screen(info)

    check_drm_ready(trt_processed, dataset_name = dataset_name)
    write.csv(trt_processed, file.path(out_path, "drc_ready.csv"), row.names = FALSE)

    cat(dataset_name, "DONE\n")
    invisible(trt_processed)
}
