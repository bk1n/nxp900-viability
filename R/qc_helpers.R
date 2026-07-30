qc_calc_cv <- function(df, well_type = "NC1") {
    df %>%
        dplyr::filter(WELL_TYPE == well_type) %>%
        group_by(CELL_LINE_NAME, PLATE_ID) %>%
        summarise(
            mean = mean(INTENSITY),
            sd = sd(INTENSITY),
            .groups = "drop"
        ) %>%
        mutate(CV = sd / mean) %>%
        dplyr::select(CELL_LINE_NAME, PLATE_ID, CV)
}

qc_calc_ncr <- function(df, nc0_well = "NC0", nc1_well = "NC1") {
    df %>%
        dplyr::filter(WELL_TYPE %in% c(nc0_well, nc1_well)) %>%
        group_by(CELL_LINE_NAME, PLATE_ID, WELL_TYPE) %>%
        summarise(mean = mean(INTENSITY), .groups = "drop") %>%
        pivot_wider(names_from = "WELL_TYPE", values_from = "mean") %>%
        mutate(NCR = .data[[nc0_well]] / .data[[nc1_well]]) %>%
        dplyr::select(CELL_LINE_NAME, PLATE_ID, NCR)
}

qc_calc_z <- function(df, nc1_well = "NC1", b_well = "B", pc1_well = "PC1D1S", pc2_well = "PC2D1S") {
    keep <- c(nc1_well, b_well, pc1_well, pc2_well)
    z <- df %>%
        dplyr::filter(WELL_TYPE %in% keep) %>%
        group_by(CELL_LINE_NAME, PLATE_ID, WELL_TYPE) %>%
        summarise(mean = mean(INTENSITY), sd = sd(INTENSITY), .groups = "drop") %>%
        pivot_wider(names_from = "WELL_TYPE", values_from = c("mean", "sd"))

    z %>%
        mutate(
            Z_B = 1 - 3 * (.data[[paste0("sd_", b_well)]] + .data[[paste0("sd_", nc1_well)]]) / (.data[[paste0("mean_", nc1_well)]] - .data[[paste0("mean_", b_well)]]),
            Z_PC1 = 1 - 3 * (.data[[paste0("sd_", pc1_well)]] + .data[[paste0("sd_", nc1_well)]]) / (.data[[paste0("mean_", nc1_well)]] - .data[[paste0("mean_", pc1_well)]]),
            Z_PC2 = 1 - 3 * (.data[[paste0("sd_", pc2_well)]] + .data[[paste0("sd_", nc1_well)]]) / (.data[[paste0("mean_", nc1_well)]] - .data[[paste0("mean_", pc2_well)]])
        ) %>%
        dplyr::select(CELL_LINE_NAME, PLATE_ID, Z_B, Z_PC1, Z_PC2)
}

qc_calc_z_single <- function(df, high_well = "NC0", low_well = "B", out_col = "Z_B") {
    z <- df %>%
        dplyr::filter(WELL_TYPE %in% c(high_well, low_well)) %>%
        group_by(CELL_LINE_NAME, PLATE_ID, WELL_TYPE) %>%
        summarise(mean = mean(INTENSITY), sd = sd(INTENSITY), .groups = "drop") %>%
        pivot_wider(names_from = "WELL_TYPE", values_from = c("mean", "sd"))

    z[[out_col]] <- 1 - 3 * (z[[paste0("sd_", low_well)]] + z[[paste0("sd_", high_well)]]) / (z[[paste0("mean_", high_well)]] - z[[paste0("mean_", low_well)]])
    z %>% dplyr::select(CELL_LINE_NAME, PLATE_ID, dplyr::all_of(out_col))
}

qc_append_standard_columns <- function(df) {
    for (col in c("CV", "NCR", "Z_B", "Z_PC1", "Z_PC2")) {
        if (!col %in% colnames(df)) df[[col]] <- NA_real_
    }
    df %>% dplyr::select(-tidyselect::any_of(c("cv", "nc_ratio", "z_B_gi50", "cv_gi50")), everything())
}
