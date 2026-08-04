library(tidyverse)
library(drc)

pbapply::pboptions(type = "timer")

.rmse <- function(model) {
    if (!is.null(model)) {
        sqrt(mean(residuals(model)^2))
    } else {
        NA
    }
}

# computes area above the curve by integrating across dose-range
# see Hafner, 2016 on why AOC is better than AUC
.aoc <- function(model, bottom_dose, top_dose) {
    if (!is.null(model)) {
        conc_min <- bottom_dose
        conc_max <- top_dose

        # integrate on the log10 dose axis: a linear integral is dominated by
        # the top decade of the dose range and barely separates curves
        f <- function(logconc) {
            nd <- data.frame(CONC = 10^logconc)
            100 - predict(model, newdata = nd)
        }
        tryCatch(
            {
                return(
                    integrate(f, lower = log10(conc_min), upper = log10(conc_max))$value
                )
            },
            error = function(e) {
                print("integrate failed, returning NA")
                return(NA)
            }
        )
    } else {
        return(
            NA
        )
    }
}

# computes DSS
# t: minimum activity level, proportion between 0.1 and 1, corresponds to Amin in paper
# normalisation_method: one of 'GI50', 'IC50', 'GR50', for GI50 and GR50 the range to calculate t percentage will be set to 100, -100 whilst for IC50 the range will be set to 100, 0
# see Yadav et al., 2014
.dss <- function(model, normalisation_method, t = 0.1, bottom_dose, top_dose) {
    if (is.null(model)) {
        return(NA)
    }

    if (!normalisation_method %in% c("GI50", "GR50", "IC50")) {
        stop("normalisation_method must be one of GI50, GR50, or IC50")
    } else if (normalisation_method == "GI50") {
        # DSS is undefined for GI50, as DSS requires scaling to fixed area and GI50 can go to -Inf
        return(
            data.frame(
                "DSS1" = NA,
                "DSS2" = NA,
                "DSS3" = NA
            )
        )
    } else {
        if (normalisation_method == "GR50") {
            t. <- 100 - (200 * t)
            # activity is 100 - response, so a -100 response is 200% activity
            amax <- 200
        } else {
            t. <- 100 - (100 * t)
            amax <- 100
        }
        me <- 100 - .extract_coefs(model)[, "c"]
    }

    conc_min <- bottom_dose
    conc_max <- top_dose

    nd <- data.frame(
        CONC = 10^(seq(log10(conc_min), log10(conc_max), length.out = 1000))
    )

    preds <- predict(model, nd)

    # calculate conc at which t. is reached
    if (all(preds > t.)) {
        return(
            data.frame(
                "DSS1" = 0,
                "DSS2" = 0,
                "DSS3" = 0
            )
        )
    } else {
        t_conc <- nd$CONC[min(which(preds <= t.))]
    }

    # integrate from t_conc to top_dose, on the log10 dose axis (Yadav, 2014)
    f <- function(logconc) {
        nd <- data.frame(CONC = 10^logconc)
        100 - predict(model, newdata = nd)
    }

    aoc_total <- tryCatch(
        {
            integrate(f, lower = log10(t_conc), upper = log10(conc_max))$value
        },
        error = function(e) {
            print("Error in integrate, returning NA")
            return(NA)
        }
    )

    if (is.na(aoc_total)) {
        dss1 <- dss2 <- dss3 <- NA
    } else {
        # all areas on the log10 dose axis, to match aoc_total
        lconc_min <- log10(conc_min)
        lconc_max <- log10(conc_max)
        lt_conc <- log10(t_conc)

        aoc_below_t <- aoc_total - (amax * t) * (lconc_max - lt_conc)
        # amplitude is (max activity - threshold activity), not the threshold
        max_area <- (amax - (amax * t)) * (lconc_max - lconc_min)

        # DSS1
        dss1 <- aoc_below_t / max_area

        # DSS2
        # scale me between zero and 1 for GR50
        # DSS2 cannot be defined for ME of 0 or flat curves (these will be set to zero anyway)
        dss2_scaler <- me / amax
        dss2 <- if (is.na(dss2_scaler) || dss2_scaler == 0 || dss2_scaler >= 1) NA else dss1 / -log10(dss2_scaler)

        # DSS3
        dss3 <- dss2 * (lconc_max - lt_conc) / (lconc_max - lconc_min)
    }

    data.frame(
        "DSS1" = dss1,
        "DSS2" = dss2,
        "DSS3" = dss3
    )
}

# get model AIC
.aic <- function(model) {
    if (!is.null(model)) {
        AIC(model)
    } else {
        NA
    }
}

.extract_coefs <- function(model) {
    if (!is.null(model)) {
        coefs <- coef(model)
        names(coefs) <- gsub(":\\(Intercept\\)", "", names(coefs))
        coefs_se <- summary(model)$coefficients[, "Std. Error"]
        names(coefs_se) <- paste0(gsub(":\\(Intercept\\)", "", names(coefs_se)), "_se")

        df <- cbind(
            t(data.frame(coefs)),
            t(data.frame(coefs_se))
        )
        return(df)
    } else {
        return(NA)
    }
}

# returns curve stats at specified response levels with columns as levels
.extract_stats <- function(model, levels = c(10, 50, 90), type = "relative") {
    if (!is.null(model)) {
        stats <- drc::ED(model, respLev = levels, type = type, display = F)
        stats <- stats %>%
            as.data.frame() %>%
            rownames_to_column("level") %>%
            dplyr::select(-`Std. Error`) %>%
            pivot_wider(names_from = level, values_from = Estimate) %>%
            as.data.frame()
        colnames(stats) <- paste0(type, str_split_i(colnames(stats), ":", 3))
        return(stats)
    } else {
        return(NA)
    }
}

# given data.frame, checks if curves are increasing
.increasing <- function(model) {
    if (is.null(model)) {
        return(NA)
    }
    if (sign(coef(model)[1]) < 0) TRUE else FALSE
}

# set parameters for increasing or flat curves
.fix_response <- function(drm_df, top_dose, bottom_dose) {
    drm_df %>%
        # extrapolation
        # mutate(across(c(matches("10|50|90")), function(x) if_else(x > top_dose, top_dose, x))) %>%
        # mutate(across(c(matches("10|50|90")), function(x) if_else(x < bottom_dose, bottom_dose, x))) %>%
        # increasing curves
        mutate(across(c(matches("10|50|90")), function(x) if_else(incr, top_dose, x))) %>%
        mutate(across(c(matches("^aoc$|^aoct$")), function(x) if_else(incr, 0, x))) %>%
        mutate(across(c(matches("^c$")), function(x) if_else(incr, 100, x))) %>%
        # flat curves
        mutate(across(c(matches("10|50|90")), function(x) if_else(isflat, top_dose, x))) %>%
        mutate(across(c(matches("^aoc$|^aoct$")), function(x) if_else(isflat, 0, x))) %>%
        mutate(across(c(matches("^c$")), function(x) if_else(isflat, 100, x)))
}

# data should have at least three columns: CELL_LINE_NAME, TRT_INTENSITY and CONC
drm_fit <- function(data, lwr_limit = -Inf) {
    # empty vars
    model <- NULL
    flat <- FALSE

    # checks
    stopifnot(all(c("CELL_LINE_NAME", "TRT_INTENSITY", "CONC") %in% colnames(data)))
    # cat(sprintf("Fitting models for %s\n", unique(data$CELL_LINE_NAME)))

    lower_e <- min(data$CONC)
    upper_e <- max(data$CONC) * 10^2

    upper_d <- 100

    # drmc
    control <- drmc(maxIt = 1000, relTol = 1e-6, noMessage = TRUE, errorm = FALSE, rmNA = TRUE)
    # fit LL4
    model <- tryCatch(
        drc::drm(
            TRT_INTENSITY ~ CONC,
            data = data,
            fct = drc::LL.3u(upper = 100),
            start = c(2, 50, median(data$CONC, na.rm = T)),
            lowerl = c(-Inf, lwr_limit, lower_e),
            upperl = c(Inf, upper_d, upper_e),
            control = control
        ),
        error = function(e) {
            cat(sprintf("Fit failed: %s\n", conditionMessage(e)))
            NULL
        }
    )

    # with errorm = FALSE drc returns a bare list (convergence = FALSE) instead
    # of a "drc" object on failure; normalise it so the is.null() guards fire
    if (!inherits(model, "drc")) {
        return(list(
            data = data,
            model = NULL,
            flat = flat
        ))
    }

    # F-test for the significance of the sigmoidal fit (Hafner, 2016)
    Npara <- 3 # 4PL w/ d fixed
    Npara_flat <- 1 # flat
    RSS2 <- sum(stats::residuals(model)^2, na.rm = TRUE)
    RSS1 <- sum((data$TRT_INTENSITY - mean(data$TRT_INTENSITY, na.rm = TRUE))^2,
        na.rm = TRUE
    )
    df1 <- (Npara - Npara_flat)
    df2 <- (length(na.omit(data$TRT_INTENSITY)) - Npara + 1)
    f_value <- ((RSS1 - RSS2) / df1) / (RSS2 / df2)
    f_pval <- stats::pf(f_value, df1, df2, lower.tail = FALSE)

    if (f_pval >= 0.05) {
        flat <- TRUE
    }

    return(list(
        data = data,
        model = model,
        flat = flat
    ))
}

# Applies dose-response modelling to each cell line in standard data
# data should be in standard format (see `test_viability_data.R`)
# Args:
# - name: name of dataset, columns will be prefixed with name e.g. GI50_d
drm_dataframe <- function(data, models, name, top_dose, bottom_dose) {
    g <- data %>%
        as_tibble() %>%
        dplyr::select(DepMapID, CELL_LINE_NAME) %>%
        distinct() %>%
        arrange(match(CELL_LINE_NAME, names(models)))

    # models
    drm <- lapply(models, function(x) x$model)

    # anova
    isflat <- lapply(models, function(x) x$flat)
    isflat <- do.call(rbind, isflat)

    # stats
    stats <- lapply(drm, .extract_stats, type = "relative")
    stats <- do.call(rbind, stats)
    statsabs <- lapply(drm, .extract_stats, type = "absolute")
    statsabs <- do.call(rbind, statsabs)
    stats <- cbind(stats, statsabs)
    rownames(stats) <- names(drm)

    # coefficients
    model_coefs <- lapply(drm, .extract_coefs)
    coefs <- do.call(rbind, model_coefs)
    rownames(coefs) <- names(model_coefs)

    # aoc
    aoc <- unlist(lapply(drm, .aoc, top_dose = top_dose, bottom_dose = bottom_dose))
    aoc <- data.frame(aoc = aoc)
    aoct <- unlist(lapply(drm, .aoc, top_dose = 0.5, bottom_dose = bottom_dose))
    aoct <- data.frame(aoct = aoct)

    # dss
    dss <- lapply(drm, .dss, normalisation_method = name, t = 0.1, bottom_dose = bottom_dose, top_dose = top_dose)
    dss <- do.call(rbind, dss)

    # rmse
    rmse <- unlist(lapply(drm, .rmse))
    rmse <- data.frame(rmse = rmse)

    # AIC
    aic <- unlist(lapply(drm, .aic))
    aic <- data.frame(aic = aic)

    # increasing
    incr <- unlist(lapply(drm, .increasing))

    stopifnot(
        g$CELL_LINE_NAME == rownames(coefs) &
            rownames(coefs) == rownames(stats) &
            rownames(coefs) == rownames(aoc) &
            rownames(coefs) == rownames(rmse) &
            rownames(coefs) == rownames(aic)
    )
    g <- cbind(g, coefs, stats, aoc, aoct, dss, rmse, aic, isflat, incr)

    # DON'T FIX RESPONSES YET
    g <- .fix_response(g, top_dose, bottom_dose)

    # rename df
    g <- g %>%
        rename_with(.fn = function(x) paste0(name, "_", x), .cols = c(-DepMapID, -CELL_LINE_NAME))

    return(g)
}

# plotting ----
plot_drm <- function(drm_df, models, name, outpath) {
    # otherwise, plot DRC
    colnames(drm_df) <- gsub("IC50_|GI50_|GR50_", "", colnames(drm_df))

    stopifnot(rownames(drm_df) == names(models))

    ylim <- ifelse(name == "IC50", list(c(0, 100)), list(c(-100, 100)))[[1]]

    models <- models[order(drm_df$incr, drm_df$isflat)]

    png(file.path(outpath, paste0(name, "_drc.png")), width = 30, height = 30, units = "in", res = 300)
    par(
        mfrow = c(ceiling(sqrt(length(models))), ceiling(sqrt(length(models)))),
        mar = c(0, 0, 0, 0)
    )
    for (n in names(models)) {
        d <- models[[n]]$data
        m <- models[[n]]$model
        fl <- models[[n]]$flat
        if (is.null(m)) {
            plot(d$CONC, y = d$TRT_INTENSITY, log = "x", xlab = "", xaxt = "n", ylab = "", yaxt = "n", ylim = ylim)
            text(1, 0, "FTF", col = "red", cex = 1)
        } else {
            plot(m, type = "all", xlab = "", ylab = "", xt = 0, yt = -500, ylim = ylim)
            if (fl) {
                text(1, 50, "Flat", col = "red", cex = 1)
            }
            if (drm_df$incr[rownames(drm_df) == n]) {
                text(1, 20, "Incr.", col = "red", cex = 1)
            }
        }
    }
    dev.off()
}

# orchestration ----
# Standard lower-bound for the lower asymptote (c) per response transform.
# IC50 floor at 0 (viability never goes negative); GI50 / GR50 can be
# negative when treatment kills below day-0 cell count.
.default_lwr_limit <- function(response) {
    switch(response,
        IC50 = 0,
        GI50 = -Inf,
        GR50 = -100,
        stop("unknown response: ", response)
    )
}

# Fit one response transform across all cell lines in a screen.
# Expects screen_data to contain TRT_INTENSITY_<response>, CONC, CELL_LINE_NAME.
.fit_response <- function(screen_data, response) {
    trt_col <- paste0("TRT_INTENSITY_", response)
    if (!trt_col %in% colnames(screen_data)) {
        stop("response ", response, " requires column ", trt_col)
    }
    d <- screen_data %>%
        dplyr::select(CELL_LINE_NAME, TRT_INTENSITY = dplyr::all_of(trt_col), CONC)
    d <- split(d, d$CELL_LINE_NAME)
    pbapply::pblapply(d, drm_fit, lwr_limit = .default_lwr_limit(response))
}

# Inner-join per-response drm dataframes on DepMapID; CELL_LINE_NAME kept once.
.combine_responses <- function(dfs) {
    g <- dfs[[1]]
    for (i in seq_along(dfs)[-1]) {
        g <- inner_join(g, dfs[[i]] %>% dplyr::select(-CELL_LINE_NAME), by = "DepMapID")
    }
    g %>% relocate(CELL_LINE_NAME, .after = DepMapID)
}

# Pivot a wide drm_dataframe output (DepMapID, CELL_LINE_NAME, <RESPONSE>_<param>...)
# into the canonical long format: DepMapID, dataset, halfmax, type, value.
# Booleans (isflat, incr) are coerced to integer 0/1 so the value column stays
# numeric after pivot_longer.
.to_long <- function(wide_df, dataset) {
    response_pattern <- "^(IC50|GI50|GR50)_"
    value_cols <- grep(response_pattern, colnames(wide_df), value = TRUE)
    stopifnot(length(value_cols) > 0)

    wide_df %>%
        dplyr::select(DepMapID, dplyr::all_of(value_cols)) %>%
        mutate(across(where(is.logical), as.integer)) %>%
        tidyr::pivot_longer(
            cols = dplyr::all_of(value_cols),
            names_to = c("halfmax", "type"),
            names_pattern = "^(IC50|GI50|GR50)_(.+)$",
            values_to = "value"
        ) %>%
        mutate(dataset = dataset) %>%
        dplyr::select(DepMapID, dataset, halfmax, type, value)
}

# Standard per-screen DRC workflow: load drc_ready.csv, fit (or load) DRC
# models per response, build the wide drc_metrics dataframe, optionally plot,
# then write drc_metrics.csv into cfg$path. Filenames are fixed (see README):
#   in:  drc_ready.csv
#   out: <RESPONSE>_models.qs2 (one per response in cfg$responses)
#        drc_metrics.csv       (long: DepMapID, dataset, halfmax, type, value)
#        <RESPONSE>_drc.png    (if plot_drc)
#
# cfg is one entry of drm.R's SCREENS list — see that file for the contract.
process_screen <- function(cfg, refit = FALSE, plot_drc = FALSE) {
    cat("=== ", basename(cfg$path), " ===\n", sep = "")
    screen_data <- read.csv(file.path(cfg$path, "drc_ready.csv"))

    if (isTRUE(cfg$grouped)) {
        .process_screen_grouped(cfg, screen_data, refit, plot_drc)
        return(invisible(NULL))
    }

    models <- list()
    for (response in cfg$responses) {
        qs_path <- file.path(cfg$path, sprintf("%s_models.qs2", response))
        if (refit) {
            cat("Fitting", response, "models\n")
            models[[response]] <- .fit_response(screen_data, response)
            cat("Saving", response, "models ->", qs_path, "\n")
            qs2::qs_save(models[[response]], qs_path)
        } else {
            models[[response]] <- qs2::qs_read(qs_path)
        }
    }

    dfs <- lapply(cfg$responses, function(response) {
        df <- drm_dataframe(
            screen_data, models[[response]], response,
            top_dose = cfg$top_dose, bottom_dose = cfg$bottom_dose
        )
        if (plot_drc) plot_drm(df, models[[response]], response, outpath = cfg$path)
        df
    })
    names(dfs) <- cfg$responses

    stopifnot(length(unique(lapply(dfs, function(d) d$DepMapID))) == 1)
    g <- .combine_responses(dfs)
    long <- .to_long(g, dataset = cfg$combine_name)

    out_csv <- file.path(cfg$path, "drc_metrics.csv")
    write.csv(long, out_csv, row.names = FALSE)
    cat(basename(cfg$path), "DONE ->", out_csv, "\n")
}

# Grouped variant: screens (temps) where different cell lines were dosed at
# different ranges. One shared fit per cell line, but AOC / DSS post-processing
# uses each group's (min, max) CONC.
.process_screen_grouped <- function(cfg, screen_data, refit, plot_drc) {
    stopifnot(length(cfg$responses) == 1)
    response <- cfg$responses
    qs_path <- file.path(cfg$path, sprintf("%s_models.qs2", response))

    if (refit) {
        cat("Fitting", response, "models (grouped screen)\n")
        models_all <- .fit_response(screen_data, response)
        cat("Saving", response, "models ->", qs_path, "\n")
        qs2::qs_save(models_all, qs_path)
    } else {
        models_all <- qs2::qs_read(qs_path)
    }

    groups <- screen_data %>%
        group_by(CELL_LINE_NAME) %>%
        summarise(min = min(CONC), max = max(CONC), .groups = "drop") %>%
        group_by(min, max) %>%
        reframe(cell = list(CELL_LINE_NAME))

    per_group <- lapply(seq_len(nrow(groups)), function(i) {
        cells <- groups$cell[[i]]
        gmodels <- models_all[match(cells, names(models_all))]
        gdata <- screen_data %>%
            dplyr::filter(CELL_LINE_NAME %in% cells) %>%
            mutate(CELL_LINE_NAME = factor(CELL_LINE_NAME, levels = cells)) %>%
            arrange(CELL_LINE_NAME) %>%
            mutate(CELL_LINE_NAME = as.character(CELL_LINE_NAME))
        list(
            df = drm_dataframe(gdata, gmodels, response,
                top_dose = groups$max[i], bottom_dose = groups$min[i]
            ),
            models = gmodels
        )
    })

    models <- do.call(c, lapply(per_group, `[[`, "models"))
    df <- do.call(rbind, lapply(per_group, `[[`, "df"))
    stopifnot(rownames(df) == names(models))

    if (plot_drc) plot_drm(df, models, response, outpath = cfg$path)
    long <- .to_long(df, dataset = cfg$combine_name)

    out_csv <- file.path(cfg$path, "drc_metrics.csv")
    write.csv(long, out_csv, row.names = FALSE)
    cat(basename(cfg$path), "DONE ->", out_csv, "\n")
}

# Combine per-screen drc_metrics.csv files into one long viability table.
# Each per-screen csv is already long (DepMapID, dataset, halfmax, type, value)
# and carries its own dataset label, so combining is just a row-bind. Screens
# sharing combine_name (e.g. oncolines + oncolines_v2 -> ONCO) will contribute
# rows under the same dataset value — duplicate (DepMapID, dataset, halfmax,
# type) keys are errored on, matching the old "no duplicated DepMapID" guard.
combine_drcfits <- function(screens, out_path) {
    parts <- lapply(screens, function(s) read.csv(file.path(s$path, "drc_metrics.csv")))
    combined <- do.call(rbind, parts)

    key_cols <- c("DepMapID", "dataset", "halfmax", "type")
    if (any(duplicated(combined[, key_cols]))) {
        stop("duplicated (", paste(key_cols, collapse = ", "), ") rows in combined output")
    }

    if (!dir.exists(dirname(out_path))) {
        dir.create(dirname(out_path), recursive = TRUE)
    }
    write.csv(combined, out_path, row.names = FALSE)
    cat("Full DRCs combined ->", out_path, "\n")
}
