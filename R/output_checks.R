# Helpers to validate output schema for dose-response fitting

check_drm_ready <- function(df, dataset_name = "dataset", require_all_trt = FALSE) {
    required_initial <- c("DepMapID", "CONC", "INTENSITY", "PLATE_ID", "PLATE_TYPE", "WELL_TYPE", "CELL_LINE_NAME")

    if (!is.data.frame(df)) stop(dataset_name, ": must be a data.frame")

    missing <- setdiff(required_initial, colnames(df))
    if (length(missing) > 0) stop(dataset_name, ": missing required initial columns: ", paste(missing, collapse = ", "))

    # Check types
    if (!is.numeric(df$CONC) && !all(is.na(df$CONC))) stop(dataset_name, ": CONC must be numeric or all NA")
    if (!is.numeric(df$INTENSITY)) stop(dataset_name, ": INTENSITY must be numeric")

    # Check PLATE_TYPE values if present
    if ("PLATE_TYPE" %in% colnames(df)) {
        ok_types <- c("MAIN", "GROWTH")
        present <- unique(na.omit(as.character(df$PLATE_TYPE)))
        bad <- setdiff(present, ok_types)
        if (length(bad) > 0) stop(dataset_name, ": unexpected PLATE_TYPE values: ", paste(bad, collapse = ", "))
    }

    # Check TRT columns - require at least one unless require_all_trt
    trt_cols <- c("TRT_INTENSITY_IC50", "TRT_INTENSITY_GI50", "TRT_INTENSITY_GR50")
    present_trt <- trt_cols[trt_cols %in% colnames(df)]
    if (require_all_trt) {
        missing_trt <- setdiff(trt_cols, present_trt)
        if (length(missing_trt) > 0) stop(dataset_name, ": missing TRT columns: ", paste(missing_trt, collapse = ", "))
    } else {
        if (length(present_trt) == 0) stop(dataset_name, ": no TRT_INTENSITY_* columns found; at least one required for fitting")
    }

    # For any present TRT column, check rows with non-NA TRT have non-NA CONC and non-NA CELL_LINE_NAME and PLATE_ID
    for (col in present_trt) {
        bad_rows <- which(!is.na(df[[col]]) & (is.na(df$CONC) | is.na(df$CELL_LINE_NAME) | is.na(df$PLATE_ID)))
        if (length(bad_rows) > 0) stop(dataset_name, ": rows with non-NA ", col, " but missing CONC/CELL_LINE_NAME/PLATE_ID. Example row: ", bad_rows[1])
    }

    # Check no empty strings in CELL_LINE_NAME
    if (any(df$CELL_LINE_NAME == "" , na.rm = TRUE)) stop(dataset_name, ": CELL_LINE_NAME contains empty string(s)")

    message("[check] ", dataset_name, ": drc-ready validation passed. rows=", nrow(df), " cols=", ncol(df))
    invisible(TRUE)
}
