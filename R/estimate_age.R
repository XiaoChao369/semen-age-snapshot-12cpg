
# ============================================================
# estimate_age.R
# Estimates chronological age from semen DNA methylation data
# using the 12-CpG SNaPshot multiplex assay models.
# 
# Author: Chao Xiao
# License: MIT
# ============================================================

#' Estimate age from semen DNA methylation (12-CpG SNaPshot panel)
#'
#' @param methylation_data A data.frame with 12 columns corresponding to
#'   methylation ratios (0–1) at the 12 CpG markers. Column names must
#'   match the R-safe marker identifiers (see Details).
#' @param model_dir Path to the directory containing model .rds files.
#'   Default: "models/" relative to script location.
#' @param method One of "SVR" (recommended), "MLR", or "MQR".
#'   SVR achieved the lowest MAE (2.82 years) in independent validation.
#' @param estimation_interval Logical. If TRUE and method = "MQR", returns
#'   80% estimation intervals (tau = 0.10 to 0.90). Default: FALSE.
#'
#' @return A data.frame with columns:
#'   - sample_id: row names from input
#'   - estimated_age: point estimate (years)
#'   - lower_80PI: lower bound of 80% PI (MQR only, if requested)
#'   - upper_80PI: upper bound of 80% PI (MQR only, if requested)
#'   - method: model used
#'
#' @details
#' Required column names (R-safe format):
#'   cg19998819, cg11262154, cg12277678, cg04123357,
#'   cg18037145, cg25187042, cg20602007, cg21843517,
#'   cg13872326, chr1.19339447, chr13.93039685, chr19.18610678
#'
#' All values should be methylation ratios between 0 and 1.
#' Input DNA should be bisulfite-converted sperm genomic DNA
#' quantified by the SNaPshot single-base extension protocol.
#'
#' @examples
#' \dontrun{
#' dat <- read.csv("example/example_input.csv", row.names = 1)
#' results <- estimate_semen_age(dat, model_dir = "models/", method = "SVR")
#' print(results)
#' }
#'
#' @export
estimate_semen_age <- function(methylation_data,
                              model_dir = "models/",
                              method = c("SVR", "MLR", "MQR"),
                              estimation_interval = FALSE) {
  
  method <- match.arg(method)
  
  # --- Required markers ---
  required_cols <- c("cg19998819", "cg11262154", "cg12277678", "cg04123357",
                     "cg18037145", "cg25187042", "cg20602007", "cg21843517",
                     "cg13872326", "chr1.19339447", "chr13.93039685",
                     "chr19.18610678")
  
  # --- Input validation ---
  missing_cols <- setdiff(required_cols, colnames(methylation_data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         "\n  See ?estimate_semen_age for required column names.")
  }
  
  # Subset and reorder to match model expectations
  input_data <- methylation_data[, required_cols, drop = FALSE]
  
  # Range check
  vals <- unlist(input_data)
  if (any(vals < 0 | vals > 1, na.rm = TRUE)) {
    warning("Some methylation values are outside [0, 1]. ",
            "Ensure values are ratios, not percentages.")
  }
  if (any(is.na(vals))) {
    stop("Input contains NA values. All 12 markers must have valid measurements.")
  }
  
  # --- Load model ---
  model_file <- switch(method,
                       "MLR" = file.path(model_dir, "mlr_full.rds"),
                       "SVR" = file.path(model_dir, "svr_full.rds"),
                       "MQR" = file.path(model_dir, "mqr_full.rds")
  )
  
  if (!file.exists(model_file)) {
    stop("Model file not found: ", model_file,
         "\n  Ensure model_dir points to the correct directory.")
  }
  
  # --- estimate ---
  if (method == "MLR") {
    model <- readRDS(model_file)
    preds <- predict(model, newdata = input_data)
    result <- data.frame(
      sample_id     = rownames(methylation_data),
      estimated_age = round(as.numeric(preds), 2),
      method        = "MLR",
      stringsAsFactors = FALSE
    )
    
  } else if (method == "SVR") {
    if (!requireNamespace("e1071", quietly = TRUE)) {
      stop("Package 'e1071' is required for SVR estimation. ",
           "Install with: install.packages('e1071')")
    }
    model <- readRDS(model_file)
    preds <- predict(model, newdata = input_data)
    result <- data.frame(
      sample_id     = rownames(methylation_data),
      estimated_age = round(as.numeric(preds), 2),
      method        = "SVR",
      stringsAsFactors = FALSE
    )
    
  } else if (method == "MQR") {
    if (!requireNamespace("quantreg", quietly = TRUE)) {
      stop("Package 'quantreg' is required for MQR estimation. ",
           "Install with: install.packages('quantreg')")
    }
    mqr_bundle <- readRDS(model_file)
    pred_050 <- predict(mqr_bundle$tau_050, newdata = input_data)
    
    result <- data.frame(
      sample_id     = rownames(methylation_data),
      estimated_age = round(as.numeric(pred_050), 2),
      method        = "MQR",
      stringsAsFactors = FALSE
    )
    
    if (estimation_interval) {
      pred_010 <- predict(mqr_bundle$tau_010, newdata = input_data)
      pred_090 <- predict(mqr_bundle$tau_090, newdata = input_data)
      result$lower_80PI <- round(as.numeric(pred_010), 2)
      result$upper_80PI <- round(as.numeric(pred_090), 2)
    }
  }
  
  rownames(result) <- NULL
  return(result)
}


#' Print model metadata and performance summary
#'
#' @param model_dir Path to directory containing model_metadata.rds
#' @export
print_model_info <- function(model_dir = "models/") {
  meta_file <- file.path(model_dir, "model_metadata.rds")
  if (!file.exists(meta_file)) {
    stop("Metadata file not found: ", meta_file)
  }
  meta <- readRDS(meta_file)
  
  cat("============================================================\n")
  cat("  ", meta$assay_name, "\n")
  cat("============================================================\n")
  cat("  Tissue:      ", meta$target_tissue, "\n")
  cat("  Population:  ", meta$population, "\n")
  cat("  Training N:  ", meta$n_training, "\n")
  cat("  Age range:   ", meta$age_range[1], "-", meta$age_range[2], "years\n")
  cat("  Markers:     ", meta$n_markers, "CpG sites\n")
  cat("------------------------------------------------------------\n")
  cat("  Performance (independent test set, n = 76):\n")
  cat("    SVR MAE:   ", meta$performance$svr_test_MAE, "years\n")
  cat("    MLR MAE:   ", meta$performance$mlr_test_MAE, "years\n")
  cat("    MQR MAE:   ", meta$performance$mqr_test_MAE, "years\n")
  cat("  Cross-validated (100x10-fold, n = 260):\n")
  cat("    SVR MAE:   ", meta$performance$svr_cv_MAE, "years\n")
  cat("    MLR MAE:   ", meta$performance$mlr_cv_MAE, "years\n")
  cat("    MQR MAE:   ", meta$performance$mqr_cv_MAE, "years\n")
  cat("------------------------------------------------------------\n")
  cat("  SVR hyperparameters:\n")
  cat("    Kernel:    ", meta$svr_hyperparams$kernel, "\n")
  cat("    Cost:      ", meta$svr_hyperparams$cost, "\n")
  cat("    Gamma:     ", meta$svr_hyperparams$gamma, "\n")
  cat("    Epsilon:   ", meta$svr_hyperparams$epsilon, "\n")
  cat("============================================================\n")
  cat("  Reference: ", meta$reference, "\n")
  cat("============================================================\n")
  
  invisible(meta)
}
