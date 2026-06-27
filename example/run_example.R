
# ============================================================
# run_example.R
# Demonstrates how to use the estimate_semen_age() function
# 
# Author: Chao Xiao @HUST
# License: MIT
# ============================================================

# Set working directory to the repository root
# setwd("/path/to/semen-age-snapshot-12cpg")

# Source the estimation function
source("R/estimate_age.R")

# --- View model information ---
print_model_info(model_dir = "models/")

# --- Load example data ---
dat <- read.csv("example/example_input.csv", row.names = 1)
cat("\nInput data (first 3 rows):\n")
print(head(dat, 3))

# --- Estimate with SVR (recommended) ---
cat("\n--- SVR Estimations ---\n")
results_svr <- estimate_semen_age(dat, model_dir = "models/", method = "SVR")
print(results_svr)

# --- Estimate with MLR ---
cat("\n--- MLR Estimations ---\n")
results_mlr <- estimate_semen_age(dat, model_dir = "models/", method = "MLR")
print(results_mlr)

# --- Estimate with MQR + 80% estimation intervals ---
cat("\n--- MQR Estimations with 80% PI ---\n")
results_mqr <- estimate_semen_age(dat, model_dir = "models/", method = "MQR",
                                 estimation_interval = TRUE)
print(results_mqr)

# --- Using your own data ---
# Prepare a data.frame with columns named exactly:
#   cg19998819, cg11262154, cg12277678, cg04123357,
#   cg18037145, cg25187042, cg20602007, cg21843517,
#   cg13872326, chr1.19339447, chr13.93039685, chr19.18610678
#
# Each value = methylation ratio (0 to 1) from SNaPshot assay
# my_data <- read.csv("my_semen_methylation_data.csv", row.names = 1)
# my_results <- estimate_semen_age(my_data, model_dir = "models/", method = "SVR")
