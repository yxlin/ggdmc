#!/usr/bin/env Rscript
# Comprehensive benchmark for Hierarchical Bayesian models
# Compares three versions:
# 1. Original run() - Baseline
# 2. run_fast() - With cached likelihood and timing
# 3. run_fast_notiming() - With cached likelihood, no timing overhead

rm(list = ls())
library(ggdmc)

use_rbenchmark <- requireNamespace("rbenchmark", quietly = TRUE)
if (!use_rbenchmark) {
    cat("Note: rbenchmark package not found. Using system.time() instead.\n")
    cat("      Install: install.packages('rbenchmark')\n\n")
}

cat("\n============================================================\n")
cat("Hierarchical Bayesian Model Performance Comparison\n")
cat("============================================================\n")

cat("This benchmark compares three versions of HB sampling:\n")
cat("1. Original run() - Recalculates hyper-likelihood every time\n")
cat("2. run_fast() - Caches hyper-likelihood + timing instrumentation\n")
cat("3. run_fast_notiming() - Caches hyper-likelihood, no timing overhead\n\n")

cat("Expected speedup depends on:\n")
cat("  • Number of subjects (more subjects = bigger speedup)\n")
cat("  • Number of hyper-level iterations\n")
cat("  • Problem complexity\n\n")

cat("------------------------------------------------------------\n")
cat("USAGE:\n")
cat("------------------------------------------------------------\n")
cat("1. Set up your HB model (config, dmis, samples)\n")
cat("2. Source this script or modify paths below\n")
cat("3. Review results and choose best version\n\n")

# ===========================================================================
# CONFIGURATION - Modify these paths for your test data
# ===========================================================================

# Example: Uncomment and modify for your test case
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
test_dir <- file.path(home_dir, "Group14")
data_file <- file.path(home_dir, "Group0_gen_data/data/lba_data0.rda")
hb_data <- file.path(home_dir, "Group14/hb_benchmark.rda")


if (!file.exists(data_file)) {
    cat("ERROR: Data file not found:", data_file, "\n")
    cat("Please modify the paths in this script to point to your HB test data.\n")
    quit(status = 1)
}

setwd(test_dir)
load(data_file)
load(hb_data)



# Check if required objects exist
if (!exists("config") || !exists("dmis") || !exists("samples")) {
    cat("\n============================================================\n")
    cat("Setup Required\n")
    cat("============================================================\n")
    cat("This script requires the following objects to be defined:\n")
    cat("  • config - Configuration object\n")
    cat("  • dmis - List of DMI objects (one per subject)\n")
    cat("  • samples - List with phi and subject_theta\\n\n")
    cat("To use this benchmark:\n")
    cat("1. Load or create your HB test data\n")
    cat("2. Modify the configuration section in this script\n")
    cat("3. Re-run the script\\n\n")
    quit(save = "no")
}

# ===========================================================================
# Validate data structure
# ===========================================================================

cat("\n============================================================\n")
cat("Data Validation\n")

n_subjects <- length(dmis)
cat(sprintf("Number of subjects: %d\n", n_subjects))

# ===========================================================================
# Define wrapper functions
# ===========================================================================


run_original <- function() {
    result <- ggdmc::run(config, dmis, samples)
    return(result)
}

run_cached_with_timing <- function() {
    result <- ggdmc::run_fast(config, dmis, samples)
    return(result)
}

run_cached_notiming <- function() {
    result <- ggdmc::run_fast_notiming(config, dmis, samples)
    return(result)
}

# ===========================================================================
# Run Benchmark
# ===========================================================================
if (use_rbenchmark) {
    cat("Using rbenchmark for performance comparison...\n")

    results <- rbenchmark::benchmark(
        original = run_original(),
        fast_with_timing = run_cached_with_timing(),
        fast_notiming = run_cached_notiming(),
        replications = 3,
        columns = c("test", "replications", "elapsed", "relative")
    )
}
results
result <- ggdmc::run_fast_notiming(config, dmis, samples)
phi <- result$phi
subject_theta <- result$subject_theta
options(digits = 2)
est_phi <- compare(phi, ps = true_vector, start = phi@nmc * 0.5)

DT <- prepare_thetas_data(subject_theta, start = phi@nmc * 0.5)
p1 <- plot_thetas(DT)
