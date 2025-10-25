#!/usr/bin/env Rscript
# Generate benchmark data for HB performance testing
# This script sets up the config, dmis, and samples needed for benchmarking

rm(list = ls())
library(ggdmc)

# ============================================================================
# Configuration
# ============================================================================

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
test_dir <- file.path(home_dir, "Group14")
data_file <- file.path(home_dir, "Group0_gen_data/data/lba_data0.rda")
save_file <- file.path(test_dir, "hb_benchmark.rda")

cat("============================================================\n")
cat("HB Benchmark Data Generator\n")
cat("============================================================\n\n")

# Check source data exists
if (!file.exists(data_file)) {
    cat("ERROR: Source data file not found:", data_file, "\n")
    cat("Please ensure the test data exists.\n")
    quit(status = 1)
}

setwd(test_dir)
cat("Working directory:", getwd(), "\n")
cat("Loading source data...\n")
load(data_file)

# ============================================================================
# Set up HB model configuration
# ============================================================================

cat("\nConfiguring HB model...\n")

# DE parameters
pop_migration_prob <- 0.00
sub_migration_prob <- 0.00
gamma_precursor <- 2.38
rp <- 0.001
is_hblocked <- FALSE
is_pblocked <- FALSE

# Get model parameters
nparameter <- pop_priors@nparameter
pnames <- pop_priors@pnames
nchain <- ggdmc:::.get_nchain(pop_priors, NULL)

# Debug flags
pop_debug <- FALSE
sub_debug <- FALSE

# Create DE input
de_input <- setDEInput(
    pop_migration_prob = pop_migration_prob,
    sub_migration_prob = sub_migration_prob,
    gamma_precursor = gamma_precursor,
    rp = rp,
    is_hblocked = is_hblocked,
    is_pblocked = is_pblocked,
    nparameter = as.integer(nparameter),
    nchain = as.integer(nchain),
    pop_debug = pop_debug,
    sub_debug = sub_debug
)

# Other settings
ncore <- 1L
seed <- 123

# Create configuration
config_list <- set_configs(
    prior = pop_priors,
    theta_input = pop_theta_input,
    de_input = de_input,
    ncore = ncore,
    seed = seed
)

# Extract objects needed for benchmark
config <- config_list[[1]]
dmis <- pop_dmis
samples <- pop_samples

# ============================================================================
# Validate and summarize
# ============================================================================

cat("\nValidating configuration...\n")

n_subjects <- length(dmis)
n_iterations <- samples$phi@nmc
n_chains <- samples$phi@nchain
n_params <- length(samples$subject_theta[[1]]@pnames)

cat("\nBenchmark Configuration Summary:\n")
cat("------------------------------------------------------------\n")
cat(sprintf("Number of subjects:       %d\n", n_subjects))
cat(sprintf("MCMC iterations:          %d\n", n_iterations))
cat(sprintf("Number of chains:         %d\n", n_chains))
cat(sprintf("Parameters per subject:   %d\n", n_params))
cat(sprintf("Hyper-parameters:         %d\n", nparameter))

# ============================================================================
# Save benchmark data
# ============================================================================

cat("\nSaving benchmark data...\n")
save(config, dmis, samples, file = save_file)

if (file.exists(save_file)) {
    cat(sprintf("✓ Saved: %s\n", save_file))
} else {
    cat("✗ ERROR: Failed to save\n")
    quit(status = 1)
}

cat("\n============================================================\n")
cat("Ready! Run: Rscript hb_benchmark.R\n")
cat("============================================================\n")
