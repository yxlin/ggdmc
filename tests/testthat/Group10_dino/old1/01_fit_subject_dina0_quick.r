#!/usr/bin/env Rscript
# Quick testing version - tests only a few small datasets for development/debugging
# q(save = "no")
cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║       Quick Test - DINA Model (Small Datasets Only)            ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Setup ====================================================================
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
pkg_ok <- sapply(pkg, require, character.only = TRUE)

if (!all(pkg_ok)) {
  stop("Failed to load required packages: ", names(pkg_ok)[!pkg_ok])
}

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"

# Test only small datasets for quick feedback
# sample_sizes <- c(100, 500, 1000) # Modify this to test different sizes
sample_sizes <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000)


# Create results directory
results_dir <- file.path(home_dir, "Group10_cdm_subjects/results_quick")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Reduced sampling for faster testing
SAMPLING_CONFIG <- list(
  nmc = 50, # Fewer iterations for quick test
  nchain = NULL,
  migration_prob = 0.05,
  thin_stage1 = 4, # Less thinning
  thin_stage2 = 2,
  is_pblocked = FALSE,
  seed = 9032,
  ncore = 3
)

cat("Quick test configuration:\n")
cat("  Sample sizes:", paste(sample_sizes, collapse = ", "), "\n")
cat("  Iterations:", SAMPLING_CONFIG$nmc, "\n")
cat("  Chains:", SAMPLING_CONFIG$nchain, "\n\n")

# Storage for results
all_results <- list()
# n <- sample_sizes[1]

# Main Testing Loop ========================================================
for (n in sample_sizes) {
  cat("────────────────────────────────────────────────────────────────\n")
  cat(sprintf("Testing N = %d\n", n))
  cat("────────────────────────────────────────────────────────────────\n")

  # Construct file paths
  data_file <- sprintf("subject_dina0_N%06d.rda", n)
  data_path <- file.path(home_dir, "Group9_gen_cdm/data", data_file)

  if (!file.exists(data_path)) {
    cat("⚠ Data file not found, skipping\n\n")
    next
  }

  cat("Loading:", data_path, "\n")
  load(data_path)

  # Stage 1: Initial sampling
  cat("Stage 1: Initial sampling...\n")
  fits0 <- StartSampling_subject(
    sub_dmis[[1]], sub_priors,
    ncore = SAMPLING_CONFIG$ncore,
    nmc = SAMPLING_CONFIG$nmc,
    nchain = SAMPLING_CONFIG$nchain,
    sub_migration_prob = SAMPLING_CONFIG$migration_prob,
    thin = SAMPLING_CONFIG$thin_stage1,
    is_pblocked = SAMPLING_CONFIG$is_pblocked,
    seed = SAMPLING_CONFIG$seed
  )


  # Stage 2: Restart sampling
  cat("Stage 2: Restart sampling...\n")
  fits1 <- ggdmc:::RestartSampling_subject(
    fits0,
    sub_migration_prob = 0.00,
    thin = SAMPLING_CONFIG$thin_stage2,
    is_pblocked = SAMPLING_CONFIG$is_pblocked,
    seed = SAMPLING_CONFIG$seed
  )

  fits <- fits1

  # Quick diagnostics
  fit <- RebuildPosterior(fits)
  hat_overall <- gelman(fit)
  cat(sprintf(
    "\nOverall MPSRF = %.4f %s\n",
    hat_overall$mpsrf,
    ifelse(hat_overall$mpsrf < 1.1, "✓", "⚠")
  ))

  # Quick estimation
  # true_vector
  # est_theta <- ggdmc::summary(fit)
  # est_theta

  options(digits = 4)
  est_theta <- ggdmc::compare(fit, ps = true_vector)
  cat("\nParameter estimates:\n")
  # print(est_theta)
  # str(est_theta)
  # row.names(est_theta)

  # Store minimal results
  bias <- est_theta["Median-True", ]

  all_results[[as.character(n)]] <- list(
    n = n,
    mpsrf = hat_overall$mpsrf,
    bias = bias,
    mae = mean(abs(bias))
  )

  cat(sprintf("\nMean Absolute Error: %.4f\n\n", mean(abs(bias))))

  # Clean up
  rm(fits0, fits1, fits, fit, sub_dmis, sub_priors)
  gc()
}

# Quick Summary ============================================================
cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                     QUICK TEST SUMMARY                         ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

if (length(all_results) > 0) {
  summary_df <- data.frame(
    MPSRF = sapply(all_results, function(x) x$mpsrf),
    Converged = sapply(all_results, function(x) x$mpsrf < 1.1),
    Mean_MAE = sapply(all_results, function(x) x$mae)
  )

  print(summary_df)

  cat("\n")
  if (all(summary_df$Converged)) {
    cat("✓ All chains converged\n")
  } else {
    cat("⚠ Some chains did not converge\n")
  }

  if (nrow(summary_df) > 1) {
    cat(sprintf(
      "\nMAE trend: %.4f → %.4f (N: %d → %d)\n",
      summary_df$Mean_MAE[1],
      summary_df$Mean_MAE[nrow(summary_df)],
      summary_df$N[1],
      summary_df$N[nrow(summary_df)]
    ))
  }
} else {
  cat("No results to summarize\n")
}

cat("\n✓ Quick test completed\n\n")
