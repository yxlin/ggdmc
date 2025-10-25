#!/usr/bin/env Rscript
# Benchmark script comparing original vs fast DE optimizer
# This version runs twice as required for correct results

rm(list = ls())
library(ggdmc)

# Check if rbenchmark is available, if not use system.time()
use_rbenchmark <- requireNamespace("rbenchmark", quietly = TRUE)
if (!use_rbenchmark) {
    cat("Note: rbenchmark package not found. Using system.time() instead.\n")
    cat("      Install rbenchmark for more robust timing: install.packages('rbenchmark')\n\n")
}

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
wkdir <- file.path(home_dir, "Group13")
data_file <- file.path(wkdir, "simple_benchmark.rda")
setwd(wkdir)

cat("\n============================================================\n")
cat("Performance Benchmark: Original vs Fast DE Optimizer\n")
cat("============================================================\n")
cat("Working directory: ", getwd(), "\n")

# Load test data
load(data_file)
cat("Loaded: config, dmi, samples\n")
cat("Number of samples: ", samples@nmc, "\n")
cat("Number of chains: ", length(samples@nchain), "\n\n")

# ============================================================
# Define wrapper functions that run twice (as required)
# ============================================================

run_original_twice <- function() {
    fit0 <- ggdmc:::run_subject(config, dmi, samples)
    fit1 <- ggdmc:::run_subject(config, dmi, samples = fit0)
    return(fit1)
}

run_fast_twice <- function() {
    fit0 <- ggdmc:::run_subject_fast(config, dmi, samples)
    fit1 <- ggdmc:::run_subject_fast(config, dmi, samples = fit0$posterior)
    return(fit1)
}

# ============================================================
# Run Benchmark
# ============================================================

if (use_rbenchmark) {
    cat("Using rbenchmark for performance comparison...\n")
    cat("Running each version 3 times (replications = 3)\n")
    cat("Each replication includes 2 sequential runs as required\n\n")

    results <- rbenchmark::benchmark(
        original = run_original_twice(),
        fast = run_fast_twice(),
        replications = 3,
        columns = c("test", "replications", "elapsed", "relative", "user.self", "sys.self")
    )

    cat("============================================================\n")
    cat("Benchmark Results (rbenchmark)\n")
    cat("============================================================\n")
    print(results)
    cat("\n")

    # Calculate speedup
    elapsed_original <- results[results$test == "original", "elapsed"]
    elapsed_fast <- results[results$test == "fast", "elapsed"]
    speedup <- elapsed_original / elapsed_fast
    time_saved <- elapsed_original - elapsed_fast

    cat("------------------------------------------------------------\n")
    cat("Summary:\n")
    cat("  Original version: ", sprintf("%.3f seconds", elapsed_original), "\n")
    cat("  Fast version:     ", sprintf("%.3f seconds", elapsed_fast), "\n")
    cat("  Speedup:          ", sprintf("%.3fx", speedup), "\n")
    cat("  Time saved:       ", sprintf(
        "%.3f seconds (%.1f%%)",
        time_saved,
        100 * (1 - elapsed_fast / elapsed_original)
    ), "\n")
    cat("------------------------------------------------------------\n\n")
} else {
    # Fallback to system.time() with multiple runs
    cat("Using system.time() for performance comparison...\n")
    cat("Running each version 3 times for stability\n\n")

    cat("Running ORIGINAL version (3 times)...\n")
    times_original <- numeric(3)
    for (i in 1:3) {
        cat("  Run", i, "... ")
        timing <- system.time({
            fit_original <- run_original_twice()
        })
        times_original[i] <- timing["elapsed"]
        cat(sprintf("%.3f seconds\n", timing["elapsed"]))
    }

    cat("\nRunning FAST version (3 times)...\n")
    times_fast <- numeric(3)
    for (i in 1:3) {
        cat("  Run", i, "... ")
        timing <- system.time({
            fit_fast <- run_fast_twice()
        })
        times_fast[i] <- timing["elapsed"]
        cat(sprintf("%.3f seconds\n", timing["elapsed"]))
    }

    # Calculate statistics
    mean_original <- mean(times_original)
    sd_original <- sd(times_original)
    mean_fast <- mean(times_fast)
    sd_fast <- sd(times_fast)
    speedup <- mean_original / mean_fast
    time_saved <- mean_original - mean_fast

    cat("\n============================================================\n")
    cat("Benchmark Results (system.time)\n")
    cat("============================================================\n")
    cat(sprintf("Original: %.3f ± %.3f seconds (mean ± sd)\n", mean_original, sd_original))
    cat(sprintf("Fast:     %.3f ± %.3f seconds (mean ± sd)\n", mean_fast, sd_fast))
    cat("------------------------------------------------------------\n")
    cat("Summary:\n")
    cat("  Speedup:          ", sprintf("%.3fx", speedup), "\n")
    cat("  Time saved:       ", sprintf(
        "%.3f seconds (%.1f%%)",
        time_saved,
        100 * (1 - mean_fast / mean_original)
    ), "\n")
    cat("------------------------------------------------------------\n\n")
}

# ============================================================
# Get detailed C++ timing from last fast run
# ============================================================
cat("============================================================\n")
cat("C++ Internal Timing Breakdown (Fast Version)\n")
cat("============================================================\n")

# Run once more to get timing details
fit_fast_final <- ggdmc:::run_subject_fast(config, dmi, samples)
fit_fast_final <- ggdmc:::run_subject_fast(config, dmi, samples = fit_fast_final$posterior)
timing <- fit_fast_final$timing

cat(sprintf("Total C++ time:      %.3f seconds (100.0%%)\n", timing$total_time))
cat(sprintf(
    "  Crossover:         %.3f seconds (%.1f%%)\n",
    timing$crossover_time,
    100 * timing$crossover_time / timing$total_time
))
cat(sprintf(
    "  Migration:         %.3f seconds (%.1f%%)\n",
    timing$migration_time,
    100 * timing$migration_time / timing$total_time
))
cat(sprintf(
    "  Likelihood:        %.3f seconds (%.1f%%)\n",
    timing$likelihood_time,
    100 * timing$likelihood_time / timing$total_time
))

other_time <- timing$total_time - (timing$crossover_time + timing$migration_time)
cat(sprintf(
    "  Other (MH, etc):   %.3f seconds (%.1f%%)\n",
    other_time,
    100 * other_time / timing$total_time
))

cat("\nOperation counts:\n")
cat(sprintf("  Crossovers:        %d\n", timing$n_crossover))
cat(sprintf("  Migrations:        %d\n", timing$n_migration))

cat("\n------------------------------------------------------------\n")
cat("Interpretation:\n")
cat("------------------------------------------------------------\n")
if (timing$likelihood_time / timing$total_time > 0.7) {
    cat("• Likelihood computation dominates (>70% of time)\n")
    cat("  → Next optimization should focus on likelihood functions\n")
} else if (timing$crossover_time / timing$total_time > 0.5) {
    cat("• Crossover operations are significant (>50% of time)\n")
    cat("  → Chain shuffling optimization is important\n")
}

cat("\n============================================================\n")
cat("Numerical Accuracy Check\n")
cat("============================================================\n")
cat("Both versions should produce statistically similar results.\n")
cat("(Check MCMC diagnostics manually if needed)\n\n")

cat("To verify numerical accuracy:\n")
cat("  fit_old  <- run_original_twice()\n")
cat("  fit_fast <- run_fast_twice()$posterior\n")
cat("  ggdmc::compare(fit_old, ps = p_vector)\n")
cat("  ggdmc::compare(fit_fast, ps = p_vector)\n")

cat("\n============================================================\n")
cat("Benchmark Complete!\n")
cat("============================================================\n")
