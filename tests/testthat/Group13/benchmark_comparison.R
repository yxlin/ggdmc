#!/usr/bin/env Rscript
# Comprehensive benchmark comparing three versions:
# 1. Original (baseline)
# 2. Fast with timing (has chrono overhead)
# 3. Fast without timing (pure optimization)

rm(list = ls())
library(ggdmc)

use_rbenchmark <- requireNamespace("rbenchmark", quietly = TRUE)
if (!use_rbenchmark) {
    cat("Note: rbenchmark package not found. Using system.time() instead.\n")
    cat("      Install: install.packages('rbenchmark')\n\n")
}

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
wkdir <- file.path(home_dir, "Group13")
data_file <- file.path(wkdir, "simple_benchmark.rda")
setwd(wkdir)

cat("\n============================================================\n")
cat("Comprehensive Performance Comparison\n")
cat("============================================================\n")
cat("Working directory: ", getwd(), "\n")

# Load test data
load(data_file)
cat("Loaded: config, dmi, samples\n")
cat("Number of samples: ", samples@nmc, "\n")
cat("Number of chains: ", samples@nchain, "\n\n")

# ============================================================
# Define wrapper functions
# ============================================================

run_original_twice <- function() {
    fit0 <- ggdmc:::run_subject(config, dmi, samples)
    fit1 <- ggdmc:::run_subject(config, dmi, samples = fit0)
    return(fit1)
}

run_fast_with_timing_twice <- function() {
    fit0 <- ggdmc:::run_subject_fast(config, dmi, samples)
    fit1 <- ggdmc:::run_subject_fast(config, dmi, samples = fit0$posterior)
    return(fit1)
}

run_fast_notiming_twice <- function() {
    fit0 <- ggdmc:::run_subject_fast_notiming(config, dmi, samples)
    fit1 <- ggdmc:::run_subject_fast_notiming(config, dmi, samples = fit0)
    return(fit1)
}

# ============================================================
# Run Benchmark
# ============================================================

if (use_rbenchmark) {
    cat("Using rbenchmark for performance comparison...\n")
    cat("Running each version 5 times for stability\n\n")

    results <- rbenchmark::benchmark(
        original = run_original_twice(),
        fast_with_timing = run_fast_with_timing_twice(),
        fast_notiming = run_fast_notiming_twice(),
        replications = 5,
        columns = c("test", "replications", "elapsed", "relative")
    )

    cat("============================================================\n")
    cat("Benchmark Results\n")
    cat("============================================================\n")
    print(results)
    cat("\n")

    # Extract times
    time_original <- results[results$test == "original", "elapsed"]
    time_with_timing <- results[results$test == "fast_with_timing", "elapsed"]
    time_notiming <- results[results$test == "fast_notiming", "elapsed"]

    cat("------------------------------------------------------------\n")
    cat("Summary:\n")
    cat("------------------------------------------------------------\n")
    cat(sprintf("Original:              %.3f seconds  (baseline)\n", time_original))
    cat(sprintf(
        "Fast with timing:      %.3f seconds  (%.2fx, %+.1f%%)\n",
        time_with_timing,
        time_original / time_with_timing,
        100 * (time_with_timing - time_original) / time_original
    ))
    cat(sprintf(
        "Fast without timing:   %.3f seconds  (%.2fx, %+.1f%%)\n",
        time_notiming,
        time_original / time_notiming,
        100 * (time_notiming - time_original) / time_original
    ))

    cat("\n")
    timing_overhead <- time_with_timing - time_notiming
    cat(sprintf(
        "Timing overhead:       %.3f seconds  (%.1f%% of fast version)\n",
        timing_overhead,
        100 * timing_overhead / time_with_timing
    ))
} else {
    # Fallback to system.time()
    cat("Using system.time() for performance comparison...\n")
    cat("Running each version 5 times\n\n")

    run_and_time <- function(func, name, n = 5) {
        cat(sprintf("Running %s (%d times)...\n", name, n))
        times <- numeric(n)
        for (i in 1:n) {
            cat("  Run", i, "... ")
            timing <- system.time({
                result <- func()
            })
            times[i] <- timing["elapsed"]
            cat(sprintf("%.3f seconds\n", timing["elapsed"]))
        }
        return(list(mean = mean(times), sd = sd(times), times = times))
    }

    orig <- run_and_time(run_original_twice, "ORIGINAL", 5)
    fast_timing <- run_and_time(run_fast_with_timing_twice, "FAST with timing", 5)
    fast_notiming <- run_and_time(run_fast_notiming_twice, "FAST no timing", 5)

    cat("\n============================================================\n")
    cat("Benchmark Results\n")
    cat("============================================================\n")
    cat(sprintf("Original:            %.3f ± %.3f seconds\n", orig$mean, orig$sd))
    cat(sprintf("Fast with timing:    %.3f ± %.3f seconds\n", fast_timing$mean, fast_timing$sd))
    cat(sprintf("Fast without timing: %.3f ± %.3f seconds\n", fast_notiming$mean, fast_notiming$sd))

    cat("\n------------------------------------------------------------\n")
    cat("Summary:\n")
    cat("------------------------------------------------------------\n")
    cat(sprintf("Original:              %.3f seconds  (baseline)\n", orig$mean))
    cat(sprintf(
        "Fast with timing:      %.3f seconds  (%.2fx, %+.1f%%)\n",
        fast_timing$mean,
        orig$mean / fast_timing$mean,
        100 * (fast_timing$mean - orig$mean) / orig$mean
    ))
    cat(sprintf(
        "Fast without timing:   %.3f seconds  (%.2fx, %+.1f%%)\n",
        fast_notiming$mean,
        orig$mean / fast_notiming$mean,
        100 * (fast_notiming$mean - orig$mean) / orig$mean
    ))

    timing_overhead <- fast_timing$mean - fast_notiming$mean
    cat(sprintf(
        "\nTiming overhead:       %.3f seconds  (%.1f%% of fast version)\n",
        timing_overhead,
        100 * timing_overhead / fast_timing$mean
    ))

    time_original <- orig$mean
    time_with_timing <- fast_timing$mean
    time_notiming <- fast_notiming$mean
}

# ============================================================
# Detailed Analysis
# ============================================================
cat("\n============================================================\n")
cat("Detailed Analysis\n")
cat("============================================================\n")

# Get C++ timing from one run
fit_timing <- run_fast_with_timing_twice()
timing <- fit_timing$timing

cat("\nC++ Internal Timing (from instrumented version):\n")
cat(sprintf("  Total C++ time:      %.3f seconds\n", timing$total_time))
cat(sprintf(
    "  Likelihood time:     %.3f seconds (%.1f%%)\n",
    timing$likelihood_time,
    100 * timing$likelihood_time / timing$total_time
))
cat(sprintf(
    "  Crossover time:      %.3f seconds (%.1f%%)\n",
    timing$crossover_time,
    100 * timing$crossover_time / timing$total_time
))
cat(sprintf(
    "  Migration time:      %.3f seconds (%.1f%%)\n",
    timing$migration_time,
    100 * timing$migration_time / timing$total_time
))

cat("\nWhere the time goes:\n")
cat(sprintf(
    "  C++ computation:     %.3f seconds (%.1f%% of total)\n",
    timing$total_time,
    100 * timing$total_time / time_with_timing
))
cat(sprintf(
    "  R overhead:          %.3f seconds (%.1f%% of total)\n",
    time_with_timing - timing$total_time,
    100 * (time_with_timing - timing$total_time) / time_with_timing
))

cat("\n------------------------------------------------------------\n")
cat("Interpretation:\n")
cat("------------------------------------------------------------\n")

if (time_notiming < time_original) {
    speedup <- 100 * (time_original - time_notiming) / time_original
    cat(sprintf("✓ Optimization SUCCESSFUL: %.1f%% faster\n", speedup))
    if (speedup > 10) {
        cat("  → Significant improvement! Pre-allocated buffers help.\n")
    } else if (speedup > 5) {
        cat("  → Moderate improvement. Pre-allocated buffers have measurable impact.\n")
    } else {
        cat("  → Small improvement. Other bottlenecks may dominate.\n")
    }
} else {
    slowdown <- 100 * (time_notiming - time_original) / time_original
    cat(sprintf("✗ Optimization SLOWER: %.1f%% slower\n", slowdown))
    cat("  → Pre-allocated buffers don't help for this problem size.\n")
}

if (timing$likelihood_time / timing$total_time > 0.7) {
    cat("\n• Likelihood dominates (>70% of C++ time)\n")
    cat("  → Next optimization should focus on likelihood functions\n")
    cat("  → Chain shuffling is NOT the bottleneck\n")
} else {
    cat("\n• DE algorithm overhead is significant\n")
    cat("  → Chain shuffling optimization may help more with larger problems\n")
}

if (time_with_timing > time_notiming * 1.02) {
    cat("\n• Timing instrumentation adds measurable overhead (>2%%)\n")
    cat("  → Use _notiming versions for production\n")
    cat("  → Use _fast versions only for profiling\n")
}

cat("\n============================================================\n")
cat("Recommendation\n")
cat("============================================================\n")

if (time_notiming < time_original) {
    cat("Use run_subject_fast_notiming() for production runs.\n")
    if (time_notiming < time_original * 0.95) {
        cat("Consider proceeding with Step 2 optimizations.\n")
    } else {
        cat("Focus on likelihood function optimization instead.\n")
    }
} else {
    cat("Stick with the original run_subject() function.\n")
    cat("Focus optimization efforts on likelihood functions.\n")
}

cat("\n============================================================\n")
cat("Benchmark Complete!\n")
cat("============================================================\n")
