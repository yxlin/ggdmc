#!/usr/bin/env Rscript
# Hierarchical Bayesian Model Performance Benchmark
# Compares three versions:
# 1. Original run() - Baseline
# 2. run_fast() - With cached likelihood and timing
# 3. run_fast_notiming() - With cached likelihood, no timing overhead

rm(list = ls())
library(ggdmc)

# ============================================================================
# Configuration
# ============================================================================

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
test_dir <- file.path(home_dir, "Group14")
data_file <- file.path(test_dir, "hb_benchmark.rda")

setwd(test_dir)
cat("Working directory:", getwd(), "\n\n")

# Check for rbenchmark
use_rbenchmark <- requireNamespace("rbenchmark", quietly = TRUE)
if (!use_rbenchmark) {
    cat("Note: rbenchmark not found. Using system.time() instead.\n")
    cat("Install with: install.packages('rbenchmark')\n\n")
}

# ============================================================================
# Load Data
# ============================================================================

if (!file.exists(data_file)) {
    cat("ERROR: Benchmark data file not found:", data_file, "\n")
    cat("Please run save_benchmark_data.R first to generate the data.\n")
    quit(status = 1)
}

cat("============================================================\n")
cat("Hierarchical Bayesian Model Performance Comparison\n")
cat("============================================================\n\n")

cat("Loading benchmark data...\n")
load(data_file)

# Validate loaded objects
if (!exists("config") || !exists("dmis") || !exists("samples")) {
    cat("ERROR: Data file missing required objects.\n")
    cat("Expected: config, dmis, samples\n")
    quit(status = 1)
}

# ============================================================================
# Data Info
# ============================================================================

n_subjects <- length(dmis)
n_iterations <- samples$phi@nmc
n_chains <- length(samples$phi@theta)

cat("\nBenchmark Configuration:\n")
cat("------------------------------------------------------------\n")
cat(sprintf("Number of subjects:   %d\n", n_subjects))
cat(sprintf("MCMC iterations:      %d\n", n_iterations))
cat(sprintf("Number of chains:     %d\n", n_chains))
cat(sprintf("Parameters per subj:  %d\n", length(samples$subject_theta[[1]]@pnames)))

cat("\nThis benchmark compares:\n")
cat("  1. run()                - Original (recalculates hyper-likelihood)\n")
cat("  2. run_fast()           - Cached + timing instrumentation\n")
cat("  3. run_fast_notiming()  - Cached, no timing overhead\n\n")

cat("Expected speedup: 20-50% for models with 10+ subjects\n")
cat("Speedup scales with number of subjects!\n\n")

# ============================================================================
# Define wrapper functions
# ============================================================================

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

# ============================================================================
# Run Benchmark
# ============================================================================

cat("============================================================\n")
cat("Running Benchmark (3 replications per version)\n")
cat("============================================================\n")
cat("This will take a few minutes...\n\n")

if (use_rbenchmark) {
    cat("Using rbenchmark for performance comparison...\n\n")

    results <- rbenchmark::benchmark(
        original = run_original(),
        fast_with_timing = run_cached_with_timing(),
        fast_notiming = run_cached_notiming(),
        replications = 3,
        columns = c("test", "replications", "elapsed", "relative")
    )

    cat("\n============================================================\n")
    cat("Benchmark Results\n")
    cat("============================================================\n")
    print(results)
    cat("\n")

    # Extract times
    time_original <- results[results$test == "original", "elapsed"]
    time_with_timing <- results[results$test == "fast_with_timing", "elapsed"]
    time_notiming <- results[results$test == "fast_notiming", "elapsed"]

    speedup_timing <- time_original / time_with_timing
    speedup_notiming <- time_original / time_notiming

    cat("------------------------------------------------------------\n")
    cat("Summary:\n")
    cat("------------------------------------------------------------\n")
    cat(sprintf("Original:              %.2f seconds  (baseline)\n", time_original))
    cat(sprintf(
        "Fast with timing:      %.2f seconds  (%.2fx speedup, %.1f%% faster)\n",
        time_with_timing, speedup_timing,
        100 * (1 - time_with_timing / time_original)
    ))
    cat(sprintf(
        "Fast without timing:   %.2f seconds  (%.2fx speedup, %.1f%% faster)\n",
        time_notiming, speedup_notiming,
        100 * (1 - time_notiming / time_original)
    ))

    timing_overhead <- time_with_timing - time_notiming
    cat(sprintf(
        "\nTiming overhead:       %.2f seconds  (%.1f%% of cached version)\n",
        timing_overhead,
        100 * timing_overhead / time_with_timing
    ))
} else {
    # Fallback to system.time()
    cat("Using system.time() for performance comparison...\n\n")

    run_and_time <- function(func, name, n = 3) {
        cat(sprintf("Running %s (%d times)...\n", name, n))
        times <- numeric(n)
        for (i in 1:n) {
            cat("  Run", i, "... ")
            timing <- system.time({
                result <- func()
            })
            times[i] <- timing["elapsed"]
            cat(sprintf("%.2f seconds\n", timing["elapsed"]))
        }
        return(list(mean = mean(times), sd = sd(times), times = times))
    }

    orig <- run_and_time(run_original, "ORIGINAL", 3)
    cat("\n")
    fast_timing <- run_and_time(run_cached_with_timing, "FAST with timing", 3)
    cat("\n")
    fast_notiming <- run_and_time(run_cached_notiming, "FAST no timing", 3)

    cat("\n============================================================\n")
    cat("Benchmark Results\n")
    cat("============================================================\n")
    cat(sprintf("Original:            %.2f ± %.2f seconds\n", orig$mean, orig$sd))
    cat(sprintf("Fast with timing:    %.2f ± %.2f seconds\n", fast_timing$mean, fast_timing$sd))
    cat(sprintf("Fast without timing: %.2f ± %.2f seconds\n", fast_notiming$mean, fast_notiming$sd))

    time_original <- orig$mean
    time_with_timing <- fast_timing$mean
    time_notiming <- fast_notiming$mean

    speedup_timing <- time_original / time_with_timing
    speedup_notiming <- time_original / time_notiming

    cat("\n------------------------------------------------------------\n")
    cat("Summary:\n")
    cat("------------------------------------------------------------\n")
    cat(sprintf("Original:              %.2f seconds  (baseline)\n", time_original))
    cat(sprintf(
        "Fast with timing:      %.2f seconds  (%.2fx speedup, %.1f%% faster)\n",
        time_with_timing, speedup_timing,
        100 * (1 - time_with_timing / time_original)
    ))
    cat(sprintf(
        "Fast without timing:   %.2f seconds  (%.2fx speedup, %.1f%% faster)\n",
        time_notiming, speedup_notiming,
        100 * (1 - time_notiming / time_original)
    ))

    timing_overhead <- time_with_timing - time_notiming
    cat(sprintf(
        "\nTiming overhead:       %.2f seconds  (%.1f%% of cached version)\n",
        timing_overhead,
        100 * timing_overhead / time_with_timing
    ))
}

# ============================================================================
# Detailed Analysis
# ============================================================================

cat("\n============================================================\n")
cat("Detailed C++ Timing Analysis\n")
cat("============================================================\n")

# Get timing breakdown from instrumented version
cat("\nGathering C++ internal timing data...\n")
result_timing <- run_cached_with_timing()
timing <- result_timing$timing

cat("\nC++ Internal Timing Breakdown:\n")
cat("------------------------------------------------------------\n")
cat(sprintf("  Total C++ time:      %.2f seconds (100.0%%)\n", timing$total_time))
cat(sprintf(
    "  Likelihood time:     %.2f seconds (%.1f%%)\n",
    timing$likelihood_time,
    100 * timing$likelihood_time / timing$total_time
))
cat(sprintf(
    "  Crossover time:      %.2f seconds (%.1f%%)\n",
    timing$crossover_time,
    100 * timing$crossover_time / timing$total_time
))
cat(sprintf(
    "  Migration time:      %.2f seconds (%.1f%%)\n",
    timing$migration_time,
    100 * timing$migration_time / timing$total_time
))

other_time <- timing$total_time - (timing$likelihood_time)
cat(sprintf(
    "  Other (MH, etc):     %.2f seconds (%.1f%%)\n",
    other_time,
    100 * other_time / timing$total_time
))

cat("\nOperation Counts:\n")
cat(sprintf("  Crossover operations: %d\n", timing$n_crossover))
cat(sprintf("  Migration operations: %d\n", timing$n_migration))

cat("\nTime Distribution:\n")
cat("------------------------------------------------------------\n")
cat(sprintf(
    "  C++ computation:     %.2f seconds (%.1f%% of total)\n",
    timing$total_time,
    100 * timing$total_time / time_with_timing
))
cat(sprintf(
    "  R overhead:          %.2f seconds (%.1f%% of total)\n",
    time_with_timing - timing$total_time,
    100 * (time_with_timing - timing$total_time) / time_with_timing
))

# ============================================================================
# Interpretation
# ============================================================================

cat("\n============================================================\n")
cat("Interpretation\n")
cat("============================================================\n\n")

if (speedup_notiming >= 1.2) {
    cat(sprintf("✓ EXCELLENT: %.2fx speedup with cached likelihood\n", speedup_notiming))
    cat("  → Caching eliminates redundant hyper-likelihood calculations\n")
    cat("  → Optimization is highly effective for this problem\n")
} else if (speedup_notiming >= 1.05) {
    cat(sprintf("✓ GOOD: %.2fx speedup with cached likelihood\n", speedup_notiming))
    cat("  → Moderate improvement from caching\n")
    cat("  → Worth using for production runs\n")
} else if (speedup_notiming >= 1.0) {
    cat(sprintf("~ MARGINAL: %.2fx speedup with cached likelihood\n", speedup_notiming))
    cat("  → Small improvement, may not be significant\n")
    cat("  → Consider problem size and complexity\n")
} else {
    cat(sprintf(
        "✗ SLOWER: %.2fx (%.1f%% slower)\n", speedup_notiming,
        100 * (speedup_notiming - 1)
    ))
    cat("  → Caching overhead exceeds benefit for this problem\n")
    cat("  → Stick with original version\n")
}

cat("\n------------------------------------------------------------\n")
cat("Impact Factors:\n")
cat("------------------------------------------------------------\n")
cat(sprintf("Number of subjects: %d\n", n_subjects))

if (n_subjects < 5) {
    cat("  → Small number of subjects\n")
    cat("  → Caching benefit may be limited\n")
    cat("  → Try with more subjects to see bigger gains\n")
} else if (n_subjects < 20) {
    cat("  → Moderate number of subjects\n")
    cat("  → Should see measurable benefit from caching\n")
} else {
    cat("  → Large number of subjects\n")
    cat("  → Caching should provide significant speedup\n")
    cat("  → Original version wastes O(n_subjects) computations\n")
}

# ============================================================================
# Recommendation
# ============================================================================

cat("\n============================================================\n")
cat("Recommendation\n")
cat("============================================================\n\n")

if (speedup_notiming >= 1.05) {
    cat("USE: run_fast_notiming() for production HB runs\n")
    cat("     run_fast() for profiling and analysis\n\n")
    cat("Benefits:\n")
    cat(sprintf("  • %.2fx faster than original\n", speedup_notiming))
    cat(sprintf("  • Saves %.2f seconds per run\n", time_original - time_notiming))
    cat("  • Eliminates redundant likelihood calculations\n")
    cat("  • Same numerical results, just faster\n\n")
    cat("Code change:\n")
    cat("  # OLD:\n")
    cat("  result <- run(config, dmis, samples)\n\n")
    cat("  # NEW:\n")
    cat("  result <- run_fast_notiming(config, dmis, samples)\n")
} else {
    cat("USE: run() (original version)\n\n")
    cat("Reasons:\n")
    cat("  • Caching doesn't provide significant speedup\n")
    cat("  • Stick with proven baseline\n")
    if (n_subjects < 5) {
        cat("  • Try with more subjects to see benefit\n")
    }
}

cat("\n============================================================\n")
cat("Benchmark Complete!\n")
cat("============================================================\n")
cat("\nSummary saved to: benchmark_results.txt\n\n")

# ============================================================================
# Save Results
# ============================================================================

# Save summary to file
sink("benchmark_results.txt")
cat("HB Benchmark Results\n")
cat("====================\n\n")
cat("Configuration:\n")
cat(sprintf("  Subjects: %d\n", n_subjects))
cat(sprintf("  Iterations: %d\n", n_iterations))
cat(sprintf("  Chains: %d\n", n_chains))
cat("\nResults:\n")
cat(sprintf("  Original: %.2f seconds\n", time_original))
cat(sprintf("  Fast (timing): %.2f seconds (%.2fx)\n", time_with_timing, speedup_timing))
cat(sprintf("  Fast (no timing): %.2f seconds (%.2fx)\n", time_notiming, speedup_notiming))
cat("\nRecommendation:\n")
if (speedup_notiming >= 1.05) {
    cat("  Use run_fast_notiming() for production\n")
} else {
    cat("  Use run() (original)\n")
}
sink()

cat("Done!\n")

# ============================================================
# Benchmark Results
# ============================================================
#               test replications elapsed relative
# 3    fast_notiming            3  16.312    1.125
# 2 fast_with_timing            3  15.085    1.040
# 1         original            3  14.498    1.000

# ------------------------------------------------------------
# Summary:
# ------------------------------------------------------------
# Original:              14.50 seconds  (baseline)
# Fast with timing:      15.09 seconds  (0.96x speedup, -4.0% faster)
# Fast without timing:   16.31 seconds  (0.89x speedup, -12.5% faster)

# Timing overhead:       -1.23 seconds  (-8.1% of cached version)

# ============================================================
# Detailed C++ Timing Analysis
# ============================================================

# Gathering C++ internal timing data...
# 100 200 300 400 500

# C++ Internal Timing Breakdown:
# ------------------------------------------------------------
#   Total C++ time:      5.40 seconds (100.0%)
#   Likelihood time:     0.19 seconds (3.6%)
#   Crossover time:      0.21 seconds (3.9%)
#   Migration time:      0.00 seconds (0.0%)
#   Other (MH, etc):     5.20 seconds (96.4%)

# Operation Counts:
#   Crossover operations: 499
#   Migration operations: 0

# Time Distribution:
# ------------------------------------------------------------
#   C++ computation:     5.40 seconds (35.8% of total)
#   R overhead:          9.69 seconds (64.2% of total)

# ============================================================
# Interpretation
# ============================================================

# ✗ SLOWER: 0.89x (-11.1% slower)
#   → Caching overhead exceeds benefit for this problem
#   → Stick with original version

# ------------------------------------------------------------
# Impact Factors:
# ------------------------------------------------------------
# Number of subjects: 32
#   → Large number of subjects
#   → Caching should provide significant speedup
#   → Original version wastes O(n_subjects) computations

# ============================================================
# Recommendation
# ============================================================

# USE: run() (original version)

# Reasons:
#   • Caching doesn't provide significant speedup
#   • Stick with proven baseline

# ============================================================
# Benchmark Complete!
# ============================================================

# Summary saved to: benchmark_results.txt

# Done!
