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
save_file <- file.path(home_dir, "Group14/hb_benchmark.rda")
if (!file.exists(data_file)) {
    cat("ERROR: Data file not found:", data_file, "\n")
    cat("Please modify the paths in this script to point to your HB test data.\n")
    quit(status = 1)
}

setwd(test_dir)
load(data_file)

cat("Please set up your hierarchical model before running this benchmark.\n")
cat("Required objects: config, dmis (list), samples (list with phi and subject_theta)\n\n")

cat("Example minimal setup:\n")
cat("------------------------------------------------------------\n")
cat("# Assume you have your HB model ready\n")


theta_input <- pop_theta_input
pop_migration_prob = 0.00
sub_migration_prob = 0.00

gamma_precursor = 2.38
rp = 0.001
is_hblocked = FALSE
is_pblocked = FALSE

nchain = NULL
thin = 1L
nparameter <- pop_priors@nparameter
pnames <- pop_priors@pnames
nchain <- ggdmc:::.get_nchain(pop_priors, nchain)
pop_debug = FALSE
sub_debug = FALSE
de_input <- setDEInput(
        pop_migration_prob = pop_migration_prob, 
        sub_migration_prob = sub_migration_prob, gamma_precursor = gamma_precursor,
        rp = rp, is_hblocked = is_hblocked, is_pblocked = is_pblocked,
        nparameter = as.integer(nparameter), nchain = as.integer(nchain),
        pop_debug = pop_debug, sub_debug = sub_debug
    )

ncore = 1L
seed = 123


config_list <- set_configs(
        prior = pop_priors, theta_input = theta_input, de_input = de_input,
        ncore = ncore, seed = seed
    )
dmis <- pop_dmis # N subjects
samples <- pop_samples
config <- config_list[[1]]

save(config, dmis, samples, file = save_file)

# res <- run(config_list[[1]], pop_dmis, pop_samples)
# phi <- res$phi
# subject_theta <- res$subject_theta
# options(digits = 2)
# est_phi <- compare(phi, ps = true_vector, start = phi@nmc * 0.5)

# DT <- prepare_thetas_data(subject_theta, start = phi@nmc * 0.5)
# p1 <- plot_thetas(DT)

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

# result <- ggdmc::run(config, dmis, samples)

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
    

    cat(\"\\n============================================================\\n\")
    cat(\"Benchmark Results\\n\")
    cat(\"============================================================\\n\")
    print(results)
    cat(\"\\n\")

    # Extract times
    time_original <- results[results$test == \"original\", \"elapsed\"]
    time_with_timing <- results[results$test == \"fast_with_timing\", \"elapsed\"]
    time_notiming <- results[results$test == \"fast_notiming\", \"elapsed\"]

    speedup_timing <- time_original / time_with_timing
    speedup_notiming <- time_original / time_notiming

    cat(\"------------------------------------------------------------\\n\")
    cat(\"Summary:\\n\")
    cat(\"------------------------------------------------------------\\n\")
    cat(sprintf(\"Original:              %.2f seconds  (baseline)\\n\", time_original))
    cat(sprintf(\"Fast with timing:      %.2f seconds  (%.2fx speedup, %.1f%% faster)\\n\",
                time_with_timing, speedup_timing,
                100 * (1 - time_with_timing/time_original)))
    cat(sprintf(\"Fast without timing:   %.2f seconds  (%.2fx speedup, %.1f%% faster)\\n\",
                time_notiming, speedup_notiming,
                100 * (1 - time_notiming/time_original)))

    timing_overhead <- time_with_timing - time_notiming
    cat(sprintf(\"\\nTiming overhead:       %.2f seconds  (%.1f%% of cached version)\\n\",
                timing_overhead,
                100 * timing_overhead / time_with_timing))

} else {
    # Fallback to system.time()
    cat(\"Using system.time() for performance comparison...\\n\\n\")

    run_and_time <- function(func, name, n = 3) {
        cat(sprintf(\"Running %s (%d times)...\\n\", name, n))
        times <- numeric(n)
        for (i in 1:n) {
            cat(\"  Run\", i, \"... \")
            timing <- system.time({ result <- func() })
            times[i] <- timing[\"elapsed\"]
            cat(sprintf(\"%.2f seconds\\n\", timing[\"elapsed\"]))
        }
        return(list(mean = mean(times), sd = sd(times), times = times))
    }

    orig <- run_and_time(run_original, \"ORIGINAL\", 3)
    cat(\"\\n\")
    fast_timing <- run_and_time(run_cached_with_timing, \"FAST with timing\", 3)
    cat(\"\\n\")
    fast_notiming <- run_and_time(run_cached_notiming, \"FAST no timing\", 3)

    cat(\"\\n============================================================\\n\")
    cat(\"Benchmark Results\\n\")
    cat(\"============================================================\\n\")
    cat(sprintf(\"Original:            %.2f ± %.2f seconds\\n\", orig$mean, orig$sd))
    cat(sprintf(\"Fast with timing:    %.2f ± %.2f seconds\\n\", fast_timing$mean, fast_timing$sd))
    cat(sprintf(\"Fast without timing: %.2f ± %.2f seconds\\n\", fast_notiming$mean, fast_notiming$sd))

    time_original <- orig$mean
    time_with_timing <- fast_timing$mean
    time_notiming <- fast_notiming$mean

    speedup_timing <- time_original / time_with_timing
    speedup_notiming <- time_original / time_notiming

    cat(\"\\n------------------------------------------------------------\\n\")
    cat(\"Summary:\\n\")
    cat(\"------------------------------------------------------------\\n\")
    cat(sprintf(\"Original:              %.2f seconds  (baseline)\\n\", time_original))
    cat(sprintf(\"Fast with timing:      %.2f seconds  (%.2fx speedup, %.1f%% faster)\\n\",
                time_with_timing, speedup_timing,
                100 * (1 - time_with_timing/time_original)))
    cat(sprintf(\"Fast without timing:   %.2f seconds  (%.2fx speedup, %.1f%% faster)\\n\",
                time_notiming, speedup_notiming,
                100 * (1 - time_notiming/time_original)))

    timing_overhead <- time_with_timing - time_notiming
    cat(sprintf(\"\\nTiming overhead:       %.2f seconds  (%.1f%% of cached version)\\n\",
                timing_overhead,
                100 * timing_overhead / time_with_timing))
}

# ===========================================================================
# Detailed Analysis
# ===========================================================================

cat(\"\\n============================================================\\n\")
cat(\"Detailed Analysis\\n\")
cat(\"============================================================\\n\")

# Get timing breakdown from instrumented version
cat(\"\\nRunning one more iteration with timing instrumentation...\\n\")
result_timing <- run_cached_with_timing()
timing <- result_timing$timing

cat(\"\\nC++ Internal Timing (from instrumented version):\\n\")
cat(sprintf(\"  Total C++ time:      %.2f seconds (100.0%%)\\n\", timing$total_time))
cat(sprintf(\"  Likelihood time:     %.2f seconds (%.1f%%)\\n\",
            timing$likelihood_time,
            100 * timing$likelihood_time / timing$total_time))
cat(sprintf(\"  Crossover time:      %.2f seconds (%.1f%%)\\n\",
            timing$crossover_time,
            100 * timing$crossover_time / timing$total_time))
cat(sprintf(\"  Migration time:      %.2f seconds (%.1f%%)\\n\",
            timing$migration_time,
            100 * timing$migration_time / timing$total_time))
cat(sprintf(\"  Operations: %d crossovers, %d migrations\\n\",
            timing$n_crossover, timing$n_migration))

cat(\"\\nWhere the time goes:\\n\")
cat(sprintf(\"  C++ computation:     %.2f seconds (%.1f%% of total)\\n\",
            timing$total_time,
            100 * timing$total_time / time_with_timing))
cat(sprintf(\"  R overhead:          %.2f seconds (%.1f%% of total)\\n\",
            time_with_timing - timing$total_time,
            100 * (time_with_timing - timing$total_time) / time_with_timing))

# ===========================================================================
# Interpretation and Recommendations
# ===========================================================================

cat(\"\\n============================================================\\n\")
cat(\"Interpretation\\n\")
cat(\"============================================================\\n\")

if (speedup_notiming >= 1.2) {
    cat(sprintf(\"✓ EXCELLENT: %.1fx speedup with cached likelihood\\n\", speedup_notiming))
    cat(\"  → Caching eliminates redundant hyper-likelihood calculations\\n\")
    cat(\"  → Optimization is highly effective for this problem\\n\")
} else if (speedup_notiming >= 1.05) {
    cat(sprintf(\"✓ GOOD: %.1fx speedup with cached likelihood\\n\", speedup_notiming))
    cat(\"  → Moderate improvement from caching\\n\")
    cat(\"  → Worth using for production runs\\n\")
} else if (speedup_notiming >= 1.0) {
    cat(sprintf(\"~ MARGINAL: %.1fx speedup with cached likelihood\\n\", speedup_notiming))
    cat(\"  → Small improvement, may not be significant\\n\")
    cat(\"  → Consider problem size and complexity\\n\")
} else {
    cat(sprintf(\"✗ SLOWER: %.1fx (%.1f%% slower)\\n\", speedup_notiming,
                100 * (1 - speedup_notiming)))
    cat(\"  → Caching overhead exceeds benefit for this problem\\n\")
    cat(\"  → Stick with original version\\n\")
}

cat(\"\\n------------------------------------------------------------\\n\")
cat(\"Impact Factors:\\n\")
cat(\"------------------------------------------------------------\\n\")
cat(sprintf(\"Number of subjects: %d\\n\", n_subjects))

if (n_subjects < 5) {
    cat(\"  → Small number of subjects\\n\")
    cat(\"  → Caching benefit may be limited\\n\")
    cat(\"  → Try with more subjects to see bigger gains\\n\")
} else if (n_subjects < 20) {
    cat(\"  → Moderate number of subjects\\n\")
    cat(\"  → Should see measurable benefit from caching\\n\")
} else {
    cat(\"  → Large number of subjects\\n\")
    cat(\"  → Caching should provide significant speedup\\n\")
    cat(\"  → Original version wastes O(n_subjects) computations\\n\")
}

cat(\"\\n============================================================\\n\")
cat(\"Recommendation\\n\")
cat(\"============================================================\\n\")

if (speedup_notiming >= 1.05) {
    cat(\"\\nUSE: run_fast_notiming() for production HB runs\\n\")
    cat(\"     run_fast() for profiling and analysis\\n\\n\")
    cat(\"Benefits:\\n\")
    cat(sprintf(\"  • %.1fx faster than original\\n\", speedup_notiming))
    cat(\"  • Eliminates redundant likelihood calculations\\n\")
    cat(\"  • Same numerical results, just faster\\n\")
} else {
    cat(\"\\nUSE: run() (original version)\\n\\n\")
    cat(\"Reasons:\\n\")
    cat(\"  • Caching doesn't provide significant speedup\\n\")
    cat(\"  • Stick with proven baseline\\n\")
    cat(\"  • May be more effective with larger models\\n\")
}

cat(\"\\n============================================================\\n\")
cat(\"Benchmark Complete!\\n\")
cat(\"============================================================\\n\")
