#!/usr/bin/env Rscript
# Benchmark script to compare original vs fast DE optimizer versions
# This script compares the performance of the pre-allocated chain indices optimization

library(ggdmc)

cat("============================================================\n")
cat("Benchmark: Original vs Fast DE Optimizer (Chain Shuffling)\n")
cat("============================================================\n\n")

# Load or create a test dataset
# This assumes you have a working test case in your test suite
# Adjust the path as needed

cat("Setting up test data...\n")

# Example: Load a simple LBA model test case
# You'll need to adjust this based on your actual test data structure
tryCatch({
    # Try to load from test data if available
    source("tests/testthat/Group1/data/0_lba_model0.r")
    cat("Loaded test data successfully\n")
}, error = function(e) {
    cat("Could not load test data. Please run this script from your test directory.\n")
    cat("Error:", conditionMessage(e), "\n")
    quit(status = 1)
})

cat("\n============================================================\n")
cat("Running ORIGINAL version...\n")
cat("============================================================\n")

# Time the original version
system.time({
    result_original <- ggdmc:::run_subject(config, dmi, samples)
}) -> time_original

cat("\nOriginal version completed\n")
cat("User time:  ", time_original["user.self"], "seconds\n")
cat("System time:", time_original["sys.self"], "seconds\n")
cat("Elapsed:    ", time_original["elapsed"], "seconds\n")

cat("\n============================================================\n")
cat("Running FAST version (with pre-allocated buffers)...\n")
cat("============================================================\n")

# Time the fast version
system.time({
    result_fast <- ggdmc:::run_subject_fast(config, dmi, samples)
}) -> time_fast

cat("\nFast version completed\n")
cat("User time:  ", time_fast["user.self"], "seconds\n")
cat("System time:", time_fast["sys.self"], "seconds\n")
cat("Elapsed:    ", time_fast["elapsed"], "seconds\n")

cat("\n============================================================\n")
cat("C++ Internal Timing Statistics (Fast Version)\n")
cat("============================================================\n")
cat("Total time:      ", result_fast$timing$total_time, "seconds\n")
cat("Crossover time:  ", result_fast$timing$crossover_time, "seconds\n")
cat("Migration time:  ", result_fast$timing$migration_time, "seconds\n")
cat("Likelihood time: ", result_fast$timing$likelihood_time, "seconds\n")
cat("N crossovers:    ", result_fast$timing$n_crossover, "\n")
cat("N migrations:    ", result_fast$timing$n_migration, "\n")

cat("\n============================================================\n")
cat("Performance Comparison\n")
cat("============================================================\n")

speedup <- time_original["elapsed"] / time_fast["elapsed"]
cat("Speedup factor:    ", sprintf("%.3fx", speedup), "\n")
cat("Time saved:        ", sprintf("%.3f seconds", time_original["elapsed"] - time_fast["elapsed"]), "\n")
cat("Percentage faster: ", sprintf("%.2f%%", (1 - time_fast["elapsed"]/time_original["elapsed"]) * 100), "\n")

cat("\n============================================================\n")
cat("Time Breakdown (Fast Version)\n")
cat("============================================================\n")

total_cpp <- result_fast$timing$total_time
total_elapsed <- time_fast["elapsed"]
overhead <- total_elapsed - total_cpp

cat("Total R elapsed:     ", sprintf("%.3f seconds", total_elapsed), "\n")
cat("Total C++ work:      ", sprintf("%.3f seconds (%.1f%%)", total_cpp, 100*total_cpp/total_elapsed), "\n")
cat("R/C++ overhead:      ", sprintf("%.3f seconds (%.1f%%)", overhead, 100*overhead/total_elapsed), "\n")
cat("\nC++ Time Breakdown:\n")
cat("  Crossover:         ", sprintf("%.3f seconds (%.1f%%)",
                                     result_fast$timing$crossover_time,
                                     100*result_fast$timing$crossover_time/total_cpp), "\n")
cat("  Migration:         ", sprintf("%.3f seconds (%.1f%%)",
                                     result_fast$timing$migration_time,
                                     100*result_fast$timing$migration_time/total_cpp), "\n")
cat("  Likelihood:        ", sprintf("%.3f seconds (%.1f%%)",
                                     result_fast$timing$likelihood_time,
                                     100*result_fast$timing$likelihood_time/total_cpp), "\n")

other_cpp <- total_cpp - (result_fast$timing$crossover_time +
                          result_fast$timing$migration_time)
cat("  Other (MH, etc):   ", sprintf("%.3f seconds (%.1f%%)", other_cpp, 100*other_cpp/total_cpp), "\n")

cat("\n============================================================\n")
cat("Verification: Results should be statistically similar\n")
cat("============================================================\n")
cat("(Check that both versions produce valid MCMC chains)\n")
cat("Original posterior class:", class(result_original), "\n")
cat("Fast posterior class:    ", class(result_fast$posterior), "\n")

cat("\nBenchmark complete!\n")
