#!/usr/bin/env Rscript
# Quick test to verify the fast functions are accessible and work

cat("Testing fast DE optimizer functions...\n\n")

library(ggdmc)

# Check if functions are exported
cat("1. Checking if functions are available...\n")
has_subject_fast <- exists("run_subject_fast", where = "package:ggdmc", mode = "function")
has_hyper_fast <- exists("run_hyper_fast", where = "package:ggdmc", mode = "function")

if (has_subject_fast) {
    cat("   ✓ run_subject_fast is available\n")
} else {
    cat("   ✗ run_subject_fast is NOT available\n")
}

if (has_hyper_fast) {
    cat("   ✓ run_hyper_fast is available\n")
} else {
    cat("   ✗ run_hyper_fast is NOT available\n")
}

if (!has_subject_fast || !has_hyper_fast) {
    cat("\nError: Fast functions not found. Did the package install correctly?\n")
    cat("Try: R CMD INSTALL --preclean --no-multiarch --with-keep.source .\n")
    quit(status = 1)
}

cat("\n2. Checking function signatures...\n")
cat("   run_subject_fast signature:\n")
cat("     ", paste(names(formals(ggdmc:::run_subject_fast)), collapse = ", "), "\n")
cat("   run_hyper_fast signature:\n")
cat("     ", paste(names(formals(ggdmc:::run_hyper_fast)), collapse = ", "), "\n")

cat("\n3. Package info...\n")
cat("   Package version:", as.character(packageVersion("ggdmc")), "\n")
cat("   Package path:", system.file(package = "ggdmc"), "\n")

cat("\n✓ All basic checks passed!\n")
cat("\nTo run a full test with actual data:\n")
cat("  1. cd tests/testthat/Group1\n")
cat("  2. R CMD BATCH 0_5param_hyper.r\n")
cat("  3. Or modify the script to use run_subject_fast()\n")

cat("\nTo benchmark performance:\n")
cat("  source('benchmark_fast_version.R')\n")
cat("  (requires test data to be loaded first)\n")
