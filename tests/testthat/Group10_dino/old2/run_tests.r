#!/usr/bin/env Rscript
# Interactive launcher for systematic DINA testing

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║        DINA Model Asymptotic Behavior Testing Suite           ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Available tests:\n")
cat("  1. Quick Test      - Test small datasets only (5-10 min)\n")
cat("  2. Systematic Test - Test all datasets (several hours)\n")
cat("  3. Analyze Results - Analyze and visualize results\n")
cat("  4. Full Workflow   - Run systematic test + analysis\n")
cat("  5. Exit\n\n")

# Get user choice
if (interactive()) {
  choice <- readline(prompt = "Enter your choice (1-5): ")
  choice <- as.integer(choice)
} else {
  # If running non-interactively, check command line args
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    choice <- as.integer(args[1])
  } else {
    cat("Usage: Rscript run_tests.r <choice>\n")
    cat("  where <choice> is 1-5\n\n")
    cat("Running in non-interactive mode without args - defaulting to quick test\n\n")
    choice <- 1
  }
}

# Validate choice
if (is.na(choice) || choice < 1 || choice > 5) {
  cat("Invalid choice. Exiting.\n\n")
  quit(save = "no", status = 1)
}

# Execute based on choice
switch(choice,
  # 1. Quick Test
  {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("  Running Quick Test\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    source("01_fit_subject_dina0_quick.r")
  },

  # 2. Systematic Test
  {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("  Running Systematic Test (this will take time!)\n")
    cat("═══════════════════════════════════════════════════════════════\n")

    if (interactive()) {
      confirm <- readline(prompt = "This may take several hours. Continue? (yes/no): ")
      if (tolower(confirm) != "yes") {
        cat("Cancelled.\n\n")
        quit(save = "no")
      }
    }

    source("01_fit_subject_dina0_systematic.r")
  },

  # 3. Analyze Results
  {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("  Analyzing Results\n")
    cat("═══════════════════════════════════════════════════════════════\n")

    # Check if results exist
    home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
    results_dir <- file.path(home_dir, "Group10_cdm_subjects/results")
    summary_path <- file.path(results_dir, "systematic_test_summary.rda")

    if (!file.exists(summary_path)) {
      cat("\n⚠ WARNING: No results found!\n")
      cat("Please run systematic test first (option 2)\n\n")
      quit(save = "no", status = 1)
    }

    source("02_analyze_asymptotic_results.r")
  },

  # 4. Full Workflow
  {
    cat("\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("  Running Full Workflow\n")
    cat("═══════════════════════════════════════════════════════════════\n")

    if (interactive()) {
      confirm <- readline(prompt = "This may take several hours. Continue? (yes/no): ")
      if (tolower(confirm) != "yes") {
        cat("Cancelled.\n\n")
        quit(save = "no")
      }
    }

    cat("\nStep 1/2: Systematic Testing\n")
    source("01_fit_subject_dina0_systematic.r")

    cat("\nStep 2/2: Analysis\n")
    source("02_analyze_asymptotic_results.r")

    cat("\n")
    cat("╔════════════════════════════════════════════════════════════════╗\n")
    cat("║              FULL WORKFLOW COMPLETED                           ║\n")
    cat("╚════════════════════════════════════════════════════════════════╝\n\n")
  },

  # 5. Exit
  {
    cat("\nExiting.\n\n")
    quit(save = "no")
  }
)

cat("\n✓ Done\n\n")
