#!/usr/bin/env Rscript
# Demo: Converting ECPE wide format data to long format

# Setup
library(dplyr)
library(tidyr)

# Source helper functions
source("../Group9_gen_cdm/00_helpers.r")

cat("\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║           Wide to Long Format Conversion Demo                 ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Example 1: Simple test case
cat("Example 1: Simple test case\n")
cat("════════════════════════════════════════════════════════════════\n\n")

wide_test <- data.frame(
  id = 1:5,
  E1 = c(1, 0, 1, 1, 0),
  E2 = c(1, 1, 0, 1, 1),
  E3 = c(0, 1, 1, 0, 1),
  E4 = c(1, 1, 1, 1, 0)
)

cat("Wide format (5 students × 4 items):\n")
print(wide_test)

# Convert to long using tidyverse
long_test <- wide2long_tidyverse(wide_test)

cat("\nLong format (20 rows = 5 students × 4 items):\n")
print(long_test)

cat("\nStructure:\n")
str(long_test)

# Example 2: ECPE-like data
cat("\n\nExample 2: ECPE-like data structure\n")
cat("════════════════════════════════════════════════════════════════\n\n")

# Create sample ECPE-like data (first 10 students, 28 items)
set.seed(123)
n_students <- 10
n_items <- 28

ecpe_wide <- data.frame(id = 1:n_students)
for (i in 1:n_items) {
  ecpe_wide[[paste0("E", i)]] <- sample(0:1, n_students, replace = TRUE)
}

cat("Wide ECPE data (first 5 students, first 8 items):\n")
print(ecpe_wide[1:5, 1:9])

# Convert to long
ecpe_long <- wide2long_tidyverse(ecpe_wide)

cat("\nLong ECPE data (first 30 rows):\n")
print(head(ecpe_long, 30))

cat("\nDimensions:\n")
cat("  Wide:  ", nrow(ecpe_wide), "rows ×", ncol(ecpe_wide), "columns\n")
cat("  Long:  ", nrow(ecpe_long), "rows ×", ncol(ecpe_long), "columns\n")
cat("  Expected:", n_students * n_items, "rows (", n_students, "students ×",
    n_items, "items)\n")

# Example 3: Round-trip conversion test
cat("\n\nExample 3: Round-trip conversion test\n")
cat("════════════════════════════════════════════════════════════════\n\n")

cat("Testing: wide → long → wide\n\n")

# Start with wide
original_wide <- wide_test

# Convert to long
intermediate_long <- wide2long_tidyverse(original_wide)

# Convert back to wide
reconstructed_wide <- long2wide_tidyverse(intermediate_long)

cat("Original wide:\n")
print(original_wide)

cat("\nReconstructed wide:\n")
print(reconstructed_wide)

cat("\nAre they equal?",
    all.equal(original_wide, as.data.frame(reconstructed_wide)), "\n")

# Example 4: Handling missing values
cat("\n\nExample 4: Handling missing values\n")
cat("════════════════════════════════════════════════════════════════\n\n")

wide_with_na <- data.frame(
  id = 1:4,
  E1 = c(1, NA, 1, 0),
  E2 = c(NA, 1, 0, 1),
  E3 = c(0, 1, NA, 1),
  E4 = c(1, 0, 1, NA)
)

cat("Wide format with missing values:\n")
print(wide_with_na)

long_with_na <- wide2long_tidyverse(wide_with_na)

cat("\nLong format:\n")
print(long_with_na)

cat("\nMissing values are preserved as NA\n")

# Example 5: Using with actual ECPE data (if available)
cat("\n\nExample 5: Using with actual ECPE data\n")
cat("════════════════════════════════════════════════════════════════\n\n")

# Try to load ECPE data if it exists
ecpe_data_path <- "data/ecpe_wide.rda"

if (file.exists(ecpe_data_path)) {
  cat("Loading ECPE data from:", ecpe_data_path, "\n")
  load(ecpe_data_path)

  cat("\nOriginal ECPE wide data dimensions:\n")
  cat("  ", nrow(ecpe_wide_data), "students ×", ncol(ecpe_wide_data) - 1, "items\n")

  # Convert to long
  ecpe_long_data <- wide2long_tidyverse(ecpe_wide_data)

  cat("\nConverted to long format:\n")
  cat("  ", nrow(ecpe_long_data), "rows\n")

  cat("\nFirst 20 rows:\n")
  print(head(ecpe_long_data, 20))

} else {
  cat("ECPE data file not found at:", ecpe_data_path, "\n")
  cat("This example shows how you would use it with your actual data:\n\n")

  cat("# Assuming you have ecpe_data as a tibble/data.frame:\n")
  cat("# ecpe_data has columns: id, E1, E2, ..., E28\n\n")
  cat("ecpe_long <- wide2long_tidyverse(ecpe_data)\n\n")
  cat("# Result will have columns: id, student, item, C\n")
  cat("# With nrow = n_students * n_items\n")
}

# Usage summary
cat("\n\n╔════════════════════════════════════════════════════════════════╗\n")
cat("║                      USAGE SUMMARY                             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("Available functions:\n\n")

cat("1. wide2long_base(wide_df)\n")
cat("   - Uses base R reshape()\n")
cat("   - No package dependencies\n\n")

cat("2. wide2long_tidyverse(wide_df)\n")
cat("   - Uses tidyr::pivot_longer()\n")
cat("   - Requires: dplyr, tidyr\n")
cat("   - Returns tibble\n\n")

cat("Input format:\n")
cat("  - Column 'id': student identifier\n")
cat("  - Columns E1, E2, ..., EN: item responses (0/1 or NA)\n\n")

cat("Output format:\n")
cat("  - Column 'id': student identifier\n")
cat("  - Column 'student': same as id (for compatibility)\n")
cat("  - Column 'item': item number (1, 2, 3, ...)\n")
cat("  - Column 'C': response (0, 1, or NA)\n\n")

cat("Example usage:\n\n")
cat("  # Load your ECPE data\n")
cat("  ecpe_wide <- read.csv('ecpe_data.csv')\n\n")
cat("  # Convert to long format\n")
cat("  ecpe_long <- wide2long_tidyverse(ecpe_wide)\n\n")
cat("  # Use with ggdmc\n")
cat("  # (now ecpe_long has the structure expected by your models)\n\n")

cat("✓ Demo completed successfully\n\n")
