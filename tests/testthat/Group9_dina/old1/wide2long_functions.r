#!/usr/bin/env Rscript
# Functions to convert wide format CDM data to long format
# Companion to long2wide functions in 00_helpers.r

#' Convert Wide Format to Long Format (Base R)
#'
#' @param wide_df Data frame in wide format with columns:
#'   - id: subject/student identifier
#'   - E1, E2, ..., EN: item response columns (0/1)
#'
#' @return Data frame in long format with columns:
#'   - id: subject identifier
#'   - student: same as id (for compatibility)
#'   - item: item number (1, 2, 3, ...)
#'   - C: response (0/1 or NA)
#'
#' @examples
#' \dontrun{
#' # Example wide data
#' wide_df <- data.frame(
#'   id = 1:3,
#'   E1 = c(1, 0, 1),
#'   E2 = c(1, 1, 0),
#'   E3 = c(0, 1, 1)
#' )
#' long_df <- wide2long_base(wide_df)
#' print(long_df)
#' }
wide2long_base <- function(wide_df) {
  # Extract id column
  if (!"id" %in% names(wide_df)) {
    stop("Wide data must have 'id' column")
  }

  # Identify item columns (those starting with E followed by digits)
  item_cols <- grep("^E[0-9]+$", names(wide_df), value = TRUE)

  if (length(item_cols) == 0) {
    stop("No item columns found (expected E1, E2, etc.)")
  }

  # Extract item numbers from column names
  item_numbers <- as.integer(sub("^E", "", item_cols))

  # Use reshape to convert from wide to long
  long_df <- reshape(
    wide_df,
    varying = item_cols,
    v.names = "C",
    timevar = "item",
    times = item_numbers,
    idvar = "id",
    direction = "long"
  )

  # Clean up
  rownames(long_df) <- NULL
  long_df$student <- long_df$id  # Add student column for compatibility

  # Reorder columns: id, student, item, C
  long_df <- long_df[, c("id", "student", "item", "C")]

  # Sort by id, then item
  long_df <- long_df[order(long_df$id, long_df$item), , drop = FALSE]

  # Ensure proper types
  long_df$id <- as.integer(long_df$id)
  long_df$student <- as.integer(long_df$student)
  long_df$item <- as.integer(long_df$item)
  long_df$C <- as.integer(long_df$C)

  return(long_df)
}


#' Convert Wide Format to Long Format (Tidyverse)
#'
#' @param wide_df Data frame or tibble in wide format with columns:
#'   - id: subject/student identifier
#'   - E1, E2, ..., EN: item response columns (0/1)
#'
#' @return Tibble in long format with columns:
#'   - id: subject identifier
#'   - student: same as id (for compatibility)
#'   - item: item number (1, 2, 3, ...)
#'   - C: response (0/1 or NA)
#'
#' @examples
#' \dontrun{
#' library(dplyr)
#' library(tidyr)
#'
#' # Example wide data
#' wide_df <- tibble(
#'   id = 1:3,
#'   E1 = c(1, 0, 1),
#'   E2 = c(1, 1, 0),
#'   E3 = c(0, 1, 1)
#' )
#' long_df <- wide2long_tidyverse(wide_df)
#' print(long_df)
#' }
wide2long_tidyverse <- function(wide_df) {
  require(dplyr)
  require(tidyr)

  # Check for id column
  if (!"id" %in% names(wide_df)) {
    stop("Wide data must have 'id' column")
  }

  long_df <- wide_df %>%
    # Pivot all E* columns to long format
    pivot_longer(
      cols = matches("^E[0-9]+$"),
      names_to = "item",
      values_to = "C",
      names_prefix = "E",
      names_transform = list(item = as.integer)
    ) %>%
    # Add student column (same as id for compatibility)
    mutate(
      student = as.integer(id),
      item = as.integer(item),
      C = as.integer(C)
    ) %>%
    # Reorder columns
    select(id, student, item, C) %>%
    # Sort by id and item
    arrange(id, item)

  return(long_df)
}


#' Convert Wide Format to Long Format (Tidyverse, keep item names)
#'
#' Alternative version that keeps item names as "E1", "E2", etc.
#' instead of converting to numeric.
#'
#' @param wide_df Data frame or tibble in wide format
#'
#' @return Tibble in long format with item as character (E1, E2, etc.)
wide2long_tidyverse_named <- function(wide_df) {
  require(dplyr)
  require(tidyr)

  if (!"id" %in% names(wide_df)) {
    stop("Wide data must have 'id' column")
  }

  long_df <- wide_df %>%
    pivot_longer(
      cols = matches("^E[0-9]+$"),
      names_to = "item",
      values_to = "C"
    ) %>%
    mutate(
      student = as.integer(id),
      C = as.integer(C)
    ) %>%
    select(id, student, item, C) %>%
    arrange(id, item)

  return(long_df)
}


# ==============================================================================
# DEMONSTRATION AND TESTING
# ==============================================================================

if (FALSE) {
  # Example 1: Small test case
  cat("\n=== Example 1: Small test case ===\n")

  wide_test <- data.frame(
    id = 1:5,
    E1 = c(1, 0, 1, 1, 0),
    E2 = c(1, 1, 0, 1, 1),
    E3 = c(0, 1, 1, 0, 1),
    E4 = c(1, 1, 1, 1, 0)
  )

  cat("\nWide format:\n")
  print(wide_test)

  # Base R version
  long_base <- wide2long_base(wide_test)
  cat("\nLong format (base R):\n")
  print(head(long_base, 12))

  # Tidyverse version
  library(dplyr)
  library(tidyr)
  long_tidy <- wide2long_tidyverse(wide_test)
  cat("\nLong format (tidyverse):\n")
  print(head(long_tidy, 12))

  # Test round-trip conversion
  cat("\n=== Testing round-trip conversion ===\n")
  source("00_helpers.r")

  # wide -> long -> wide
  wide_reconstructed <- long2wide_base(long_base)
  cat("\nRound-trip successful:", all.equal(wide_test, wide_reconstructed), "\n")

  # Example 2: With missing values
  cat("\n=== Example 2: With missing values ===\n")

  wide_with_na <- data.frame(
    id = 1:3,
    E1 = c(1, NA, 1),
    E2 = c(NA, 1, 0),
    E3 = c(0, 1, NA)
  )

  cat("\nWide format with NAs:\n")
  print(wide_with_na)

  long_with_na <- wide2long_tidyverse(wide_with_na)
  cat("\nLong format:\n")
  print(long_with_na)

  # Example 3: Your ECPE data structure
  cat("\n=== Example 3: ECPE-like data ===\n")

  # Simulate your data structure
  set.seed(123)
  n_students <- 10
  n_items <- 28

  ecpe_wide <- data.frame(
    id = 1:n_students
  )

  for (i in 1:n_items) {
    ecpe_wide[[paste0("E", i)]] <- sample(0:1, n_students, replace = TRUE)
  }

  cat("\nWide ECPE data (first 5 students, first 8 items):\n")
  print(ecpe_wide[1:5, 1:9])

  ecpe_long <- wide2long_tidyverse(ecpe_wide)
  cat("\nLong ECPE data (first 20 rows):\n")
  print(head(ecpe_long, 20))

  cat("\nDimensions:\n")
  cat("  Wide:", nrow(ecpe_wide), "×", ncol(ecpe_wide), "\n")
  cat("  Long:", nrow(ecpe_long), "×", ncol(ecpe_long), "\n")
  cat("  Expected long rows:", n_students * n_items, "\n")
}
