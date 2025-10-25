long2wide_base <- function(df) {
    # keep only needed columns
    d <- df[, c("student", "item", "C")]

    # coerce & clean: id = integer(student), item = "E<integer(item)>", C = integer
    id <- as.integer(as.character(d$student))
    item <- paste0("E", as.integer(as.character(d$item)))
    C <- as.integer(d$C)
    d <- data.frame(id = id, item = item, C = C, stringsAsFactors = FALSE)

    # de-duplicate by (id, item): keep first
    # (order by id and numeric part of item so the "first" is stable)
    ord <- order(d$id, suppressWarnings(as.integer(sub("^E", "", d$item))))
    d <- d[ord, , drop = FALSE]
    keep <- !duplicated(d[c("id", "item")])
    d <- d[keep, , drop = FALSE]

    # ensure rectangular layout: all id × item combos present
    ids <- sort(unique(d$id), na.last = NA)
    items <- paste0("E", sort(unique(as.integer(sub("^E", "", d$item)))))
    d$item <- factor(d$item, levels = items)

    grid <- expand.grid(
        id = ids, item = items,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    out <- merge(grid, d, by = c("id", "item"), all.x = TRUE, sort = FALSE)

    # pivot wide (base reshape produces columns C.E1, C.E2, ...; then rename)
    wide <- reshape(out, idvar = "id", timevar = "item", direction = "wide")
    names(wide) <- sub("^C\\.", "", names(wide))

    # ensure integer storage (except id)
    if (ncol(wide) > 1L) {
        for (j in 2:ncol(wide)) storage.mode(wide[[j]]) <- "integer"
    }

    # order by id
    wide[order(wide$id), , drop = FALSE]
}


# long: columns = student (fct/int/chr), item (fct/int/chr), C (0/1), s (ignore)

long2wide_tidyverse <- function(df) {
    require(dplyr)
    require(tidyr)
    require(stringr)
    wide <- df %>%
        # keep only what we need
        select(student, item, C) %>%
        # make sure ids and items are clean, ordered, and item names look like E1..E28
        mutate(
            id   = as.integer(as.character(student)),
            item = paste0("E", as.integer(as.character(item))),
            C    = as.integer(C)
        ) %>%
        select(id, item, C) %>%
        # if there could be duplicates per (id,item), pick one; change to mean/max if needed
        distinct(id, item, .keep_all = TRUE) %>%
        # ensure a full rectangular layout (so missing item administrations show up as NA)
        complete(id, item) %>%
        # pivot to wide
        pivot_wider(names_from = item, values_from = C, values_fill = NA_integer_) %>%
        arrange(id)

    wide
}


# ==============================================================================
# WIDE TO LONG CONVERSION FUNCTIONS
# ==============================================================================

#' Convert Wide Format to Long Format (Base R)
#'
#' Reverse of long2wide_base. Converts wide format CDM data back to long format.
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
    long_df$student <- long_df$id # Add student column for compatibility

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
#' Reverse of long2wide_tidyverse. Converts wide format CDM data back to long format.
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


compare_profile_probs <- function(mat1, mat2, prefix = "pi_", reorder = TRUE) {
    # Remove prefix from mat1 rownames
    rownames1 <- gsub(paste0("^", prefix), "", rownames(mat1))

    # If reorder is TRUE, reorder mat2 to match mat1's order
    if (reorder) {
        # Find common profiles
        common_profiles <- intersect(rownames1, rownames(mat2))

        # Reorder mat2 to match mat1
        mat2_reordered <- mat2[common_profiles, , drop = FALSE]
        mat1_subset <- mat1[paste0(prefix, common_profiles), , drop = FALSE]
        rownames(mat1_subset) <- common_profiles

        # Calculate differences
        diff_mat <- mat1_subset - mat2_reordered

        # Print comparison
        cat("\n=== Profile Probability Comparison ===\n\n")
        cat("Matrix 1 (", deparse(substitute(mat1)), "):\n", sep = "")
        print(mat1_subset)
        cat("\nMatrix 2 (", deparse(substitute(mat2)), "):\n", sep = "")
        print(mat2_reordered)
        cat("\nDifferences (Matrix 1 - Matrix 2):\n")
        print(round(diff_mat, 4))
        cat("\nAbsolute differences:\n")
        print(round(abs(diff_mat), 4))
        cat("\nMax absolute difference:", round(max(abs(diff_mat)), 4), "\n")
        cat("Mean absolute difference:", round(mean(abs(diff_mat)), 4), "\n")

        # Return invisibly for further analysis
        invisible(list(
            mat1 = mat1_subset,
            mat2 = mat2_reordered,
            diff = diff_mat,
            max_abs_diff = max(abs(diff_mat)),
            mean_abs_diff = mean(abs(diff_mat))
        ))
    } else {
        # Just remove prefix without reordering
        mat1_renamed <- mat1
        rownames(mat1_renamed) <- rownames1

        cat("\n=== Matrices with aligned names ===\n\n")
        cat("Matrix 1:\n")
        print(mat1_renamed)
        cat("\nMatrix 2:\n")
        print(mat2)

        invisible(list(mat1 = mat1_renamed, mat2 = mat2))
    }
}

# Quick comparison function that returns just the differences
quick_compare <- function(mat1, mat2, prefix = "pi_") {
    rownames1 <- gsub(paste0("^", prefix), "", rownames(mat1))
    common_profiles <- intersect(rownames1, rownames(mat2))

    mat2_reordered <- mat2[common_profiles, , drop = FALSE]
    mat1_subset <- mat1[paste0(prefix, common_profiles), , drop = FALSE]
    rownames(mat1_subset) <- common_profiles

    mat1_subset - mat2_reordered
}
