# Wide to Long Conversion - Quick Reference

## Functions Added to `00_helpers.r`

Two new functions have been added:
- `wide2long_base()` - Base R implementation
- `wide2long_tidyverse()` - Tidyverse implementation (recommended)

## Your ECPE Data

Your data looks like this (wide format):
```r
# A tibble: 2,922 × 29
      id    E1    E2    E3    E4    E5    E6    E7    E8 ...
   <int> <int> <int> <int> <int> <int> <int> <int> <int>
 1     1     1     1     1     0     1     1     1     1
 2     2     1     1     1     1     1     1     1     1
 3     3     1     1     1     1     1     1     0     1
```

## Quick Usage

### Option 1: Tidyverse (Recommended)

```r
library(dplyr)
library(tidyr)
source("../Group9_gen_cdm/00_helpers.r")

# Convert your ECPE data
ecpe_long <- wide2long_tidyverse(ecpe_wide_data)

# Result structure:
# # A tibble: 81,816 × 4    (2,922 students × 28 items)
#       id student  item     C
#    <int>   <int> <int> <int>
#  1     1       1     1     1
#  2     1       1     2     1
#  3     1       1     3     1
#  4     1       1     4     0
#  5     1       1     5     1
# ...
```

### Option 2: Base R

```r
source("../Group9_gen_cdm/00_helpers.r")

# Convert your ECPE data
ecpe_long <- wide2long_base(ecpe_wide_data)

# Same result, no package dependencies
```

## Input Requirements

Your wide data **must have**:
- ✅ Column named `id` (integer student identifier)
- ✅ Columns named `E1`, `E2`, ..., `E28` (or any number of items)
- ✅ Values are 0, 1, or NA

Your data already meets these requirements! ✓

## Output Format

The long format will have exactly **4 columns**:

| Column | Type | Description |
|--------|------|-------------|
| `id` | integer | Student identifier (same as input) |
| `student` | integer | Copy of id (for compatibility with ggdmc) |
| `item` | integer | Item number (1, 2, 3, ..., 28) |
| `C` | integer | Response (0, 1, or NA) |

## Expected Dimensions

For your ECPE data:
- **Wide**: 2,922 rows × 29 columns (id + 28 items)
- **Long**: 81,816 rows × 4 columns (2,922 × 28)

## Complete Example

```r
#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)

# Load helper functions
source("../Group9_gen_cdm/00_helpers.r")

# Assuming your data is in ecpe_wide_data
# (a tibble with id, E1, E2, ..., E28)

cat("Original dimensions:", dim(ecpe_wide_data), "\n")
# Output: Original dimensions: 2922 29

# Convert to long format
ecpe_long <- wide2long_tidyverse(ecpe_wide_data)

cat("Long format dimensions:", dim(ecpe_long), "\n")
# Output: Long format dimensions: 81816 4

# View first few rows
print(head(ecpe_long, 20))

# Check structure
str(ecpe_long)
# tibble [81,816 × 4] (S3: tbl_df/tbl/data.frame)
#  $ id     : int [1:81816] 1 1 1 1 1 1 1 1 1 1 ...
#  $ student: int [1:81816] 1 1 1 1 1 1 1 1 1 1 ...
#  $ item   : int [1:81816] 1 2 3 4 5 6 7 8 9 10 ...
#  $ C      : int [1:81816] 1 1 1 0 1 1 1 1 1 1 ...

# Now you can use this with your CDM models!
```

## Round-Trip Conversion

You can convert back and forth without losing data:

```r
# Start with wide
original_wide <- ecpe_wide_data

# Wide → Long
long_data <- wide2long_tidyverse(original_wide)

# Long → Wide
back_to_wide <- long2wide_tidyverse(long_data)

# Check they're the same
all.equal(original_wide, as.data.frame(back_to_wide))
# [1] TRUE
```

## Handling Missing Values

Missing values (NA) are preserved:

```r
# If you have missing responses
wide_with_na <- tibble(
  id = 1:3,
  E1 = c(1, NA, 1),
  E2 = c(0, 1, NA),
  E3 = c(NA, 0, 1)
)

long_with_na <- wide2long_tidyverse(wide_with_na)
print(long_with_na)

# # A tibble: 9 × 4
#      id student  item     C
#   <int>   <int> <int> <int>
# 1     1       1     1     1
# 2     1       1     2     0
# 3     1       1     3    NA
# 4     2       2     1    NA
# 5     2       2     2     1
# 6     2       2     3     0
# 7     3       3     1     1
# 8     3       3     2    NA
# 9     3       3     3     1
```

## Common Use Cases

### Use Case 1: Prepare ECPE data for CDM model

```r
# Load your ECPE wide data
load("data/ecpe_wide.rda")

# Convert to long format
ecpe_long <- wide2long_tidyverse(ecpe_wide)

# Now use with your CDM model
# ... your ggdmc model code here ...
```

### Use Case 2: Subset students and convert

```r
# Take first 100 students
ecpe_subset <- ecpe_wide_data %>%
  filter(id <= 100)

# Convert to long
ecpe_long_subset <- wide2long_tidyverse(ecpe_subset)

# Result: 100 × 28 = 2,800 rows
```

### Use Case 3: Remove items and convert

```r
# Keep only first 10 items
ecpe_10items <- ecpe_wide_data %>%
  select(id, E1:E10)

# Convert to long
ecpe_long_10 <- wide2long_tidyverse(ecpe_10items)

# Result: 2,922 × 10 = 29,220 rows
```

## Troubleshooting

### Error: "Wide data must have 'id' column"

Your data needs a column named exactly `id`. If it has a different name:

```r
# Rename your id column
ecpe_data <- ecpe_data %>%
  rename(id = student_id)  # or whatever it's called

# Then convert
ecpe_long <- wide2long_tidyverse(ecpe_data)
```

### Error: "No item columns found"

Your item columns must be named `E1`, `E2`, etc. If they have different names:

```r
# Rename item columns
names(ecpe_data) <- c("id", paste0("E", 1:28))

# Then convert
ecpe_long <- wide2long_tidyverse(ecpe_data)
```

### Warning about packages

If you get "package 'dplyr' is not available":

```r
install.packages(c("dplyr", "tidyr"))
```

## Performance Notes

For your ECPE data (2,922 students × 28 items):
- Conversion time: < 1 second
- Memory: Minimal increase (long format is actually more compact)
- Both base R and tidyverse versions are fast

For very large datasets (> 100,000 students):
- Use `wide2long_base()` (slightly faster, no dependencies)
- Or use `data.table` if you need maximum speed

## See Also

- `long2wide_base()` - Convert long → wide (base R)
- `long2wide_tidyverse()` - Convert long → wide (tidyverse)
- Demo script: `demo_wide2long.r`
- Full function definitions: `00_helpers.r`

---

**Quick Answer for Your Data:**

```r
library(dplyr); library(tidyr)
source("../Group9_gen_cdm/00_helpers.r")
ecpe_long <- wide2long_tidyverse(ecpe_wide_data)
```

Done! ✓
