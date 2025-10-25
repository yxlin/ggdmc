# DE Optimizer Optimization - Step 1: Pre-allocated Chain Indices

## Summary

This document describes the first optimization step for the ggdmc package's Differential Evolution (DE) optimizer. This optimization focuses on eliminating repeated memory allocations and vector shuffling operations in the hot path of the MCMC sampler.

## Problem Identified

**Original Issue (High Impact, Low Complexity):**
- `get_chains()` and `get_subchains()` create full copies of m_chains and shuffle them on every call
- This happens in every crossover/migration step of the MCMC iteration
- **Impact**: O(n_chain) memory allocation + O(n_chain log n_chain) shuffle per MCMC iteration

## Solution Implemented

### 1. Pre-allocated Buffers (de.h)
Added three member variables to the `de_class`:
```cpp
arma::uvec m_chains_buffer;     // Pre-allocated buffer for chain shuffling
arma::uvec m_subchains_buffer;  // Pre-allocated buffer for subchain selection
arma::uvec m_shuffle_indices;   // Pre-allocated shuffle indices
```

These buffers are initialized once in the constructor and reused throughout execution.

### 2. Fast Chain Selection Functions (de.cpp)

#### `get_chains_fast(unsigned int k, unsigned int nsubchain)`
- Updates `m_subchains` in-place instead of returning a new vector
- Uses Fisher-Yates shuffle algorithm without memory allocation
- Excludes chain k and selects nsubchain chains efficiently

#### `get_subchains_fast()`
- Updates `m_subchains` in-place using pre-allocated buffers
- Implements Fisher-Yates shuffle without temporary allocations
- Maintains same behavior as original but with better memory efficiency

### 3. Optimized Sampler Functions (de.cpp)

#### `crossover_fast()`
- Uses `get_chains_fast()` instead of `get_chains()`
- Includes timing instrumentation for performance analysis
- Tracks crossover time and likelihood computation time separately

#### `migration_fast()`
- Uses `get_subchains_fast()` instead of `get_subchains()`
- Includes timing instrumentation
- Tracks migration time and likelihood computation time

#### `run_chains_fast()`
- Coordinates the optimized crossover and migration operations
- Collects comprehensive timing statistics
- Resets timing stats at the start of each run

### 4. Timing Instrumentation (de.h)

Added `TimingStats` structure to track:
```cpp
struct TimingStats {
    double total_time = 0.0;        // Total execution time
    double crossover_time = 0.0;     // Time in crossover operations
    double migration_time = 0.0;     // Time in migration operations
    double likelihood_time = 0.0;    // Time computing likelihoods
    unsigned int n_crossover = 0;    // Number of crossover calls
    unsigned int n_migration = 0;    // Number of migration calls
};
```

This allows precise measurement of where time is spent in the C++ code.

### 5. R Interface (de2R.cpp)

Created two new functions that wrap the fast implementations:

#### `run_subject_fast()`
- Parallel to `run_subject()` but uses `run_chains_fast()`
- Returns a list with:
  - `$posterior`: Same structure as original function
  - `$timing`: C++ timing statistics

#### `run_hyper_fast()`
- Parallel to `run_hyper()` but uses `run_chains_fast()`
- Returns a list with posterior and timing information

#### Helper function `get_timing_stats()`
- Converts C++ TimingStats to R list
- Makes timing information accessible from R

## Files Modified

1. **src/de.h**
   - Added `<chrono>` header for timing
   - Added pre-allocated buffer member variables
   - Added declarations for fast functions
   - Added `TimingStats` structure

2. **src/de.cpp**
   - Implemented `get_chains_fast()` (lines 92-116)
   - Implemented `get_subchains_fast()` (lines 118-143)
   - Implemented `crossover_fast()` (lines 310-366)
   - Implemented `migration_fast()` (lines 368-422)
   - Implemented `run_chains_fast()` (lines 424-474)
   - Buffer initialization in constructors (lines 15-18, 34-37)

3. **src/de2R.cpp**
   - Added `get_timing_stats()` helper (lines 5-14)
   - Implemented `run_subject_fast()` (lines 121-141)
   - Implemented `run_hyper_fast()` (lines 146-168)

4. **R/RcppExports.R** (auto-generated)
   - Exports `run_subject_fast()` and `run_hyper_fast()` to R

## Key Design Decisions

### 1. Parallel Implementation (Not Replacement)
- Original functions remain unchanged for backward compatibility
- New `_fast` versions run alongside original versions
- Allows easy benchmarking and validation

### 2. In-place Updates
- Fast functions update `m_subchains` member variable directly
- Eliminates return value overhead
- Reduces memory allocations

### 3. Fisher-Yates Shuffle
- Implemented in-place shuffle algorithm
- Same statistical properties as `arma::shuffle()`
- Zero additional memory allocation

### 4. Comprehensive Timing
- Measures time at C++ level (more accurate than R timing)
- Separates different operation types
- Helps identify remaining bottlenecks

## Testing and Validation

### Compilation
Package compiles successfully without warnings or errors:
```bash
R CMD INSTALL --preclean --no-multiarch --with-keep.source .
# Result: SUCCESS
```

### Benchmarking
Two benchmark scripts created:

1. **benchmark_fast_version.R** - Full comparison template
   - Compares original vs fast version
   - Reports speedup and time breakdown
   - Shows C++ vs R overhead

2. **benchmark_simple.R** - Usage guide
   - Documents how to use new functions
   - Explains timing statistics
   - Provides example code

## Expected Performance Impact

### Theoretical Gains
- **Memory allocations**: Reduced from 2-3 per MCMC iteration to 0
- **Vector copies**: Eliminated repeated copying of chain indices
- **Cache efficiency**: Better memory locality with pre-allocated buffers

### Measurement Strategy
1. **R-level timing** (system.time): Total elapsed time including overhead
2. **C++-level timing**: Precise measurement within critical sections
3. **Component breakdown**: Separate crossover, migration, likelihood times

This allows identification of:
- Time spent in R vs C++
- Time spent in different MCMC operations
- Impact of optimization on overall runtime

## Next Steps

After benchmarking this optimization:

1. **If successful** (>10% speedup): Proceed to step 2
   - Cache hyper-likelihood calculations (HB models)
   - Eliminate redundant sumloghlike() calls

2. **If minimal impact**: Investigate other bottlenecks
   - Likelihood function optimization may have bigger impact
   - Profile to find actual hot spots

3. **Documentation**: Update user-facing documentation if fast versions become default

## Usage Example

```r
library(ggdmc)

# Load your model (config, dmi, samples)
source("tests/testthat/Group1/data/0_lba_model0.r")

# Original version
result_old <- ggdmc:::run_subject(config, dmi, samples)

# Fast version with timing
result_fast <- ggdmc:::run_subject_fast(config, dmi, samples)

# Access results
posterior <- result_fast$posterior  # Same as result_old
timing <- result_fast$timing

# Print timing breakdown
cat("Total time:", timing$total_time, "seconds\n")
cat("Crossover:", timing$crossover_time, "seconds\n")
cat("Migration:", timing$migration_time, "seconds\n")
cat("Likelihood:", timing$likelihood_time, "seconds\n")
```

## Implementation Quality

- **Correctness**: Uses same algorithms as original, just more efficiently
- **Maintainability**: Clean separation between old and new code
- **Readability**: Well-commented, follows existing code style
- **Safety**: No changes to memory management beyond pre-allocation
- **Testing**: Can be validated by comparing results with original version

## Conclusion

Step 1 optimization is complete and ready for benchmarking. The implementation:
- ✅ Eliminates repeated memory allocations
- ✅ Maintains backward compatibility
- ✅ Includes comprehensive timing instrumentation
- ✅ Compiles without errors
- ✅ Ready for performance testing

The timing instrumentation will help assess the impact of this optimization and guide future optimization efforts.
