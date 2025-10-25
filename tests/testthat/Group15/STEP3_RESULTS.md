# Step 3: Type Conversion Optimization - Results and Analysis

## Date: 2025-10-15

## Implementation Summary

Attempted to optimize type conversions in likelihood evaluation by:
1. Creating `sumloglike_fast(const arma::vec&)` - accepts const reference instead of by-value
2. Creating `lba_likelihood_fast()` - optimized likelihood path
3. Adding `run_subject_typeconv()` wrapper

## Test Results

**Performance:**
- Original version: 0.026s (200 iterations)
- Typeconv version: 0.056s (200 iterations)
- **Result: 115% SLOWER** ❌

**Correctness:**
- Dimensions match: ✓
- Log-prior difference: 0.0 (identical)
- Log-likelihood difference: 2205.14 (DIFFERENT!) ⚠️

## Root Cause Analysis

### Why It's Slower

The optimization didn't actually eliminate conversions:

**Original code path:**
```cpp
sumloglike(arma::vec parameters) {          // Pass by value = 1 copy
    vector<double> params = convert(parameters);  // 1 conversion
    lba_likelihood(params);                       // No conversion
}
```

**Optimized code path:**
```cpp
sumloglike_fast(const arma::vec& parameters) {    // No copy ✓
    lba_likelihood_fast(parameters);               // Const ref ✓
        vector<double> params = convert(parameters); // Still 1 conversion!
        // Same downstream processing
}
```

**Net result:**
- Saved 1 arma::vec copy
- Still do 1 conversion
- Added extra function call overhead
- **No real benefit, just complexity**

### Why Results Differ

The large log-likelihood difference (2205.14) is concerning and suggests:
1. Possible RNG differences due to different code paths
2. OR a subtle bug in the implementation
3. Need to investigate further before deploying

## Why This Optimization Failed

### Misconception in Original Analysis

The Step 2 benchmark showed "96.4% Other time", which I attributed to type conversions. However:

1. **The original code is already efficient**: It converts once per likelihood call
2. **The real "Other" time** likely includes:
   - Parameter extraction and mapping
   - Model setup overhead
   - Armadillo matrix operations
   - Cache misses
   - NOT primarily type conversions

### Fundamental Issue

To truly eliminate conversions, we would need to:
1. Modify LBA/DDM/CDM classes to accept Armadillo types natively
2. Redesign the entire parameter extraction pipeline
3. This is a MAJOR refactoring, not a simple optimization

### Comparison to Previous Steps

| Step | Target | Result | Why |
|------|--------|--------|-----|
| 1 | Pre-allocate buffers | 2% slower | arma::shuffle already optimal |
| 2 | Cache HB likelihoods | 12% slower | Likelihood only 3.6% of time |
| 3 | Eliminate conversions | 115% slower | Conversion not the bottleneck |

**Pattern:** All three optimizations targeted the wrong bottleneck!

## Lessons Learned

### 1. Profile First, Optimize Second

The "96.4% Other" category was too vague. Should have used:
- Detailed C++ profilers (perf, gprof, valgrind)
- Measured actual conversion overhead
- Identified specific hot functions

### 2. The Original Code is Well-Optimized

- Uses Armadillo efficiently where it matters
- Conversions are minimal and necessary
- The DE-MCMC algorithm itself may be the bottleneck, not implementation details

### 3. Premature Optimization

All three steps were premature optimizations based on assumptions rather than hard profiling data.

## Recommendations

###  1. Accept Current Performance

The package is already quite fast for its problem domain. Focus on:
- Correctness
- Maintainability
- User features

### 2. If Performance is Critical

Use proper profiling tools:
```bash
# Linux: perf
R -d "perf record -g" --vanilla < benchmark.R
perf report

# Or valgrind callgrind
R -d "valgrind --tool=callgrind" --vanilla < benchmark.R
kcachegrind callgrind.out.*
```

### 3. Algorithmic Improvements

Instead of micro-optimizations, consider:
- Adaptive MCMC (fewer iterations to convergence)
- Parallel chains (already implemented with ncore > 1)
- GPU acceleration for likelihood calculations (major undertaking)

### 4. Revert Step 3 Changes

**Recommendation: Do NOT deploy the type conversion optimization**

Reasons:
- Slower performance
- Different numerical results (concerning)
- Added code complexity
- No measurable benefit

## What to Keep

From this exercise, the valuable outputs are:
1. ✅ Documentation of what was tried and why it failed
2. ✅ Better understanding of the codebase structure
3. ✅ Test infrastructure for benchmarking
4. ✅ Confirmation that original code is well-optimized

## Code Status

### Files Modified (Should Revert)

**ggdmcHeaders:**
- `inst/include/ggdmcHeaders/design.h` - Added `set_parameter_values_fast()`
- `inst/include/ggdmcHeaders/lba.h` - Added `set_parameters_fast()`
- `inst/include/ggdmcHeaders/likelihood.h` - Added `lba_likelihood_fast()`, `sumloglike_fast()`

**ggdmc:**
- `src/de.h` - Added typeconv function declarations
- `src/de.cpp` - Added `crossover_typeconv()`, etc.
- `src/de2R.cpp` - Added `run_subject_typeconv()`
- `NAMESPACE` - Exported `run_subject_typeconv`

### Files to Keep

**Documentation:**
- `OPTIMIZATION_STEP3.md` - Design document
- `STEP3_IMPLEMENTATION_SUMMARY.md` - What was built
- `STEP3_RESULTS.md` - This file

**Tests:**
- `tests/testthat/Group15/02_test_typeconv.r` - Benchmark script

## Conclusion

**Step 3 optimization should be abandoned.** The implementation:
- Makes code slower
- Produces different results
- Adds complexity without benefit
- Targets wrong bottleneck

The original ggdmc code is already well-optimized. Future performance work should:
1. Use proper profiling tools first
2. Focus on algorithmic improvements
3. Avoid micro-optimizations without evidence

---

**Status: FAILED - Do Not Deploy**
**Recommendation: Revert all Step 3 changes**
**Lesson: Always profile before optimizing**
