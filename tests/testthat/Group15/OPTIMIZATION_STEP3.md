# Step 3: Eliminate Type Conversions in Likelihood Functions

## Problem

Current implementation converts `arma::vec` → `std::vector<double>` on every likelihood evaluation:

```cpp
// Current bottleneck (likelihood.h:317-321)
double sumloglike(arma::vec parameters, bool debug = false) {
    std::vector<double> parameters_std =
        arma::conv_to<std::vector<double>>::from(parameters);  // EXPENSIVE!
    // ... then calls model-specific likelihood
}
```

**Impact**: This conversion happens hundreds/thousands of times per MCMC run, dominating runtime (96.4% "Other" time in benchmarks).

## Solution

Create optimized likelihood path that accepts `const arma::vec&` directly:

1. **Add `sumloglike_fast(const arma::vec& parameters)`** - No conversion, pass by const reference
2. **Modify model likelihood functions** to extract parameters using Armadillo's efficient `.subvec()`
3. **Use Armadillo operations throughout** - Avoid std::vector entirely in hot path

## Implementation Details

### Files to Modify

**ggdmcHeaders/inst/include/ggdmcHeaders/likelihood.h:**
- Add `sumloglike_fast(const arma::vec&)` method
- Add `lba_likelihood_fast()`, `ddm_likelihood_fast()`, `cdm_likelihood_fast()`
- Keep original versions unchanged

**ggdmcHeaders/inst/include/ggdmcHeaders/lba.h:**
- Add `set_parameters_fast(const arma::mat& param_matrix)` accepting Armadillo directly
- Use `.col()` and `.row()` for parameter extraction

**ggdmcHeaders/inst/include/ggdmcHeaders/ddm.h:**
- Add `set_parameters_fast(const arma::vec& params, bool is_minus)`
- Use direct indexing with `params(idx)` instead of vector conversion

**ggdmcHeaders/inst/include/ggdmcHeaders/cdm.h:**
- Similar optimizations for CDM parameter extraction

**src/de.cpp:**
- Add `crossover_fast()`, `migration_fast()`, `run_chains_fast()` using new likelihood path
- Add timing and no-timing versions

**src/de2R.cpp:**
- Export `run_subject_fast()`, `run_subject_fast_notiming()` to R

## Expected Performance

**Baseline**: Step 2 benchmark showed:
- Likelihood calculations: 3.6% of C++ time
- Other (conversions, parameter extraction): 96.4% of C++ time

**Target**: 20-40% overall speedup by eliminating conversion overhead

**Benefits**:
- Universal improvement (all models: LBA, DDM, CDM)
- Benefits both subject-level and hierarchical models
- No algorithmic changes - purely eliminates overhead

## Benchmark Strategy

Same three-version approach:
1. **Original** - Baseline with conversions
2. **Fast with timing** - Optimized + instrumentation
3. **Fast no timing** - Pure optimization

Compare on representative model to isolate pure optimization benefit from measurement overhead.

## Technical Notes

### Armadillo Efficiency

Using Armadillo's `.subvec()`, `.col()`, `.row()` operations:
- **Zero-copy views** when possible
- Efficient memory access patterns
- Already optimized for numerical computation

### Parameter Extraction Pattern

**Before (with conversion)**:
```cpp
std::vector<double> params_std = conv_to<std::vector<double>>::from(params);
// Then extract into std::vector<double> for each parameter
m_A = parameter_matrix[0];  // Copy entire vector
```

**After (direct Armadillo)**:
```cpp
// Use const reference, no conversion
const arma::mat& param_matrix = ...;
m_A = arma::conv_to<std::vector<double>>::from(param_matrix.row(0).t());
// Only convert what's needed, or keep as arma::vec
```

**Even Better**:
```cpp
// Keep as Armadillo vectors throughout
arma::vec m_A_vec = param_matrix.row(0).t();
// Use arma::vec operations directly
```

## Risk Assessment

**Low Risk**:
- Keep original functions intact
- New functions are parallel implementations
- Can easily revert if no benefit

**Testing**:
- Benchmark will show if optimization is effective
- Numerical results should be identical (same algorithm)
- Use existing test suite to verify correctness

## Next Steps

1. Implement optimized likelihood path in ggdmcHeaders
2. Add fast DE functions in src/de.cpp
3. Export to R via src/de2R.cpp
4. Create benchmark script
5. Test and measure performance gain

---
Created: 2025-10-15
Status: In Progress
