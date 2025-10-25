# Asymptotic Behavior Study - Ready to Run

## Status: ✓ Complete Setup

All scripts and documentation are in place for studying asymptotic behavior of DINA model parameter estimates with discrete profile probabilities.

## Quick Start

```bash
cd /media/yslin/Tui/01_Projects/ggdmc/tests/testthat/Group9_gen_cdm

# Step 1: Generate 10 datasets (N = 100 to 100,000)
# Time: ~2-5 minutes
Rscript 01_gen_subject_dina0.r

# Step 2: Fit MCMC and analyze asymptotic behavior
# Time: ~30-60 minutes (500 iterations × 3 chains × 10 datasets)
Rscript 02_analyze_asymptotic.r

# View results
ls -lh data/      # Generated datasets
ls -lh results/   # Analysis results and plots
```

## What Gets Generated

### Step 1 Output (Data Generation)
- **10 datasets**: `data/subject_dina0_N000100.rda` through `data/subject_dina0_N100000.rda`
- **Metadata**: `data/generation_metadata.rda`
- **Console output**: Empirical profile frequencies for each N

### Step 2 Output (Analysis)
- **Results**: `results/asymptotic_results.rda` (all estimates and statistics)
- **4 diagnostic plots**:
  1. `bias_vs_N.pdf` - Bias should converge to 0
  2. `se_vs_N.pdf` - SE should decay as O(1/√N)
  3. `coverage_vs_N.pdf` - 95% CI coverage should approach 0.95
  4. `estimates_vs_N.pdf` - Estimates should converge to true values
- **Console output**: Summary statistics and asymptotic normality checks

## Model Specification

### Q-Matrix (J=3 items, K=2 skills)
```
       Algebra  Geometry
Item1    1        0
Item2    0        1
Item3    1        1
```

### True Parameter Values (9 total)
```r
Guess:  guess1 = 0.10, guess2 = 0.20, guess3 = 0.30

Profile probabilities (compact, L-1 = 3):
  pi1 = 0.1578  # P(no skills)
  pi2 = 0.2630  # P(Algebra only)
  pi3 = 0.1508  # P(Geometry only)
  # pi4 = 0.4284  # P(both skills) = 1 - pi1 - pi2 - pi3

Slip:  slip1 = 0.20, slip2 = 0.40, slip3 = 0.60
```

### Parameter Layout
```
[guess1, guess2, guess3, pi1, pi2, pi3, slip1, slip2, slip3]
 |<------- J=3 -------->| |<- L-1=3 ->| |<------- J=3 -------->|
```

Total: 2*J + (L-1) = 2*3 + 3 = 9 parameters

## Sample Sizes
```r
N = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000]
```

## Expected Asymptotic Behavior

As N → ∞:

1. **Consistency**: Estimates → True values
   - Bias → 0

2. **Efficiency**: Standard errors → 0
   - SE decreases as O(1/√N)
   - On log-log plot: slope ≈ -0.5

3. **Asymptotic Normality**: √N × (estimate - true) → N(0, σ²)
   - Standardized bias should be O(1), not growing with N

4. **Coverage**: 95% CI coverage → 0.95
   - May be poor for small N
   - Should stabilize around 0.95 for large N

## Adjusting MCMC Settings

For quick testing, reduce iterations in `02_analyze_asymptotic.r` (lines 55-57):

```r
nmc <- 200      # Faster: 200 instead of 500
nchain <- 2     # Faster: 2 instead of 3
burnin <- 0.5   # Keep at 50%
```

For production runs, increase for better convergence:

```r
nmc <- 1000     # Better convergence
nchain <- 4     # More robust
burnin <- 0.5   # May need to increase if chains mix slowly
```

## Important Note: Parameter Count Fix

**Before running this study**, verify that the parameter count issue is fixed in `cdm.h`.

The `set_use_mvn()` function should recalculate `m_expected_nparameter`:

```cpp
void set_use_mvn(bool use_mvn)
{
    m_use_mvn = use_mvn;

    // Recalculate expected parameter count for DINA/DINO
    if (m_rule == Rule::DINA || m_rule == Rule::DINO)
    {
        if (m_use_mvn)
            m_expected_nparameter = 2u * m_n_item + m_n_skill + 1;
        else
            m_expected_nparameter = 2u * m_n_item + (m_n_profile - 1);
    }
}
```

**Quick test**:
```bash
Rscript debug_parameter_count.r
```

Expected output: No `[expected, input] = 6, 9` error

## Troubleshooting

### If data generation fails
- Check that all packages are installed: `ggdmc`, `ggdmcPrior`, `ggdmcModel`, `cdModel`
- Verify `cdm.h` has the discrete profile probability changes (see `DINA_DINO_DISCRETE_PROFILES.md`)
- Check parameter count fix is applied (see `FIX_PARAMETER_COUNT_ISSUE.md`)

### If MCMC doesn't converge
- For small N (100-500), convergence may be poor - this is expected
- Check Gelman-Rubin diagnostic (Rhat) in output
- Increase `nmc` or `burnin` if needed

### If memory issues occur
- N=100,000 uses ~500MB RAM per chain
- Reduce number of chains or run on subset of N values first

## Documentation

- **Full workflow guide**: `README_ASYMPTOTIC.md`
- **Implementation details**: `/media/yslin/Tui/01_Projects/ggdmcHeaders/DINA_DINO_DISCRETE_PROFILES.md`
- **Parameter count fix**: `/media/yslin/Tui/01_Projects/ggdmcHeaders/FIX_PARAMETER_COUNT_ISSUE.md`

## Next Steps After Running

1. Examine the 4 diagnostic plots to verify asymptotic properties
2. Check console output for summary statistics
3. Verify that:
   - ✓ Bias → 0 as N increases
   - ✓ SE decreases as 1/√N
   - ✓ Coverage rate → 0.95 for large N
   - ✓ No systematic bias patterns

If all checks pass, the DINA model with discrete profile probabilities is working correctly!
