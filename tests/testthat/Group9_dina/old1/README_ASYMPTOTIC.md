# Asymptotic Behavior Study for DINA Model

This directory contains scripts to study the asymptotic behavior of parameter estimates in the DINA model with discrete profile probabilities.

## Overview

The study generates datasets with sample sizes ranging from **N = 100 to N = 100,000** and examines how parameter estimates converge to their true values as N increases.

## Files

### Data Generation
- **`01_gen_subject_dina0.r`**: Generate datasets for multiple sample sizes
  - Creates 10 datasets: N = 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000
  - Uses DINA model with K=2 skills, J=3 items
  - Discrete profile probabilities (use_mvn = FALSE)
  - Saves to `data/subject_dina0_N######.rda`

### Analysis
- **`02_analyze_asymptotic.r`**: Fit MCMC and analyze convergence
  - Fits each dataset with DE-MCMC
  - Extracts parameter estimates and standard errors
  - Creates diagnostic plots
  - Checks asymptotic properties

## Workflow

### Step 1: Generate Data
```bash
cd /media/yslin/Tui/01_Projects/ggdmc/tests/testthat/Group9_gen_cdm
Rscript 01_gen_subject_dina0.r
```

**Output:**
- `data/subject_dina0_N000100.rda` through `data/subject_dina0_N100000.rda`
- `data/generation_metadata.rda` (contains info about all datasets)

**Time:** ~2-5 minutes (depends on N=100000 simulation)

### Step 2: Analyze Asymptotic Behavior
```bash
Rscript 02_analyze_asymptotic.r
```

**Output:**
- `results/asymptotic_results.rda` (all estimates and statistics)
- `results/bias_vs_N.pdf` (bias convergence to 0)
- `results/se_vs_N.pdf` (SE decay as 1/√N)
- `results/coverage_vs_N.pdf` (95% CI coverage)
- `results/estimates_vs_N.pdf` (parameter convergence)

**Time:** ~30-60 minutes (MCMC for 10 datasets)

## Model Specification

### Q-Matrix (3 items × 2 skills)
```
       Algebra  Geometry
Item1    1        0
Item2    0        1
Item3    1        1
```

### True Parameter Values
```
Guess parameters:
  guess1 = 0.10
  guess2 = 0.20
  guess3 = 0.30

Profile probabilities (compact):
  pi1 = 0.1578  [P(no skills)]
  pi2 = 0.2630  [P(Algebra only)]
  pi3 = 0.1508  [P(Geometry only)]
  pi4 = 0.4284  [P(both skills)] = 1 - pi1 - pi2 - pi3

Slip parameters:
  slip1 = 0.20
  slip2 = 0.40
  slip3 = 0.60
```

Total: **9 parameters**

## Expected Asymptotic Behavior

As N → ∞:

1. **Consistency**: Estimates → True values
   - Bias → 0
   - Plot: Estimates should converge to dashed lines

2. **Efficiency**: Standard errors → 0
   - SE decreases as O(1/√N)
   - Plot: Should be approximately linear on log-log scale

3. **Asymptotic Normality**: √N × (estimate - true) → N(0, σ²)
   - Standardized bias should be O(1), not growing with N
   - 95% credible intervals should have ~95% coverage

4. **Coverage**: 95% CI coverage → 0.95
   - May be poor for small N
   - Should stabilize around 0.95 for large N

## Diagnostic Plots

### 1. Bias vs N (`bias_vs_N.pdf`)
Shows bias (estimate - true value) for each parameter across sample sizes.

**Expected pattern:** Lines converge to 0 as N increases

### 2. Standard Error vs N (`se_vs_N.pdf`)
Shows posterior standard errors on log-log scale.

**Expected pattern:** Approximately linear with slope -0.5 (SE ∝ 1/√N)

### 3. Coverage vs N (`coverage_vs_N.pdf`)
Shows proportion of true values within 95% credible intervals.

**Expected pattern:** Approaches 0.95 (dashed line) as N increases

### 4. Estimates vs N (`estimates_vs_N.pdf`)
Shows point estimates with 95% CIs for each parameter.

**Expected pattern:**
- Estimates converge to true values (dashed lines)
- Credible intervals shrink as N increases

## Customization

### Adjust Sample Sizes
In `01_gen_subject_dina0.r`, line 110:
```r
N_values <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
```

### Adjust MCMC Settings
In `02_analyze_asymptotic.r`, lines 52-54:
```r
nmc <- 500      # iterations per chain (increase for better convergence)
nchain <- 3     # number of chains (increase for robustness)
burnin <- 0.5   # burn-in proportion (adjust if chains need longer warmup)
```

### Adjust True Parameters
In `01_gen_subject_dina0.r`, lines 70-74:
```r
p_vector <- c(
    guess1 = .10, guess2 = .20, guess3 = .30,
    pi1 = .1578, pi2 = .2630, pi3 = .1508,
    slip1 = .20, slip2 = .40, slip3 = .60
)
```

## Troubleshooting

### MCMC Not Converging (Small N)
For N=100, chains may not converge well. Consider:
- Increasing `nmc` (more iterations)
- Increasing `burnin` (discard more early samples)
- Checking Rhat (Gelman-Rubin diagnostic)

### Memory Issues (Large N)
For N=100,000:
- Data generation: Uses ~100MB RAM
- MCMC: Uses ~500MB RAM per chain
- Solution: Run analysis on subset of N values first

### Long Runtime
Expected times:
- Data generation: ~5 min for all 10 datasets
- MCMC analysis: ~3-6 min per dataset → ~30-60 min total

To speed up:
- Reduce `nmc` to 200-300 for quick tests
- Reduce `nchain` to 2
- Run only subset of N values initially

## Validation

The asymptotic study is successful if:

✓ Bias → 0 as N increases (all parameters)
✓ SE decreases approximately as 1/√N
✓ Coverage rate → 0.95 for large N
✓ No systematic bias patterns across parameters
✓ Estimates converge smoothly (no jumps or irregularities)

## References

Theory:
- Asymptotic properties of maximum likelihood estimators
- Central Limit Theorem for MCMC
- De la Torre, J. (2009). DINA model and parameter estimation. Psychometrika.

Implementation:
- ggdmc package documentation
- cdModel package for CDM structures
