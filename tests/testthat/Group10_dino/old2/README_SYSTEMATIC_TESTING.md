# Systematic DINA Model Asymptotic Behavior Testing

This directory contains scripts for systematically testing the asymptotic behavior of the DINA model across different sample sizes.

## Overview

The testing workflow consists of three main steps:
1. **Data Generation** - Generate synthetic datasets with varying sample sizes
2. **Model Fitting** - Fit DINA models to each dataset systematically
3. **Analysis** - Analyze convergence and asymptotic behavior

## Directory Structure

```
Group10_cdm_subjects/
├── 01_fit_subject_dina0.r                    # Original single-dataset fitting script
├── 01_fit_subject_dina0_systematic.r         # NEW: Systematic testing (all datasets)
├── 01_fit_subject_dina0_quick.r              # NEW: Quick test (small datasets only)
├── 02_analyze_asymptotic_results.r           # NEW: Analysis and visualization
├── data/                                      # Fitting results storage
├── results/                                   # Systematic test results
│   ├── fit_dina0_N000100.rda                 # Individual results by N
│   ├── fit_dina0_N000200.rda
│   ├── ...
│   ├── systematic_test_summary.rda           # Overall summary
│   ├── systematic_test_summary.csv           # Summary table (CSV)
│   └── plots/                                # Visualization outputs
│       ├── convergence_mpsrf.pdf
│       ├── bias_by_parameter.pdf
│       ├── mae_by_parameter.pdf
│       ├── mean_mae.pdf
│       ├── mean_mae_fitted.pdf
│       └── computational_time.pdf
└── README_SYSTEMATIC_TESTING.md              # This file
```

## Available Datasets

Datasets are located in `../Group9_gen_cdm/data/` with the following sample sizes:

| Dataset File                  | Sample Size | File Size |
|-------------------------------|-------------|-----------|
| subject_dina0_N000100.rda     | 100         | 2.2K      |
| subject_dina0_N000200.rda     | 200         | 3.0K      |
| subject_dina0_N000500.rda     | 500         | 5.3K      |
| subject_dina0_N001000.rda     | 1,000       | 11K       |
| subject_dina0_N002000.rda     | 2,000       | 21K       |
| subject_dina0_N005000.rda     | 5,000       | 53K       |
| subject_dina0_N010000.rda     | 10,000      | 110K      |
| subject_dina0_N020000.rda     | 20,000      | 203K      |
| subject_dina0_N050000.rda     | 50,000      | 481K      |
| subject_dina0_N100000.rda     | 100,000     | 929K      |

## Scripts

### 1. Quick Test (Recommended for Development)

**File:** `01_fit_subject_dina0_quick.r`

**Purpose:** Test a few small datasets quickly to verify the workflow.

**Usage:**
```r
source("01_fit_subject_dina0_quick.r")
```

**Configuration:**
- Tests N = 100, 500, 1000 by default
- Reduced iterations (50 instead of 100)
- Faster thinning
- Results saved to `results_quick/`

**When to use:**
- Testing script changes
- Debugging
- Quick verification of workflow

**Expected runtime:** 5-10 minutes

### 2. Systematic Test (Full Analysis)

**File:** `01_fit_subject_dina0_systematic.r`

**Purpose:** Test all available datasets systematically.

**Usage:**
```r
source("01_fit_subject_dina0_systematic.r")
```

**What it does:**
- Loops through all 10 datasets (N = 100 to 100,000)
- Runs two-stage sampling for each:
  - Stage 1: Migration-based exploration
  - Stage 2: Focused sampling
- Computes diagnostics (Gelman-Rubin MPSRF)
- Compares estimates to true values
- Saves individual and summary results

**Configuration (editable in script):**
```r
SAMPLING_CONFIG <- list(
  nmc = 100,                    # Number of iterations
  nchain = 3,                   # Number of chains
  migration_prob = 0.05,        # Stage 1 migration probability
  thin_stage1 = 8,              # Stage 1 thinning
  thin_stage2 = 4,              # Stage 2 thinning
  is_pblocked = FALSE,
  seed = 9032,
  ncore = 1                     # Increase for parallel processing
)
```

**Output:**
- Individual `.rda` files for each N in `results/`
- `systematic_test_summary.rda` - Complete results
- `systematic_test_summary.csv` - Summary table

**Expected runtime:** Several hours (depends on configuration)

### 3. Analysis and Visualization

**File:** `02_analyze_asymptotic_results.r`

**Purpose:** Analyze results and create visualizations.

**Usage:**
```r
source("02_analyze_asymptotic_results.r")
```

**Prerequisites:** Must run systematic test first.

**What it does:**
1. **Convergence Analysis**
   - Checks Gelman-Rubin MPSRF for each dataset
   - Plots MPSRF vs sample size

2. **Estimation Accuracy**
   - Computes bias (estimate - true value) for each parameter
   - Computes median absolute error (MAE)
   - Plots bias and MAE vs sample size

3. **Asymptotic Rate Analysis**
   - Fits power law: MAE ~ N^(-α)
   - Estimates convergence rate α
   - Compares to theoretical rate (1/√N, α = 0.5)

4. **Computational Time**
   - Analyzes runtime by sample size
   - Breaks down by sampling stage

**Output:**
All plots saved to `results/plots/`:
- `convergence_mpsrf.pdf` - Convergence diagnostics
- `bias_by_parameter.pdf` - Parameter-wise bias
- `mae_by_parameter.pdf` - Parameter-wise MAE
- `mean_mae.pdf` - Average MAE trend
- `mean_mae_fitted.pdf` - MAE with power law fit
- `computational_time.pdf` - Runtime analysis

## Workflow

### Quick Start

For quick verification:
```r
# 1. Quick test (small datasets only)
source("01_fit_subject_dina0_quick.r")
```

### Full Analysis

For complete asymptotic behavior analysis:

```r
# 1. Run systematic testing (this will take time!)
source("01_fit_subject_dina0_systematic.r")

# 2. Analyze results and create plots
source("02_analyze_asymptotic_results.r")
```

### Custom Testing

To test specific sample sizes, modify the `sample_sizes` vector:

```r
# In 01_fit_subject_dina0_systematic.r (line 18)
sample_sizes <- c(100, 1000, 5000, 10000)  # Test only these
```

## Interpreting Results

### Convergence (MPSRF)

- **MPSRF < 1.1**: Converged ✓
- **MPSRF > 1.1**: Not converged, consider:
  - More iterations (`nmc`)
  - Different thinning
  - Different migration probability

### Estimation Accuracy

**Bias:**
- Should approach 0 as N increases
- Persistent bias may indicate model misspecification

**MAE (Median Absolute Error):**
- Should decrease with N
- Typically: MAE ~ 1/√N for parametric models
- Power law fit provides convergence rate α

**Asymptotic Rate (α):**
- **α ≈ 0.5**: Consistent with parametric rate ✓
- **α > 0.3**: Good convergence
- **α < 0.3**: Slow convergence, investigate

## Customization

### Change Sampling Parameters

Edit `SAMPLING_CONFIG` in the systematic script:

```r
SAMPLING_CONFIG <- list(
  nmc = 200,           # More iterations
  nchain = 4,          # More chains
  ncore = 4            # Parallel processing
)
```

### Test Different Models

The scripts can be adapted for other CDM rules:
1. Generate data for different rules (DINO, RRUM, etc.)
2. Update `data_path` pattern in systematic script
3. Adjust parameter names in analysis script

### Save Intermediate Results

For long runs, uncomment save statements in the systematic script to save progress after each dataset.

## Troubleshooting

### "Data file not found"
- Check that datasets exist in `../Group9_gen_cdm/data/`
- Verify file naming pattern matches: `subject_dina0_N%06d.rda`

### "Package not found"
Install required packages:
```r
install.packages(c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel"))
```

### Long runtime
- Use quick test first: `01_fit_subject_dina0_quick.r`
- Reduce `nmc` in `SAMPLING_CONFIG`
- Test fewer sample sizes
- Increase `ncore` for parallel processing

### Poor convergence
- Increase `nmc` (more iterations)
- Adjust `thin_stage1` and `thin_stage2`
- Try different `migration_prob`
- Use `is_pblocked = TRUE`

## Expected Results

### Convergence
- Small N (< 1000): May have higher MPSRF
- Large N (> 5000): Should converge well (MPSRF < 1.1)

### Estimation
- Bias should decrease with N
- MAE typically decreases as 1/√N
- Parameters with larger true values may have larger absolute errors

### Timing
- Runtime increases with N
- Stage 1 typically takes longer than Stage 2

## References

- Original single-dataset script: `01_fit_subject_dina0.r`
- Data generation: `../Group9_gen_cdm/`
- Helper functions: `../Group9_gen_cdm/00_helpers.r`

## Contact

For questions or issues, check the main ggdmc documentation or consult the development team.

---

**Last Updated:** 2025-10-22
**Author:** Systematic testing framework for DINA asymptotic behavior
