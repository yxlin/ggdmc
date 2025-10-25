# Systematic Testing Framework - Summary

## What Was Created

A complete systematic testing framework for evaluating DINA model asymptotic behavior across 10 different sample sizes (N = 100 to 100,000).

### New Scripts

| File | Purpose | Runtime |
|------|---------|---------|
| `run_tests.r` | **Interactive launcher** - Menu-driven interface | Instant |
| `01_fit_subject_dina0_systematic.r` | **Systematic testing** - Tests all datasets | Hours |
| `01_fit_subject_dina0_quick.r` | **Quick test** - Tests small datasets only | 5-10 min |
| `02_analyze_asymptotic_results.r` | **Analysis & visualization** - Creates plots | 1-2 min |
| `README_SYSTEMATIC_TESTING.md` | **Documentation** - Complete user guide | - |

### Original Scripts (Preserved)

- `01_fit_subject_dina0.r` - Your original single-dataset script (unchanged)
- `01_subject_dina0.r` - Data generation script (unchanged)

## Quick Start

### Interactive Mode

```bash
cd /media/yslin/Tui/01_Projects/ggdmc/tests/testthat/Group10_cdm_subjects
Rscript run_tests.r
```

Then select:
- **Option 1**: Quick test (recommended first)
- **Option 2**: Full systematic test
- **Option 3**: Analyze existing results
- **Option 4**: Run full workflow

### Command Line Mode

```bash
# Quick test
Rscript run_tests.r 1

# Systematic test
Rscript run_tests.r 2

# Analyze results
Rscript run_tests.r 3

# Full workflow
Rscript run_tests.r 4
```

### Direct Script Execution

```r
# Quick test (5-10 minutes)
source("01_fit_subject_dina0_quick.r")

# Systematic test (several hours)
source("01_fit_subject_dina0_systematic.r")

# Analysis
source("02_analyze_asymptotic_results.r")
```

## Features

### 1. Systematic Testing (`01_fit_subject_dina0_systematic.r`)

**Automatically tests all 10 datasets:**
- Loops through N = 100, 200, 500, 1K, 2K, 5K, 10K, 20K, 50K, 100K
- Two-stage sampling (migration → focused)
- Convergence diagnostics (Gelman-Rubin MPSRF)
- Parameter estimation accuracy
- Timing analysis

**Output:**
```
results/
├── fit_dina0_N000100.rda          # Individual results
├── fit_dina0_N000200.rda
├── ...
├── systematic_test_summary.rda    # Combined results
└── systematic_test_summary.csv    # Summary table (Excel-compatible)
```

**Summary Table Includes:**
- Sample size (N)
- Overall MPSRF (convergence)
- Convergence status (TRUE/FALSE)
- Mean absolute error (MAE)
- Timing (Stage 1, Stage 2, Total)

### 2. Quick Test (`01_fit_subject_dina0_quick.r`)

**For rapid testing/debugging:**
- Tests only N = 100, 500, 1000
- Reduced iterations (50 vs 100)
- Faster thinning
- Results in `results_quick/`

**Use cases:**
- Verifying script changes
- Testing different sampling parameters
- Quick sanity checks

### 3. Analysis & Visualization (`02_analyze_asymptotic_results.r`)

**Creates 6 publication-ready plots:**

1. **convergence_mpsrf.pdf**
   - Gelman-Rubin MPSRF vs sample size
   - Shows which datasets converged

2. **bias_by_parameter.pdf**
   - Bias (estimate - true) for each parameter
   - Shows if bias → 0 as N increases

3. **mae_by_parameter.pdf**
   - Median absolute error by parameter
   - Log-log scale shows convergence rates

4. **mean_mae.pdf**
   - Average MAE across all parameters
   - Overall estimation accuracy trend

5. **mean_mae_fitted.pdf**
   - MAE with power law fit: MAE ~ N^(-α)
   - Estimates asymptotic convergence rate
   - Compares to theoretical rate (1/√N)

6. **computational_time.pdf**
   - Runtime analysis by sample size
   - Breakdown by sampling stage

**Statistical Analysis:**
- Power law fitting: MAE ~ N^(-α)
- Convergence rate estimation
- Correlation analysis
- Comparison to parametric rate (α = 0.5)

### 4. Interactive Launcher (`run_tests.r`)

**User-friendly menu interface:**
```
╔════════════════════════════════════════════════════════════════╗
║        DINA Model Asymptotic Behavior Testing Suite           ║
╚════════════════════════════════════════════════════════════════╝

Available tests:
  1. Quick Test      - Test small datasets only (5-10 min)
  2. Systematic Test - Test all datasets (several hours)
  3. Analyze Results - Analyze and visualize results
  4. Full Workflow   - Run systematic test + analysis
  5. Exit

Enter your choice (1-5):
```

**Safety features:**
- Validates user input
- Confirms before long-running tests
- Checks for required data files
- Works in both interactive and batch mode

## Key Improvements Over Original Script

### Before (01_fit_subject_dina0.r)

❌ Manual editing required to change datasets
❌ One dataset at a time
❌ Results scattered across comments
❌ No systematic comparison
❌ Manual convergence checking

### After (Systematic Framework)

✅ Automatic loop through all datasets
✅ Batch processing
✅ Structured results storage
✅ Automated comparison tables
✅ Automated convergence diagnostics
✅ Publication-ready visualizations
✅ Statistical analysis of convergence rates
✅ CSV export for further analysis

## Expected Results

### Convergence (MPSRF)

Typical pattern:
- N = 100-500: MPSRF may be 1.1-1.5 (marginal convergence)
- N = 1000-5000: MPSRF should be < 1.2
- N = 10000+: MPSRF should be < 1.1 (good convergence)

### Estimation Accuracy

Typical MAE trend:
```
N       Mean MAE    Notes
100     0.2-0.4     High variability
500     0.1-0.2     Improving
1000    0.05-0.15   Good
5000    0.02-0.08   Very good
10000+  < 0.05      Excellent
```

### Asymptotic Rate

Expected power law: **α ≈ 0.4-0.6**
- α = 0.5: Theoretical parametric rate (1/√N) ✓
- α > 0.3: Good convergence
- α < 0.3: Investigate (possible issues)

## Workflow Recommendations

### First Time

1. **Quick test** to verify everything works
   ```r
   source("01_fit_subject_dina0_quick.r")
   ```

2. **Check results** - should complete in ~10 minutes

3. **Review convergence** - adjust parameters if needed

4. **Run systematic test** (when ready for full analysis)
   ```r
   source("01_fit_subject_dina0_systematic.r")
   ```

5. **Analyze results**
   ```r
   source("02_analyze_asymptotic_results.r")
   ```

### Debugging/Development

Use quick test with modified parameters:

```r
# Edit in 01_fit_subject_dina0_quick.r
sample_sizes <- c(100, 500)  # Test only 2 datasets
SAMPLING_CONFIG$nmc <- 25     # Even faster
```

### Production Run

Use systematic test with optimized parameters:

```r
# Edit in 01_fit_subject_dina0_systematic.r
SAMPLING_CONFIG$nmc <- 200    # More iterations for accuracy
SAMPLING_CONFIG$ncore <- 4    # Parallel processing
```

## Customization

### Test Different Sample Sizes

```r
# In 01_fit_subject_dina0_systematic.r, line 18
sample_sizes <- c(1000, 5000, 10000)  # Custom selection
```

### Adjust Sampling Parameters

```r
# In either systematic or quick script
SAMPLING_CONFIG <- list(
  nmc = 150,              # More/fewer iterations
  nchain = 4,             # More chains
  migration_prob = 0.10,  # Higher migration
  thin_stage1 = 10,       # More thinning
  ncore = 4               # Parallel processing
)
```

### Add More Diagnostics

Edit `01_fit_subject_dina0_systematic.r` after line 123 to add:
- Effective sample size
- Autocorrelation checks
- Posterior predictive checks
- Additional parameter comparisons

## Troubleshooting

### "Data file not found"

Check data directory:
```bash
ls ../Group9_gen_cdm/data/subject_dina0_N*.rda
```

If files missing, run data generation script first.

### Poor Convergence

Try:
1. Increase `nmc` (more iterations)
2. Adjust `thin_stage1` and `thin_stage2`
3. Try different `migration_prob`
4. Use `is_pblocked = TRUE`
5. Increase `nchain` for better diagnostics

### Long Runtime

Options:
1. Use quick test instead
2. Reduce `nmc`
3. Test fewer sample sizes
4. Increase `ncore` for parallel processing
5. Run overnight for large datasets

### Memory Issues

For very large datasets (N > 50K):
1. Increase system memory
2. Test one dataset at a time
3. Clean up intermediate objects
4. Use more aggressive thinning

## File Sizes

Approximate result file sizes:

| N       | Result File Size |
|---------|------------------|
| 100     | ~200 KB          |
| 1,000   | ~500 KB          |
| 10,000  | ~2 MB            |
| 50,000  | ~8 MB            |
| 100,000 | ~15 MB           |

Total for all 10 datasets: ~30-50 MB

## Next Steps

After running systematic tests:

1. **Review plots** in `results/plots/`
2. **Check convergence** - rerun if MPSRF > 1.1
3. **Examine bias patterns** - are they decreasing with N?
4. **Verify asymptotic rate** - is α ≈ 0.5?
5. **Export results** - use CSV for external analysis
6. **Publication** - plots are ready for papers/presentations

## Support

- **Documentation**: See `README_SYSTEMATIC_TESTING.md`
- **Examples**: Check script comments
- **Issues**: Consult ggdmc package documentation

---

**Created:** 2025-10-22
**Framework Version:** 1.0
**Compatible with:** ggdmc, ggdmcPrior, ggdmcModel, cdModel
