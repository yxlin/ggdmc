# Quick Reference Card - DINA Systematic Testing

## One-Liners

```bash
# Interactive menu (recommended)
Rscript run_tests.r

# Quick test (5-10 min)
Rscript run_tests.r 1

# Full systematic test (hours)
Rscript run_tests.r 2

# Analyze results
Rscript run_tests.r 3

# Complete workflow
Rscript run_tests.r 4
```

## Or in R

```r
# Quick test
source("01_fit_subject_dina0_quick.r")

# Systematic test
source("01_fit_subject_dina0_systematic.r")

# Analysis
source("02_analyze_asymptotic_results.r")
```

## What Gets Tested

| Dataset | Sample Size | Expected Runtime |
|---------|-------------|------------------|
| Quick   | 100, 500, 1K | 5-10 min |
| Systematic | 100-100K (10 files) | Several hours |

## Key Outputs

### Systematic Test
```
results/
├── fit_dina0_N000100.rda          ← Individual results
├── fit_dina0_N000200.rda
├── ...
├── systematic_test_summary.rda    ← All results
└── systematic_test_summary.csv    ← Summary table
```

### Analysis
```
results/plots/
├── convergence_mpsrf.pdf          ← Gelman-Rubin diagnostics
├── bias_by_parameter.pdf          ← Parameter bias
├── mae_by_parameter.pdf           ← Estimation accuracy
├── mean_mae.pdf                   ← Overall MAE trend
├── mean_mae_fitted.pdf            ← Power law fit
└── computational_time.pdf         ← Runtime analysis
```

## Interpreting Results

### Convergence (MPSRF)
- **< 1.1**: ✓ Converged
- **1.1-1.2**: ⚠ Marginal (acceptable for exploration)
- **> 1.2**: ✗ Not converged (rerun with more iterations)

### Estimation Accuracy (MAE)
- Should **decrease** with N
- Typical: MAE ~ 1/√N
- Check if MAE < 0.05 for large N

### Asymptotic Rate (α)
From power law: MAE ~ N^(-α)
- **α ≈ 0.5**: ✓ Theoretical parametric rate
- **α > 0.3**: Good convergence
- **α < 0.3**: Investigate

## Customization

### Test Fewer Datasets

Edit line 18 in `01_fit_subject_dina0_systematic.r`:
```r
sample_sizes <- c(100, 1000, 10000)  # Test only 3 sizes
```

### More Iterations

Edit SAMPLING_CONFIG:
```r
SAMPLING_CONFIG <- list(
  nmc = 200,    # More iterations
  nchain = 4,   # More chains
  ncore = 4     # Parallel
)
```

### Faster Testing

Edit in quick script:
```r
sample_sizes <- c(100, 500)  # Only 2 datasets
SAMPLING_CONFIG$nmc <- 25    # Fewer iterations
```

## Common Issues

| Problem | Solution |
|---------|----------|
| Data not found | Check `../Group9_gen_cdm/data/` exists |
| Poor convergence | Increase `nmc`, adjust thinning |
| Too slow | Use quick test, reduce `nmc`, increase `ncore` |
| Memory issues | Test one dataset at a time, increase thinning |
| No plots | Run systematic test before analysis |

## Files

| File | Purpose |
|------|---------|
| `run_tests.r` | 🚀 **Start here** - Interactive menu |
| `01_fit_subject_dina0_systematic.r` | Full testing |
| `01_fit_subject_dina0_quick.r` | Quick testing |
| `02_analyze_asymptotic_results.r` | Analysis |
| `README_SYSTEMATIC_TESTING.md` | 📖 Full documentation |
| `SYSTEMATIC_TESTING_SUMMARY.md` | Complete guide |
| `QUICK_REFERENCE.md` | This file |

## Workflow

```
┌─────────────────┐
│  Quick Test     │  (5-10 min, verify setup)
│  run_tests.r 1  │
└────────┬────────┘
         │
         ↓ OK?
         │
┌────────┴─────────┐
│ Systematic Test  │  (hours, test all datasets)
│  run_tests.r 2   │
└────────┬─────────┘
         │
         ↓
         │
┌────────┴─────────┐
│  Analyze         │  (1-2 min, create plots)
│  run_tests.r 3   │
└──────────────────┘
         │
         ↓
   Review plots in results/plots/
```

## Typical Runtime

| N | Iterations | Chains | Runtime |
|---|------------|--------|---------|
| 100 | 100 | 3 | ~5 min |
| 1,000 | 100 | 3 | ~10 min |
| 10,000 | 100 | 3 | ~30 min |
| 50,000 | 100 | 3 | ~60 min |
| All 10 | 100 | 3 | ~3-5 hours |

*Times vary by hardware and configuration*

## Getting Help

1. See `README_SYSTEMATIC_TESTING.md` for details
2. See `SYSTEMATIC_TESTING_SUMMARY.md` for features
3. Check script comments for inline documentation
4. Review original script: `01_fit_subject_dina0.r`

---

**Quick Reference v1.0** | Created: 2025-10-22
