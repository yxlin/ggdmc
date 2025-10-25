# Group14: HB Optimization Benchmark

This directory contains benchmarking scripts for the Step 2 optimization (cached hyper-likelihood for HB models).

## Files

- **`setup_benchmark_data.R`** - Generates benchmark data from Group0
- **`hb_benchmark.R`** - Main benchmark script (clean, production-ready)
- **`hb_benchmark.rda`** - Generated benchmark data (config, dmis, samples)
- **`benchmark_results.txt`** - Saved results summary (generated after benchmark)

## Quick Start

### 1. Generate Benchmark Data

```bash
Rscript setup_benchmark_data.R
```

This will:
- Load test data from Group0_gen_data
- Configure HB model parameters
- Save config, dmis, and samples to `hb_benchmark.rda`

### 2. Run Benchmark

```bash
Rscript hb_benchmark.R
```

This will:
- Compare 3 versions (original, fast with timing, fast without timing)
- Run 3 replications of each
- Show detailed timing breakdown
- Provide recommendations
- Save results to `benchmark_results.txt`

## From R Console

```r
# Generate data
source("setup_benchmark_data.R")

# Run benchmark
source("hb_benchmark.R")
```

## Expected Runtime

Depends on:
- Number of subjects
- Number of iterations
- Model complexity

Typical: 5-15 minutes for full benchmark (3 versions × 3 replications)

## Interpreting Results

The benchmark will show:

1. **Speedup factor** - How much faster the optimized versions are
2. **Time breakdown** - Where time is spent (C++ vs R, likelihood vs sampling)
3. **Recommendation** - Which version to use for your model

### Expected Speedup

| Subjects | Expected Speedup |
|----------|------------------|
| 5-10     | 1.2-1.3x (20-30% faster) |
| 10-20    | 1.3-1.5x (30-40% faster) |
| 20+      | 1.5-1.8x (40-50% faster) |

## Output Files

After running benchmark:
- **Console output** - Detailed results and interpretation
- **benchmark_results.txt** - Summary for reference

## Troubleshooting

**"Data file not found"**
→ Run `setup_benchmark_data.R` first

**"rbenchmark not found"**
→ Install with: `install.packages('rbenchmark')`
→ Or script will fall back to `system.time()`

**Results seem wrong**
→ Check that Group0 data exists
→ Verify setup completed successfully

## What's Being Tested

The optimization eliminates redundant hyper-likelihood calculations:

**Original**: Recalculates sumloghlike(ALL subjects) hundreds of times
**Optimized**: Caches values, only recalculates when subjects change

**Benefit scales with number of subjects!**

## Documentation

See main package documentation:
- `../../HB_OPTIMIZATION_STEP2.md` - Technical details
- `../../QUICK_REFERENCE.md` - User guide
- `../../OPTIMIZATION_SUMMARY.md` - Complete overview
