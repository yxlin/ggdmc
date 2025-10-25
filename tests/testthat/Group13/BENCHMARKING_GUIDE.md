# Benchmarking Guide: How to Compare Original vs Fast Versions

## Quick Answer: Which Timing Method to Use?

### For Your Use Case: **rbenchmark** (Recommended)

**Why rbenchmark is best here:**
- ✅ Simple and straightforward API
- ✅ Handles multiple replications automatically
- ✅ Built-in statistical summary (mean, min, max, relative)
- ✅ Good for comparing 2-3 variants
- ✅ Automatically handles warm-up runs
- ✅ Perfect for MCMC benchmarks (seconds to minutes runtime)

**When to use alternatives:**
- `system.time()`: Quick one-off timing (no statistics)
- `microbenchmark`: For very fast operations (microseconds to milliseconds)
- `bench`: Modern alternative with memory profiling (but overkill here)

## Comparison of Timing Methods

### 1. system.time() - Basic but Limited

```r
system.time({
    result <- run_subject(config, dmi, samples)
})
```

**Pros:**
- Built into base R (no dependencies)
- Simple to use
- Shows user, system, and elapsed time

**Cons:**
- ✗ No automatic replications
- ✗ No statistical summary
- ✗ Manual warm-up needed
- ✗ More code to compare multiple versions

**Best for:** Quick one-off timing checks

---

### 2. rbenchmark::benchmark() - Recommended for This Case

```r
library(rbenchmark)

benchmark(
    original = run_original_twice(),
    fast = run_fast_twice(),
    replications = 3,
    columns = c("test", "elapsed", "relative")
)
```

**Output example:**
```
      test replications elapsed relative
1 original            3  15.234    1.156
2     fast            3  13.181    1.000
```

**Pros:**
- ✅ Simple API, easy to read
- ✅ Automatic replications
- ✅ Automatic relative comparison
- ✅ Good for longer-running operations (your case)
- ✅ Shows user/system/elapsed time

**Cons:**
- ✗ Requires installation
- ✗ Less detailed than microbenchmark
- ✗ No memory profiling

**Best for:** Comparing 2-5 variants of functions that run for seconds

---

### 3. microbenchmark::microbenchmark() - Overkill Here

```r
library(microbenchmark)

microbenchmark(
    original = run_original_twice(),
    fast = run_fast_twice(),
    times = 10
)
```

**Pros:**
- ✅ Very precise (nanosecond resolution)
- ✅ Many replications (default 100)
- ✅ Nice summary statistics
- ✅ Good visualization with ggplot2

**Cons:**
- ✗ Overkill for long operations (your MCMC runs take seconds)
- ✗ More complex setup
- ✗ 100 reps would take too long for MCMC

**Best for:** Comparing very fast operations (microseconds to milliseconds)

---

### 4. bench::mark() - Modern but Heavy

```r
library(bench)

mark(
    original = run_original_twice(),
    fast = run_fast_twice(),
    iterations = 3,
    check = FALSE
)
```

**Pros:**
- ✅ Modern tidyverse-friendly
- ✅ Memory profiling included
- ✅ Detailed GC statistics
- ✅ Nice output format

**Cons:**
- ✗ More dependencies
- ✗ Overkill for simple comparisons
- ✗ More complex interpretation

**Best for:** When you need memory profiling and detailed statistics

---

## Recommended Approach for Your MCMC Benchmark

### Step 1: Install rbenchmark (if needed)

```r
install.packages("rbenchmark")
```

### Step 2: Use the Updated Script

```bash
cd /media/yslin/Tui/01_Projects/ggdmc/tests/testthat/Group13
Rscript benchmark_simple.R
```

The script automatically:
1. Detects if rbenchmark is available
2. Falls back to system.time() if not
3. Runs both versions 3 times
4. Calculates speedup and statistics
5. Shows C++ timing breakdown

### Step 3: Interpret Results

The script will show:

```
============================================================
Benchmark Results (rbenchmark)
============================================================
      test replications elapsed relative user.self sys.self
1 original            3  45.678    1.156    45.123    0.234
2     fast            3  39.512    1.000    38.987    0.198

------------------------------------------------------------
Summary:
  Original version:  45.678 seconds
  Fast version:      39.512 seconds
  Speedup:           1.156x
  Time saved:        6.166 seconds (13.5%)
------------------------------------------------------------

============================================================
C++ Internal Timing Breakdown (Fast Version)
============================================================
Total C++ time:      38.234 seconds (100.0%)
  Crossover:         12.456 seconds (32.6%)
  Migration:         3.123 seconds (8.2%)
  Likelihood:        22.655 seconds (59.2%)
  Other (MH, etc):   2.778 seconds (7.3%)
```

**Interpretation:**
- **Speedup < 1.05x**: Chain shuffling wasn't the bottleneck
- **Speedup 1.05-1.15x**: Modest improvement, focus on likelihood next
- **Speedup > 1.15x**: Significant improvement, continue with Step 2 optimizations

## Understanding the Timing Breakdown

### Two Levels of Timing

1. **R-level timing** (from rbenchmark or system.time)
   - Includes R overhead, object creation, function dispatch
   - Total elapsed time from R's perspective

2. **C++-level timing** (from TimingStats)
   - Pure C++ execution time
   - Excludes R/C++ interface overhead
   - More accurate for identifying bottlenecks

### Key Metrics

**Speedup = Original_Time / Fast_Time**
- 1.0x = No improvement
- 1.1x = 10% faster
- 2.0x = Twice as fast (50% time saved)

**Percentage breakdown:**
- Crossover + Migration = DE algorithm overhead
- Likelihood = Model-specific computation
- Other = MH acceptance, random number generation, etc.

If **Likelihood > 70%**: Focus optimization efforts on likelihood functions next
If **Crossover > 30%**: Chain shuffling optimization was important

## Best Practices

### 1. Warm-up Runs
Both rbenchmark and the manual system.time() approach in the script handle this automatically by running multiple times.

### 2. Number of Replications
- **3-5 replications**: Good for MCMC (already slow)
- **10+ replications**: Better statistics but takes longer
- **100+ replications**: Only for very fast functions

### 3. Fair Comparison
The script ensures fairness by:
- Using identical data (config, dmi, samples)
- Running both versions twice (as required for correct results)
- Running in the same R session (same state)
- Using wrapper functions (same interface)

### 4. Statistical Stability
If results vary widely between runs:
- Increase replications (5-10)
- Check for background processes
- Use `nice` to set process priority
- Consider using bench::press() for parameter studies

## Example: Running Your Benchmark

```bash
# Navigate to test directory
cd /media/yslin/Tui/01_Projects/ggdmc/tests/testthat/Group13

# Run the benchmark
Rscript benchmark_simple.R

# Or from within R
R
> source("benchmark_simple.R")
```

## Troubleshooting

**Q: "rbenchmark not found"**
A: The script automatically falls back to system.time(). Install rbenchmark for better results:
```r
install.packages("rbenchmark")
```

**Q: "Results vary too much between runs"**
A: Increase replications in the script (line 59):
```r
replications = 10,  # Increase from 3 to 10
```

**Q: "Want more detailed profiling"**
A: Use Rprof() for line-by-line profiling:
```r
Rprof("profile.out")
result <- run_fast_twice()
Rprof(NULL)
summaryRprof("profile.out")
```

**Q: "Memory usage concerns"**
A: Use bench::mark() which includes memory profiling:
```r
bench::mark(
    original = run_original_twice(),
    fast = run_fast_twice(),
    iterations = 3,
    memory = TRUE
)
```

## Summary Recommendation

**For your MCMC benchmarks:**

✅ **Use rbenchmark** (recommended)
- Simple, reliable, perfect for this use case
- The updated benchmark_simple.R uses it

✅ **Fallback to system.time()** (good enough)
- Built-in, no dependencies
- Script handles this automatically

❌ **Don't use microbenchmark** (overkill)
- Designed for fast operations (μs-ms)
- Your operations take seconds

❌ **Don't use bench unless needed** (too complex)
- Only if you need memory profiling
- More dependencies and complexity

The updated `benchmark_simple.R` script implements the best approach for your use case!
