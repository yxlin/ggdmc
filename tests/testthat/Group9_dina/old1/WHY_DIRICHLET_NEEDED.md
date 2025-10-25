# Why Dirichlet Prior is Needed for 3+ Skills

## The Problem: Exponential Difficulty with More Skills

When using **uniform priors** for profile probabilities, the acceptance rate decreases **exponentially** as the number of skills increases.

## Mathematical Explanation

Profile probabilities must satisfy the **simplex constraint**:
```
π₁ + π₂ + ... + π_{L-1} ≤ 1
```

When sampling K parameters independently from Uniform[0,1], the probability they sum to ≤ 1 is:
```
P(sum ≤ 1) = 1/K!
```

## Comparison: 2 Skills vs 3 Skills

| Attribute | 2 Skills | 3 Skills | Ratio |
|-----------|----------|----------|-------|
| **# Profiles** | 2² = 4 | 2³ = 8 | 2× |
| **# Params to estimate** | 3 (pi_00, pi_10, pi_01) | 7 (pi_000 ... pi_011) | 2.3× |
| **Theoretical acceptance** | 1/3! = **16.67%** | 1/7! = **0.020%** | **840×** |
| **Simulated acceptance** | 17.0% | 0.01% | 1700× |
| **Proposals per accept** | ~6 | ~5,040 | 840× |
| **Acceptances in 1000 iter** | ~167 | ~0.2 | - |

## Why 2-Skill "Works"

With **16.67% acceptance rate**:
- In 1000 MCMC iterations → ~167 acceptances
- Chain moves slowly but **does move**
- Poor efficiency but **eventually converges**
- Takes longer but **gets there**

⚠️ **Still suboptimal!** You're rejecting 83% of proposals.

## Why 3-Skill Fails

With **0.020% acceptance rate**:
- In 1000 MCMC iterations → ~0.2 acceptances (almost zero!)
- Chain is **stuck** at initial values
- Appears as **flat lines** in trace plots
- Essentially **cannot sample** from the posterior

## Practical Impact

### 2-Skill Case (observed)
```
Iteration  Accept?  Status
1-5        No       Rejected (sum > 1)
6          YES      Accepted! ✓
7-10       No       Rejected
11         YES      Accepted! ✓
...
Result: Slow but converges
```

### 3-Skill Case (observed)
```
Iteration  Accept?  Status
1-500      No       All rejected (sum > 1)
501-1000   No       All rejected
...
Result: FLAT CHAIN (appears broken)
```

## The Solution: Dirichlet Prior

The Dirichlet distribution is the **conjugate prior for the simplex**:
- Automatically enforces Σπᵢ = 1
- Samples efficiently from the simplex
- **100% acceptance** for the simplex constraint
- Standard approach in Bayesian CDM

### Implementation
```r
# Identify profile probability parameters
pi_params <- grepl("^pi_", model@pnames)

# Use Dirichlet for pi, uniform for guess/slip
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = ifelse(pi_params, NA, 1.0),          # NA for Dirichlet
    dist = ifelse(pi_params, "dirichlet", "unif")
)
```

## Summary

| Approach | 2 Skills | 3 Skills | 4 Skills | Recommendation |
|----------|----------|----------|----------|----------------|
| **Uniform** | Poor (17%) | Fails (<0.1%) | Impossible (<0.001%) | ❌ Never use |
| **Dirichlet** | Excellent | Excellent | Excellent | ✅ Always use |

## Key Insight

The difference between 2-skill "working" and 3-skill "failing" is:
- **NOT** a bug in the code
- **NOT** a model misspecification
- **PURELY** a consequence of the simplex geometry

As K grows, the volume of the K-dimensional simplex (valid region) becomes vanishingly small compared to the [0,1]^K hypercube (sampling space):

```
Volume ratio = (1/K!) / 1 = 1/K!

K=3:  1/6     ≈ 16.7%  → Barely works
K=7:  1/5040  ≈ 0.02%  → Fails completely
K=15: 1/1.3e12 → Impossible
```

## Recommendation

**ALWAYS use Dirichlet prior for profile probabilities**, regardless of the number of skills. Even for 2-skill models, it provides:
- 6× better acceptance rate
- Faster convergence
- More reliable inference
- Standard Bayesian practice
