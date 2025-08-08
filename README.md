# 📦 ggdmc

<!-- Badges -->
[![CRAN Status](https://www.r-pkg.org/badges/version/ggdmc)](https://cran.r-project.org/package=ggdmc)
[![Downloads](https://cranlogs.r-pkg.org/badges/ggdmc)](https://cran.r-project.org/package=ggdmc)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/yxlin/ggdmc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yxlin/ggdmc/actions/workflows/R-CMD-check.yaml)

`ggdmc` (version 0.2.8.9) provides tools for conducting **Bayesian inference** on a range of **choice response time models**, such as the **Linear Ballistic Accumulator (LBA)** and **Diffusion Decision Model (DDM)**.

## 🚀 Getting Started

### ✨ Installation

Install the package from CRAN:

```r
install.packages("ggdmc")
```

---

## 🔢 Overview

`ggdmc` supports *hierarchical Bayesian inference* with cognitive models of decision making. It includes:

- Flexible model specification with support for condition-dependent parameters
- Efficient MCMC sampling with population-level migration
- Truncated normal priors and bounded uniform priors
- Built-in tools for posterior summarisation and convergence diagnostics

Models supported:
- LBA (Linear Ballistic Accumulator)
- DDM (Diffusion Decision Model)
- Extendable to user-defined models


## 🔧 Dependencies

Requires:

- R (≥ 3.3.0)
- C++ integration: `Rcpp`, `RcppArmadillo`
- Data processing: `data.table`, `matrixStats`, `lattice`
- Core components: `ggdmcHeaders`, `ggdmcModel`, `ggdmcPrior`, `ggdmcLikelihood`
- Supported models: `lbaModel`, `ddModel`


Installation:

```r
install.packages(c(
  "Rcpp", "RcppArmadillo", "data.table", 
  "matrixStats", "lattice",
  "ggdmcHeaders", "ggdmcModel", 
  "ggdmcPrior", "ggdmcLikelihood",
  "lbaModel", "ddModel"
))
```



## 📄 Citation

If you use `ggdmc`, please cite:

- Lin, Y.-S., & Strickland, L. (2020). *Evidence accumulation models with R: A practical guide to hierarchical Bayesian methods*. The Quantitative Methods for Psychology, 16(2), 133–153. https://doi.org/10.20982/tqmp.16.2.p133. [PDF](https://www.tqmp.org/RegularArticles/vol16-2/p133/p133.pdf)
 
- Heathcote, A., Lin, Y.-S., Reynolds, A., Strickland, L., Gretton, M., & Matzke, D. (2018). Dynamic models of choice. *Behavior Research Methods*. https://doi.org/10.3758/s13428-018-1067-y


## 👨‍💼 Contributors

The early version of `ggdmc` was adapted from the Dynamic Model of Choice (Heathcote et al., 2018). Bug reports and suggestions are welcome via [(email)](mailto:yishinlin001@gmail.com) or [GitHub Issues](https://github.com/yxlin/ggdmc/issues).

### Example 1: LBA Model with Population Recovery

This example demonstrates how to define a condition-dependent LBA model, simulate data, and recover parameters at both subject and population levels.

▶️ See script under: `tests/testthat/Group1/data/`

```r
model <- BuildModel(
  p_map = list(A = "1", B = c("S", "COLOR"), t0 = "1", mean_v = c("NOISE", "M"), sd_v = "M", st0 = "1"),
  match_map = list(M = list(left = "z_key", right = "x_key")),
  factors = list(S = c("left", "right"), COLOR = c("red", "blue"), NOISE = c("high", "moderate", "low")),
  constants = c(st0 = 0, sd_v.false = 1),
  accumulators = c("z_key", "x_key"),
  type = "lba"
)
```

Run `simulate()`, `StartSampling()`, and `compare()` to analyse the recovery performance.

---

### Example 2: Minimal DDM Recovery

Defines and fits a simple DDM in a hierarchical setting.

▶️ See script under: `tests/testthat/Group6/data/` [to be updated.]

```r
model <- BuildModel(
  p_map = list(a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1", t0 = "1", st0 = "1", s = "1", precision = "1"),
  match_map = list(M = list(s1 = "r1", s2 = "r2")),
  factors = list(S = c("s1", "s2")),
  constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
  accumulators = c("r1", "r2"),
  type = "fastdm"
)
```

Use `simulate()`, `setThetaInput()`, and `StartSampling()` to conduct full recovery.

---

## 📓 Acknowledgments

- DDM functions are based on Voss & Voss's `fast-dm` and Gretton's contributions to `rtdists` (Singmann et al.). These were rewritten in C++ to support MCMC-based inference.
- Truncated normal sampling draws from:
  - Jonathan Olmsted's [RcppTN 0.1-8](https://github.com/olmjo/RcppTN)
  - Christopher Jackson's [`msm`](https://cran.r-project.org/web/packages/msm/index.html)
  - Robert, C. P. (1995). *Simulation of truncated normal variables*. *Statistics and Computing, 5*(2), 121–125. https://doi.org/10.1007/BF00143942

---


