# 📦 ggdmc

<!-- Badges -->
[![CRAN Status](https://www.r-pkg.org/badges/version/ggdmc)](https://cran.r-project.org/package=ggdmc)
[![Downloads](https://cranlogs.r-pkg.org/badges/ggdmc)](https://cran.r-project.org/package=ggdmc)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/yxlin/ggdmc/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yxlin/ggdmc/actions/workflows/R-CMD-check.yaml)

`ggdmc` (v0.2.8.9) is an R package for **Bayesian inference** on cognitive **choice response time models**, including:

- **Linear Ballistic Accumulator (LBA)**
- **Diffusion Decision Model (DDM)**

It supports **hierarchical Bayesian modelling**, efficient MCMC sampling, and flexible model definitions.


## 🚀 Getting Started

### Installation

From CRAN:
```r
install.packages("ggdmc")
```

From GitHub (development version):
```r
remotes::install_github("yxlin/ggdmc")
```


## 🔢 Overview

`ggdmc`  helps researchers fit hierarchical Bayesian models of decision-making.
Key features include:

- Flexible **model specification** with condition-dependent parameters
- Efficient **MCMC sampling** with population-level migration
- Support for Truncated normal priors and bounded uniform priors
- Built-in tools for:
  - Posterior summarisation
  - Convergence diagnostics
  - Model comparison 

### Currently supported models:
- LBA (Linear Ballistic Accumulator)
- DDM (Diffusion Decision Model)
- Regular statistical models (via hyper-only route)
- Extendable to user-defined models (e.g., pLBA, DDM with varying drift rates)


## ⚡ Quick Start — Fit a Simple DDM
Below is a step-by-step example of fitting a basic single-subject DDM.
We’ll simulate some data, fit the model, and plot the posterior distributions.


```r
# 1. Load package
library(ggdmc)

# 2. Define a simple DDM model
model <- ggdmcModel::BuildModel(
  p_map = list(a="1", v="1", z="1", sz="1", t0="1"),  # Free parameters
  match_map = list(M = list(s1="r1", s2="r2")),       # Mapping stimuli to responses
  factors = list(S = c("s1", "s2")),                  # Experimental factors
  constants = c(d=0, s=1, st0=0, sv=0, precision=3),  # Fixed parameters
  accumulators = c("r1", "r2"),
  type = "fastdm"
)

# 3. Simulate a small dataset for one subject
true_params <- c(a=1.0, sz=0.25, t0=0.15, v=2.5, z=0.38)
dat <- simulate(model, nsim=300, parameter_vector=true_params, n_subject=1)

# 4. Prepare the data for fitting
dmi <- ggdmcModel::BuildDMI(dat, model)

# 5. Set flat (uninformative) priors
pri <- ggdmcPrior::set_priors(
  ggdmcPrior::BuildPrior(
    p0=rep(0, model@npar), p1=rep(10, model@npar),
    lower=rep(NA, model@npar), upper=rep(NA, model@npar),
    dist=rep("unif", model@npar), log_p=rep(TRUE, model@npar)
  )
)

# 6. Fit the model (MCMC sampling)
fit <- ggdmc::StartSampling_subject(dmi[[1]], pri, thin=2, seed=123)

# 7. Rebuild and plot posterior distributions
post <- ggdmc::RebuildPosterior(fit)
ggdmc::plot(post, pll=FALSE, den=TRUE)
```

💡 **Tip**: For faster tests, lower `nsim` and `nmc` in `setThetaInput()`. For real analyses, increase them for better convergence.


## 🔧 Dependencies

`ggdmc` requires:

- R (≥ 3.3.0)
- C++ integration: `Rcpp`, `RcppArmadillo`
- Data handling: `data.table`, `matrixStats`
- Plotting: `lattice`
- Core `ggdmc` components: `ggdmcHeaders`, `ggdmcModel`, `ggdmcPrior`, `ggdmcLikelihood`
- Model modules: `lbaModel`, `ddModel`


Install all dependencies:

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

- Lin, Y.-S., & Strickland, L. (2020). *Evidence accumulation models with R: A practical guide to hierarchical Bayesian methods*. The Quantitative Methods for Psychology, 16(2), 133–153. [doi.org/10.20982/tqmp.16.2.p133](https://doi.org/10.20982/tqmp.16.2.p133) | [PDF](https://www.tqmp.org/RegularArticles/vol16-2/p133/p133.pdf)
 
- Heathcote, A., Lin, Y.-S., Reynolds, A., Strickland, L., Gretton, M., & Matzke, D. (2018). Dynamic models of choice. *Behavior Research Methods*. [doi.org/10.3758/s13428-018-1067-y](https://doi.org/10.3758/s13428-018-1067-y)

## 👨‍💼 Contributors
The initial version of `ggdmc` was adapted from the **Dynamic Model of Choice** (Heathcote et al., 2018). 

Bug reports and suggestions are welcome via:
- 📧 [Email](mailto:yishinlin001@gmail.com)
- 🐛 [GitHub Issues](https://github.com/yxlin/ggdmc/issues)

## 📓 Acknowledgments

- DDM functions are based on Voss & Voss's `fast-dm` and Gretton's contributions to `rtdists` (Singmann et al.). These were rewritten in C++ to support MCMC-based inference.
- Truncated normal sampling draws from:
  - Jonathan Olmsted's [RcppTN 0.1-8](https://github.com/olmjo/RcppTN)
  - Christopher Jackson's [`msm`](https://cran.r-project.org/web/packages/msm/index.html)
  - Robert, C. P. (1995). *Simulation of truncated normal variables*. *Statistics and Computing, 5*(2), 121–125. https://doi.org/10.1007/BF00143942

---

## Further Examples

### ⚡ Quick Start — Fit a Simple LBA
This example simulates a two-choice LBA for one participant, fits it, and plots the posterior. We keep drift variability fixed (sd_v = 1) and let the core parameters vary: A, B, t0, mean_v.true, mean_v.false.

```# 1. Load package
library(ggdmc)

# 2. Define a simple LBA model (two accumulators)
model <- ggdmcModel::BuildModel(
  p_map  = list(A="1", B="1", t0="1", mean_v="M", st0="1"),  # 'M' => match/mismatch (true/false)
  match_map = list(M = list(s1="r1", s2="r2")),
  factors   = list(S = c("s1","s2")),
  constants = c(st0=0, sd_v=1),                              # fix drift SD
  accumulators = c("r1","r2"),
  type = "lba"
)

# 3. Simulate data for one subject
true_params <- c(A=0.5, B=1.2, t0=0.30, mean_v.true=2.5, mean_v.false=1.5)
dat <- simulate(model, nsim=400, parameter_vector=true_params, n_subject=1)

# 4. Prepare data & priors, then fit
dmi <- ggdmcModel::BuildDMI(dat, model)
pri <- ggdmcPrior::set_priors(
  ggdmcPrior::BuildPrior(
    p0=rep(0, model@npar), p1=rep(10, model@npar),
    lower=rep(NA, model@npar), upper=rep(NA, model@npar),
    dist=rep("unif", model@npar), log_p=rep(TRUE, model@npar)
  )
)
fit  <- ggdmc::StartSampling_subject(dmi[[1]], pri, thin=2, seed=42)

# 5. Inspect posteriors
post <- ggdmc::RebuildPosterior(fit)
ggdmc::plot(post, pll=FALSE, den=TRUE)
```

**Tip**: For a quick smoke test, reduce nsim (e.g., 200) or thin; for real analyses, increase them to improve convergence and stability.


### Example 1: LBA Model with Population Recovery

This example shows how to:

1. Define a condition-dependent LBA model
2. Simulate data
3. Recover parameters at both subject and population levels

See scripts under:

- `tests/testthat/Group1/data/` – simulation data
- `tests/testthat/Group1/0_5param_hyper.r` – hyper-level only model
- `tests/testthat/Group1/1_6param_fit_subject.r` – single participant fitting
- `tests/testthat/Group1/2_v_model_multiple_level.r` – hierarchical model with varying drift rates
- `tests/testthat/Group1/3_B_model_中文.r` – hierarchical model with varying thresholds (Chinese example)

Key functions:

- `simulate()` – generate data
- `StartSampling()` – find optimised parameters
- `compare()` – analyse posterior recovery

---

### Example 2: Minimal DDM Recovery
- `tests/testthat/Group6/data/` – simulation data
- `tests/testthat/Group6//` – fitting the data to find optimised parameters 


