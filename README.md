# Bayesian Cognitive Modelling

_ggdmc_, evolving from dynamic model of choice (_DMC_, Heathcote et al., 2018),
is a generic tool for hierarchical Bayesian Computations. 

1. Instead of using Gibbs or HMC, _ggdmc_ uses population-based MCMC (pMCMC) 
samplers. A notable Gibbs example is the Python-based 
HDDM (Wiecki, Sofer & Frank, 2013), which has yet provided convenient 
interface to examine the variabilities of parameters. 

2. _ggdmc_ implements two different variants of _migration_ operator. They
are by default used during the burning period to search the target parameter
space. 

3. _ggdmc_ uses two parallel methods. First is via the _parallel_ package in R.
This facilitates the computations of fitting many participants separately.  
The second is via OpenMP library. This facilitates the computations for 
fitting hierarchical models.  

For an advanced parallel computation technique / algorithm, please see my CUDA 
C, R package, [_ppda_](https://github.com/yxlin/ppda), which implements GPU 
parallel computations.

## Getting Started
The getting-started example uses Wiener diffusion model. The modelling of the
Wiener diffusion model, hierarchical or not, can be finished by a couple of 
seconds.  For other more time-demanding examples, see my 
[tutorials site](https://yxlin.github.io/) for more details. 

The names of _R_ functions attemp to inform the user what the functions are for,
, such _BuildModel_.  As the user is usually warned in Bayesian tools, please
use with your own risk. That is, the user must always conduct model checking
after model fits.

Note the sequence of parameters in a parameter vector (i.e., p.vector) must 
follow the sequence in the _p.vector_ reported by _BuildModel_. 
Otherwise, _run_ function will throw errors to stop model fitting. 

## Fit a fixed-effect model to a participant

```
## Set up model ----
## fixing sv & sz to 0, makes "rd" (Ratcliff's diffusion model) to 
## Wiener diffusion model
require(ggdmc)
model <- BuildModel(
  p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1", 
                   st0="1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1","r2"),
  constants = c(st0 = 0, d = 0, sv = 0, sz = 0),  
  type      = "rd")   

npar <- length(GetPNames(model))
p.vector <- c(a=1, v=1.5, z=0.5, t0=.15)
dat <- simulate(model, nsim = 50, ps = p.vector)
dmi <- BuildDMI(dat, model)

p.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1=c(a=1, v=0, z=1, t0=1),
  p2=c(a=1, v=2, z=1, t0=1),
  lower = c(0, -5, rep(0, 2)),
  upper = rep(NA, npar))

## Fit model ----
fit0 <- StartNewsamples(dmi, p.prior)
fit  <- run(fit0, 5e2)

## Model checking -----------
plot(fit)
plot(fit, den = TRUE)
plot(fit, pll=FALSE)
plot(fit, pll=FALSE, den = TRUE)

## 
gelman(fit)
est <- summary(fit, recovery = TRUE, ps = p.vector, verbose = TRUE)

```

## How to fit fixed-effect and hierarchical model with multiple participants

```
library(ggdmc);

model <- BuildModel(
  p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1", 
                   st0="1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1","r2"),
  constants = c(st0 = 0, d = 0, sv = 0, sz = 0),
  type      = "rd")

npar <- length(GetPNames(model))
pop.mean  <- c(a=2,   v=4, z=0.5, t0=0.3)
pop.scale <- c(a=0.5, v=.5, z=0.1, t0=0.05)
pop.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale,
    lower = c(0,-5,  0, 0),
    upper = c(5, 7,  1, 1))

## Simulate some data
dat <- simulate(model, nsub = 50, nsim = 30, prior = pop.prior)
dmi <- BuildDMI(dat, model)
ps <- attr(dat, "parameters")

p.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale*5,
    lower = c(0,-5, 0, 0),
    upper = c(5, 7, 1, 1))

plot(p.prior, ps = ps)  ## Check if all true pvectors in the range of prior

## Sampling separately
fit0 <- StartNewsamples(dmi, p.prior, ncore=2)
fit  <- run(fit0, 5e2, ncore=2)
fit  <- run(fit, 1e2, add=TRUE, ncore=2)  ## add additional 100 samples

## Check model check -----
gelman(fit, verbose=TRUE)
plot(fit)
est0 <- summary(fit, recovery = TRUE, ps = ps, verbose =TRUE)

## Sampling hierarchically
  mu.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale*5,
    lower = c(0,-5,  0, 0),
    upper = c(5, 7,  1, 1)
    )

  sigma.prior <- BuildPrior(
    dists = rep("beta", npar),
    p1    = c(a=1, v=1, z=1, t0=1),
    p2    = rep(1, npar),
    upper = rep(1, npar))

  ## Note the names are important
  priors <- list(pprior=p.prior, location=mu.prior, scale=sigma.prior)
  names(priors)
  # [1] "pprior"   "location" "scale"

## Fit hierarchical model ----
fit0 <- StartNewsamples(dmi, priors)
fit  <- run(fit0, 5e2)

p0 <- plot(fit, hyper = TRUE)
p0 <- plot(fit, hyper = TRUE, den = TRUE, pll=FALSE)

## Check model -----------
res <- hgelman(fit, verbose = TRUE)
est0 <- summary(fit, recovery = TRUE, ps = ps, verbose = TRUE)
est1 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.mean,  type = 1, verbose = TRUE)
est2 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.scale, type = 2, verbose = TRUE)


for(i in 1:length(fit))
{
  est <- summary(fit[[i]], recovery = TRUE, ps = ps[i,], verbose=TRUE)
}

```

## List of models currently hard-wired in _ggdmc_
1. The LBA model, type = "norm",
2. The DDM, type = "rd",

## PDA-based models 
3. The Piecewise LBA model 0; CPU-based PDA likelihoods; type = "plba0",
4. The Piecewise LBA model 1; CPU-based PDA likelihoods; type = "plba1", 
5. The Piecewise LBA model 0; GPU-based PDA likelihoods; type = "plba0_gpu", 
6. The Piecewise LBA model 1; GPU-based PDA likelihoods; type = "plba1_gpu", 
7. The LBA model; GPU-based PDA likelihoods;, type = "norm_pda_gpu",
8. The correlated accumualtor model; type = "cnorm".

For the details regarding PLBA types, please see 
[Holmes, Trueblood, and Heathcote (2016)](http://dx.doi.org/10.1016/j.cogpsych.2015.11.002)

## Further information  
Please see my [tutorials site, Cognitive Model](https://yxlin.github.io/), for 
more details.

## Prerequisites
 - R (>= 3.4.0)
 - R packages: Rcpp (>= 0.12.10), RcppArmadillo (>= 0.7.700.3.0), ggplot2 (>= 2.1.0),
 rtdists (>= 0.6-6), coda (>= 0.16-1) tmvtnorm, matrixStats,
 data.table
 - Windows users need Rtools (>= 3.3.0.1959) 
 - Mac OS users need to make clang understand OpenMP flag
 - Linux/Unix users may need to install Open MPI library, if it has not 
   been installed. 
 - [Armadillo](https://CRAN.R-project.org/package=RcppArmadillo)
   may need a recent g++ compiler > 4.6

## Installation

From CRAN: 
> install.packages("ggdmc")

From source: 

> install.packages("ggdmc_0.2.5.5.tar.gz", repos = NULL, type="source")

From GitHub (you will need to install _devtools_:

> devtools::install_github(“yxlin/ggdmc”)

For Mac Users:

1. You will need to install [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS).
As to 27, Aug, 2018, please install gfortran 6.1, even you are using a macOS 
High Sierra Version 10.13.4. gfortran 6.3 may not work.

2. Install clang4-r. 
[James Balamuta](https://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/)
has created a convenient tool, [clang4-r](https://uofi.app.box.com/v/r-macos-clang-pkg).
Once you install clang4-r, your clang will then understand the OpenMP flag
in _ggdmc_. The aim is to allow macOS to understand OpenMP flag, so you may use 
other methods for that purpose, if you do not want to install clang4-r. The
clang4-r is the most straightforward we found so far. 
However we have not looked into the source code of clang4-r. Use it at your
own risk.

## How to install many supporting packages
One common hurlde in install a large R package is to install many supporting
packages.  One useful method to mitigate such hurdle is use a convenient
command line script. 

This works in Unix-like OS. First of all, you need the littler package.

> sudo apt-get install r-cran-littler

Second copy-and-paste the following script into a file, called it _install.r_.
```
#!/usr/bin/env r

if (is.null(argv) | length(argv)<1) {
  cat("Usage: installr.r pkg1 [pkg2 pkg3 ...]\n")
  q()
}

repos <- "http://cran.rstudio.com"

lib.loc <- "/usr/local/lib/R/site-library"

install.packages(argv, lib.loc, repos)

```

(Third, add x mode to install.r)
> chmod +x install.r

Finally, at the command line, enter the following:
> install.r Rcpp RcppArmadillo ggplot2 ggmcmc coda tmvtnorm matrixStats data.table

This shall install all supporting packages at once. 
This useful method is from Dirk Eddelbuettel's [littler]<http://dirk.eddelbuettel.com/code/littler.html>.

## Citation

If you use this package, please cite the software, for example:

Lin, Y.-S. (in preparation). Tutorial on Bayesian cognitive modeling. 


## Contributors
The R documentation, tutorials, C++ codes, parallel computation in OpenMP, new
genetic algorithm, and R helper functions and R packaging are 
developed by Yi-Shin Lin. DMC is developed 
by Andrew Heathcote (Heathcote et al., 2018).

_ggdmc_ evolves from DMC.  Although these two tools share some similarities, 
they have some differences.  They are designed to tackle different problems, 
although some may overlap.  _ggdmc_ implements almost all time-critical
Bayesian routines in C++, but uses similar attribute naming system, so 
samples drawn from DMC may be piped into ggdmc C++ samplers.  Please report bugs
to [me](mailto:yishinlin001@gmail.com).

## License

GPL-2 

## Acknowledgments

* The PDF, CDF and random deviates generation of DDM are derived from 
Voss & Voss's fast-dm 30.2 and rtdists 0.9-0. 
* Truncated normal functions are originally based on 
[Jonathan Olmsted's](mailto:jpolmsted@gmail.com) RcppTN 0.1-8 at
https://github.com/olmjo/RcppTN,
[Christopher Jackson's](chris.jackson@mrc-bsu.cam.ac.uk) R codes in msm package,
and Robert (1995, Statistics & Computing). 
* Armadillo is a collection of C++ library, conducting linear
algebra <http://arma.sourceforge.net/>. 
* Thanks to Matthew Gretton's sharing his know-how in rtdists. 
* Thanks to Andrew Heathcote for lending me his MacBook Air. 
_ggdmc_ works on OS X (macOS High Sierra Version 10.13.4) 
