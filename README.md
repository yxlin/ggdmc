# Modelling Cognitive Processes

_ggdmc_ is an R package for modelling cognitive processes. Although its focus
is on the challenging hierarchical Bayesian models, fitting them with Bayesian
MCMC. It can also fit cognitive models with conventional methods, such as
maximum likelihood estimation and least squares. The package uses the sampling
method of population-based Markov chain Monte Carlo (pMCMC).

## Getting Started
This example demonstrates the Wiener diffusion model.  For other models, 
see my [tutorials site](https://yxlin.github.io/). The naming of _R_ functions 
in _ggdmc_ attempts to inform the user what the functions are for. For 
example,  _BuildModel_ is to build a model object.  

As the user is often reminded in using Bayesian tools, it is always a good 
practice to check the result of a model fit. Note that the sequence of 
parameters in a parameter vector (i.e., p.vector) must follow the sequence in 
the _p.vector_ reported by _BuildModel_.  Some built-in checks will try to 
safeguard this, but they are far from bulletproof. 

## Fit a fixed-effect model to a participant

```
## Set up model ----
## Fixing sv & sz to 0 to set up a Wiener diffusion model
require(ggdmc)
model <- BuildModel(
  p.map     = list(a = "1", v="1", z="1", d="1", sz="1", sv="1", t0="1", 
                   st0="1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1","r2"),
  constants = c(st0 = 0, d = 0, sv = 0, sz = 0),  
  type      = "rd")   

npar <- model@npar   ## Note this works for version > 0.2.7.5; 
## npar <- length(GetPNames(model))   ## Use GetPNames instead in 0.2.6.0

p.vector <- c(a=1, v=1.5, z=0.5, t0=.15)
dat <- simulate(model, nsim = 50, ps = p.vector)
dmi <- BuildDMI(dat, model)

p.prior <- BuildPrior(
  dists = rep("tnorm", npar),
  p1=c(a=1, v=0, z=1, t0=1),
  p2=c(a=1, v=2, z=1, t0=1),
  lower = c(0, -5, rep(0, 2)),
  upper = rep(NA, npar))

## Fit model -------------
fit0 <- StartNewsamples(dmi, p.prior)
fit  <- run(fit0)

## Check model -----------
plot(fit)
plot(fit, den = TRUE)
plot(fit, pll = FALSE)
plot(fit, pll = FALSE, den = TRUE)

isconv <- gelman(fit)
est    <- summary(fit, recovery = TRUE, ps = p.vector, verbose = TRUE)

```

## How to fit fixed-effect and hierarchical model with multiple participants

```
require(ggdmc);
model <- BuildModel(
  p.map     = list(a = "1", v ="1", z ="1", d ="1", sz ="1", sv ="1", t0 ="1", 
                   st0 ="1"),
  match.map = list(M = list(s1 = "r1", s2 = "r2")),
  factors   = list(S = c("s1", "s2")),
  responses = c("r1","r2"),
  constants = c(st0 = 0, d = 0, sv = 0, sz = 0),
  type      = "rd")

npar <- model@npar
pop.mean  <- c(a = 2,   v = 4,  z = 0.5, t0 = 0.3)
pop.scale <- c(a = 0.5, v = .5, z = 0.1, t0 = 0.05)
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

plot(p.prior, ps = ps)  ## Check if all true values are in the range 

## Sampling separately
fit0 <- StartNewsamples(dmi, p.prior, ncore = 4)
fit  <- run(fit0, 5e2, ncore = 4)
fit  <- run(fit, 1e2, add = TRUE, ncore = 4)  ## add additional 100 samples

## Check model -----
isconv <- gelman(fit, verbose = TRUE)
plot(fit)
est0 <- summary(fit, recovery = TRUE, ps = ps, verbose = TRUE)

## Sampling hierarchically
mu.prior <- BuildPrior(
    dists = rep("tnorm", npar),
    p1    = pop.mean,
    p2    = pop.scale*5,
    lower = c(0,-5,  0, 0),
    upper = c(5, 7,  1, 1))

sigma.prior <- BuildPrior(
    dists = rep("beta", npar),
    p1    = c(a=1, v=1, z=1, t0=1),
    p2    = rep(1, npar),
    upper = rep(1, npar))

## !!!The names are important!!!
priors <- list(pprior = p.prior, location = mu.prior, scale = sigma.prior)
names(priors)
## [1] "pprior"   "location" "scale"

## Fit hierarchical model ----
fit0 <- StartNewsamples(dmi, priors)
fit  <- run(fit0, 5e2)

p0 <- plot(fit, hyper = TRUE)
p0 <- plot(fit, hyper = TRUE, den = TRUE, pll=FALSE)

## Check model -----------
## hgelman function is deprecated 
res  <- gelman(fit, verbose = TRUE)
est0 <- summary(fit, recovery = TRUE, ps = ps, verbose = TRUE)
est1 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.mean,  type = 1, verbose = TRUE)
est2 <- summary(fit, hyper = TRUE, recovery = TRUE, ps = pop.scale, type = 2, verbose = TRUE)

```

## Response time models 
1. The LBA model, type = "norm",
2. The DDM, type = "rd",
3. The Wiener diffusion, type = "rd" and set sv=0 and sz=0

## PDA-based models 
4. The Piecewise LBA model 0; CPU-based PDA likelihoods; type = "plba0",
5. The Piecewise LBA model 1; CPU-based PDA likelihoods; type = "plba1", 
6. The Piecewise LBA model 0; GPU-based PDA likelihoods; type = "plba0_gpu", 
7. The Piecewise LBA model 1; GPU-based PDA likelihoods; type = "plba1_gpu", 
8. The LBA model; GPU-based PDA likelihoods;, type = "norm_pda_gpu",
9. The leaky, competing accumulator model (Experimental!).

4 to 8 are separated from the latest version of the package. For these 
PDA-based models, see my BRM paper and associated packages there (osf project). 
9 is in a separate module, which has yet to be incorporated. See the LCA 
[tutorial](https://yxlin.github.io/cognitive-model/lca/) for its testing result
using MLE. 

For the details regarding PLBA types, please see 
[Holmes, Trueblood, and Heathcote (2016)](http://dx.doi.org/10.1016/j.cogpsych.2015.11.002)

## Experimental (untested) models 
10. 2-D/circular drift-diffusion model, type = "cddm"
11. Prospective memory model, type = "norm" (see tutorial for more details)
12. Time-varying changes in other free parameters

## Further information
One aim in designing _ggdmc_ is to read objects from DMC, which share 
similarities. They have, however, some differences. For example, in the latest
version of _ggdmc_, the dimension of theta and phi arrays are 
'npar x nchain x nmc'. DMC uses 'nchain x npar x nmc'. To reduce the 
computation time for manipulating the matrices and arrays, we change this to accommodate the Armadillo convention
change to accommodate the convention in Armadillo. Similarly, the dimension of 
the 'log_likelihoods' and 'summed_log_prior' 
matrices are 'nchain x nmc'. DMC uses 'nmc x nchain'.  Remember to transpose
them if you want to operate objects back and forth. Currently, we use
the R functions, 'aperm' and 't', to transpose matrices and arrays when operating
have to operate between DMC and _ggdmc_. The following two convenient functions
are designed for doing this operation. 

```
DMC2ggdmc <- function(x) {
  ## x is an object of posterior samples from individual subject fit
  x$theta <- aperm(x$theta, c(2, 1, 3))
  x$summed_log_prior <- t(x$summed_log_prior)
  x$log_likelihoods <- t(x$log_likelihoods)
  class(x) <- c('list', 'model')
  return(x)
}
ggdmc2DMC <- function(x) {
  ## Should change $ to @ whenoperatinge output from ggdmc's run function,
  ## because ggdmc now uses S4 class
  x$theta <- aperm(x$theta, c(2, 1, 3))
  x$summed_log_prior <- t(x$summed_log_prior)
  x$log_likelihoods <- t(x$log_likelihoods)
  return(x)
}
```

Note **Dstats.dmc** in DMC is also affected by the issue of the different array and 
matrix dimensions because Dstats.dmc calculates the means of the theta/phi  
array across the column, 

```
apply(samples$theta,2,mean)

```
_ggdmc_ provides DIC function, which uses a back-end function, 
**deviance_model** to attain the same operation.

The tutorial in [3-accumulator LBA model](https://yxlin.github.io/lba3) 
illustrates an example of doing the back-and-forth operation.

Note that we start to use S4 class after version 0.2.7.5, so switch to use "@" 
operator to extract object components (i.e., slot). 

## Prerequisites
 - R (>= 3.3.0)
 - R packages: Rcpp (>= 0.12.10), RcppArmadillo (>= 0.7.700.3.0), 
   ggplot2 (>= 2.1.0), coda (>= 0.16-1), matrixStats, data.table
 - Windows users need Rtools (>= 3.3.0.1959) 
 - ~~Mac OS users need to make Clang understand the OpenMP flag.~~
 - ~~Linux/Unix users may need to install the Open MPI library if it has not 
   been installed.~~ 
 - ~~[Armadillo](https://CRAN.R-project.org/package=RcppArmadillo)
   may need a recent g++ compiler > 4.6~~

## Installation

We now use S4 class after version 0.2.7.5. The new design enables a more 
user-friendly interface in help pages.

From CRAN (0.2.6.0): 
> install.packages("ggdmc")

From source: 

> install.packages("ggdmc_0.2.8.0.tar.gz", repos = NULL, type="source")

From GitHub (need _devtools_) (0.2.8.0):

> devtools::install_github("yxlin/ggdmc")

For Microsoft R users:

As of 06-01-2020, because Microsoft R uses R version 3.5.3, the user who wishes 
to deploy ggdmc on Microsoft R may encounter two challenges. First is 
RcppArmadillo on MRAN is behind the one on R CRAN. The RcppArmadillo on MRAN 
has yet to introduce recent Armadillo functions, for instance, randperm in C++. 
This can be resolved by installing RcppArmadillo directly from its source 
tarball, downloaded from CRAN. Secondly, the default installation process on 
Windows is to look for the package binary matching the R version on a Windows 
machine. This may result in Microsoft R looking for a version of ggdmc matching 
R 3.5.3, and thereby, it cannot find one. This can be resolved similarly by 
installing from the source tarball. 

For Mac Users:

~~1. Install [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries#MacOS).
As of 27, Aug, 2018, the gfortran version has to be 6.1, even if you are using a 
macOS High Sierra Version 10.13.4. gfortran 6.3 may not work.~~

~~2. Install clang4-r. 
[James Balamuta](https://thecoatlessprofessor.com/programming/openmp-in-r-on-os-x/)
has created a convenient tool, [clang4-r](https://uofi.app.box.com/v/r-macos-clang-pkg).
Once you install clang4-r, your clang will understand the OpenMP flag
in _ggdmc_. The aim is to allow macOS to understand the OpenMP flag, so you may use 
other methods if you do not want to install clang4-r. The
clang4-r is the most straightforward we have found so far. 
However, we have yet to look into the source code of clang4-r. You can use it at your
own risk.~~

The configure script now disables OpenMP, so macOS users should be able to 
install without encountering the OpenMP problem. 


FAQ:
1. When the error message arises, "/usr/bin/ld: cannot find -lgsl" and/or 
"/usr/bin/ld: cannot find -lgslcblas", installing libgsl-dev may resolve this 
problem (Ubuntu).

## Citation
Lin, Y.-S and Strickland, L. (2020). Evidence accumulation models with R: A 
a practical guide to hierarchical Bayesian methods. The Quantitative Methods for Psychology.

## Contributors
The R documentation, tutorials, C++ codes, parallel computations, new
genetic algorithm, R helper functions and R packaging are 
developed by Yi-Shin Lin. A substantial part of R codes for handling 
experimental designs are adapted from the DMC, developed by Andrew Heathcote 
(Heathcote et al., 2018). You could find different and more interesting 
cognitive models in DMC. 

Please report bugs to [me](mailto:yishinlin001@gmail.com) or start an issue
here.

## Correction
The help page for the function, _likelihood_ in Density.cpp states that 
it returns log-likelihood (v2.8.0). An inspection of 
the source code found that it returns likelihood, not log likelihood. (13-06-2022; v2.8.1). 
Thanks for Nachshon Meiran pointing it out.

## License

[GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt)

## Acknowledgments

* The PDF, CDF and random number generation of DDM were derived from 
Voss & Voss's fast-dm 30.2 and rtdists 0.9-0. 
* Truncated normal functions were originally based on 
[Jonathan Olmsted's](mailto:jpolmsted@gmail.com) RcppTN 0.1-8 at
https://github.com/olmjo/RcppTN,
[Christopher Jackson's](chris.jackson@mrc-bsu.cam.ac.uk) R codes in msm package,
and Robert's paper (1995, Statistics & Computing). 
* Thanks to Matthew Gretton's consultation regarding the rtdists package. 
* Thanks to Andrew Heathcote for lending me his MacBook Air. _ggdmc_ works on 
OS X (macOS High Sierra Version 10.13.4) 
* The PDF and random number generation of the 2-D diffusion/circular diffusion 
model is based on Smith (2016).

## Reference
* Heathcote, A., Lin, Y.-S., Reynolds, A., Strickland. L. Gretton, M., & Matzke, D. (2018). Dynamic models of choice, _Behavior Research Methods_. https://doi.org/10.3758/s13428-018-1067-y
* Lin, Y.-S. and Strickland, L. (2020). [Evidence accumulation models with R: A practical guide to hierarchical Bayesian methods](https://www.tascl.org/uploads/4/9/3/3/49339445/lin_strickland_2019.pdf). _The Quantitative Methods for Psychology_.
* Smith, P. (2016). Diffusion Theory of Decision Making in Continuous Report, _Psychological Review_, 123(4), 425-451. http://dx.doi.org/10.1037/rev0000023
