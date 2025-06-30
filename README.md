# ggdmc
ggdmc provides tools to conduct Bayesian inference on a range of choice response time models.

# Getting Started

## Installation

From CRAN:
```
install.packages("ggdmc")
```

## Examples
This is the core package for doing Bayesian inference on a range of choice response time models. 
The following example shows a complex LBA model design and its recovery study. 

```
cat("\n\n------------A LBA B x v model-------------")
rm(list = ls())
pkg <- c("lbaModel", "ggdmcPrior", "ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ggdmc/tests/testthat/Group1/data/"

helper_path <- paste0(wkdir, "helpers.r")
save_path <- paste0(wkdir, "lba_data6.rda")
source(helper_path)


model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = c("S", "COLOR"), t0 = "1", mean_v = c("NOISE", "M"), sd_v = "M", st0 = "1"),
    match_map = list(M = list(left = "z_key", right = "x_key")),
    factors = list(
        S = c("left", "right"),
        COLOR = c("red", "blue"),
        NOISE = c("high", "moderate", "low")
    ),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("z_key", "x_key"),
    type = "lba"
)

pop_mean <- c(
    A = 0.25,
    B.left.red = 2.33,
    B.right.red = 2.13,
    B.left.blue = 2.55,
    B.right.blue = 2.35,
    t0 = 0.05,
    mean_v.high.true = 2.1,
    mean_v.moderate.true = 3.1, mean_v.low.true = 4.1, mean_v.high.false = 1.82,
    mean_v.moderate.false = 1.70, mean_v.low.false = 1.5, sd_v.true = 0.2
)

pop_scale <- c(
    A = 0.2,
    B.left.red = .23,
    B.right.red = .23,
    B.left.blue = .255,
    B.right.blue = .235,
    t0 = 0.05,
    mean_v.high.true = .25,
    mean_v.moderate.true = .35, mean_v.low.true = .45, mean_v.high.false = .182,
    mean_v.moderate.false = .170, mean_v.low.false = .15, sd_v.true = 0.2
)


pop_dist <- ggdmcPrior::BuildPrior(
    dists = rep("tnorm", model@npar),
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    log_p = rep(FALSE, model@npar)
)

p_vector <- c(
    A = 0.25,
    B.left.red = 2.33,
    B.right.red = 2.13,
    B.left.blue = 2.55,
    B.right.blue = 2.35,
    t0 = 0.05,
    mean_v.high.true = 2.1,
    mean_v.moderate.true = 3.1, mean_v.low.true = 4.1, mean_v.high.false = 1.82,
    mean_v.moderate.false = 1.70, mean_v.low.false = 1.5, sd_v.true = 0.2
)

# ---------------------------------------
sub_model <- lbaModel::setLBA(model)
pop_model <- lbaModel::setLBA(model, population_distribution = pop_dist)

dat <- simulate(sub_model, nsim = 768, parameter_vector = p_vector, n_subject = 1)
hdat <- simulate(pop_model, nsim = 768, n_subject = 32)

sub_dmis <- ggdmcModel::BuildDMI(dat, model)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)
hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)

<!-- res <- zx_get_accuracy(dat)
res <- zx_get_accuracy(hdat) -->
ps <- attr(hdat, "parameters")

true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)

# Generate subj samples -------------------------
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

sub_priors <- set_priors(p_prior = p_prior)


sub_theta_input <- setThetaInput(pnames = model@pnames)
sub_samples <- initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]])

save(hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi, sub_priors, sub_samples, sub_theta_input,
    file = save_path
)

# Generate pop samples -------------------------
p0 <- runif(model@npar)
names(p0) <- model@pnames

model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("tnorm", model@npar),
    log_p = rep(TRUE, model@npar)
)

# Prior log likelihoods
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
l_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
s_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

pop_priors <- ggdmcPrior::set_priors(p_prior = model_likelihood, l_prior = l_prior, s_prior = s_prior)

pop_theta_input <- ggdmc::setThetaInput(pnames = pop_priors@pnames)
pop_samples <- initialise_phi(pop_theta_input, pop_priors, pop_dmis)


save(hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi,
    sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input,
    file = save_path
)

# Sampling
fits0 <- StartSampling(pop_dmis, pop_priors, sub_migration_prob = 0.06, thin = 8L, seed = 9032)
save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    sub_migration_prob = 0.00, pop_migration_prob = 0.05,
    nmc = 1000L, report_length = 200L,
    thin = 8L, seed = 9032
)
save(fits0, fits1, file = save_path)

fits2 <- RestartSampling(fits1,
    sub_migration_prob = 0.00, pop_migration_prob = 0.01,
    nmc = 1000L, report_length = 200L,
    thin = 8L, seed = 9032
)
save(fits0, fits1, fits2, file = save_path)



fits <- fits2
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)
est_phi <- compare(phi, ps = true_vector)
est_theta <- compare_many(thetas, ps = ps)

rhat <- gelman(phi)
plot(phi, pll = F, den = T)

DT <- prepare_thetas_data(thetas, start = 5000)
DT1 <- prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
DT2 <- prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[2]]$phi@nmc * 0.5)
DT3 <- prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[3]]$phi@nmc * 0.5)

p1 <- plot_thetas(DT)
p1 <- plot_thetas(DT1, start = 300, end = 400)
p1 <- plot_thetas(DT1, start = 300, end = 400, subjects = 5)
p1 <- plot_thetas(DT1, start = 300, end = 400, subjects = as.character(1:10))
p1 <- plot_thetas(DT1, start = 300, end = 400, max_subjects = 8)


```


# Prerequisites
R (>= 3.3.0), Rcpp (>= 1.0.7), RcppArmadillo (>= 0.10.7.5.0), ggdmcHeaders, ggdmcPrior, ggdmcLikelihood, data.table, matrixStats,lattice.

See DESCRIPTION for details


## Citation
Lin, Y.-S and Strickland, L. (2020). Evidence accumulation models with R: A 
a practical guide to hierarchical Bayesian methods. The Quantitative Methods for Psychology.

## Contributors
The early version of the ggdmc was developed from the Dynamic Model of Choice, (Heathcote et al., 2018). 

Please report bugs to [me](mailto:yishinlin001@gmail.com) or start an issue here.

# More Examples
The second example shows a minimal DDM design and its recovery study. More variants should
be updated at the [tutorials site](https://yxlin.github.io/). 

```
pkg <- c("ddModel", "ggdmcPrior")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/data/"
save_path <- paste0(wkdir, "ddm_data0.rda")

hyper_model <- ggdmcModel::BuildModel(
    p_map = list(
        a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
        t0 = "1", st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
    accumulators = c("r1", "r2"),
    type = "hyper",
    verbose = FALSE
)

model <- ggdmcModel::BuildModel(
    p_map = list(
        a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
        t0 = "1", st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
    accumulators = c("r1", "r2"),
    type = "fastdm"
)


pop_mean <- c(a = 1, sz = 0.25, t0 = 0.15, v = 2.5, z = 0.38)
pop_scale <- c(a = 0.05, sz = 0.01, t0 = 0.02, v = .5, z = 0.01)
pop_dist <- ggdmcPrior::BuildPrior(
    p0    = pop_mean,
    p1    = pop_scale,
    lower = c(0, 0, 0, -10, 0),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)


# ---------------------------------------
sub_model <- setDDM(model)
pop_model <- setDDM(model, population_distribution = pop_dist)

# Note the package asks the user to enter absolute starting point (z & sz),
# not the relative startning point (zr & szr)
# (ie starting point and its variability relative to a).
#
# The internal C++ backend will calculate zr and szr for you.
# The advantage of this is that the user can now set wide range of 
# values for the z and sz. 
#
# In the old version, when the user enters the zr and szr that is
# not within 0 and 1, the sampler will return low density.

p_vector <- c(a = 1, sz = 0.25, t0 = 0.15, v = 2.5, z = .38)

dat <- simulate(sub_model, nsim = 256, parameter_vector = p_vector, n_subject = 1, debug = TRUE)
hdat <- simulate(pop_model, nsim = 128, n_subject = 32)

# You may change the simulate time resolution.
# dat <- simulate(sub_model,
#     nsim = 256, parameter_vector = p_vector, n_subject = 1,
#     time_parameters = c(t_min = -0.5, tmax = 0.5, dt = 0.01),
#     debug = TRUE
# )
# hdat <- simulate(pop_model, nsim = 128, n_subject = 32, 
# time_parameters = c(t_min = -0.9, tmax = 0.9, dt = 0.01))

sub_dmis <- ggdmcModel::BuildDMI(dat, model)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)
hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)

options(digits = 2)
cat("Accuracy: \n")
c(mean(dat$C), mean(hdat$C))

ps <- attr(hdat, "parameters")
true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)


# Generate subj samples -------------------------
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

sub_priors <- set_priors(p_prior = p_prior)

nmc <- 500
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames)
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, 
sub_dmis[[1]], seed = 846671, verbose = F)

save(hyper_model, model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi, sub_priors, sub_samples, sub_theta_input,
    file = save_path
)


# Generate pop samples -------------------------
p0 <- runif(model@npar)
names(p0) <- model@pnames
model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("tnorm", model@npar),
    log_p = rep(TRUE, model@npar)
)



# Prior log likelihoods
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
l_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)


s_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

pop_priors <- ggdmcPrior::set_priors(p_prior = model_likelihood, l_prior = l_prior, s_prior = s_prior)
pop_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = pop_priors@pnames)

pop_samples <- ggdmc::initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671)

save(hyper_model, model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi,
    sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input,
    file = save_path
)


```

The following show the parameter recovery study.
```
fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.05,
    thin = 8L, pop_debug = F, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    pop_migration_prob = 0.02,
    sub_migration_prob = 0.00,
    thin = 4L, seed = 9032
)
save(fits0, fits1, file = save_path)

fits2 <- RestartSampling(fits1,
    pop_migration_prob = 0.01,
    sub_migration_prob = 0.00,
    thin = 2L, seed = 9032
)

fits3 <- RestartSampling(fits2,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    thin = 2L, seed = 9032
)

fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

est_phi <- summary(phi)
est_phi <- compare(phi, ps = true_vector)

post_summary <- summary(phi, start = 501) # Discard first 500 as burn-in
# Custom quantiles
detailed_summary <- summary(phi, probability = seq(0.1, 0.9, by = 0.1))
detailed_summary$quantiles

# Extract specific elements
posterior_means <- post_summary$statistics[, "Mean"]
credible_intervals <- post_summary$quantiles[, c("5%", "97.5%")]
result <- summary_many(thetas)
result <- summary_many(thetas, verbose = T)

options(digits = 2)
result <- compare_many(thetas, ps = ps)
result <- compare_many(thetas, ps = ps, verbose = TRUE)



# The PSRF calculate is more straightforward.
hat <- gelman(phi)

hat <- gelman(fits[[1]]$phi)
hat <- gelman(fits[[2]]$phi)
hat <- gelman(fits[[3]]$phi)

DT <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
DT <- ggdmc::prepare_thetas_data(thetas, start = 5000)

p1 <- plot_thetas(DT)
p1 <- plot_thetas(DT, start = 300, end = 400)
p1 <- plot_thetas(DT, start = 300, end = 400, subjects = 5)
p1 <- plot_thetas(DT, start = 300, end = 400, subjects = as.character(1:10))
p1 <- plot_thetas(DT, start = 300, end = 400, max_subjects = 8)

```

## Acknowledgments
* The early version of the PDF, CDF, and random number generation of DDM were derived from 
Voss & Voss's fast-dm 30.2 and rtdists 0.9-0. They have been rewritten in the latest version. 

* Truncated normal functions were originally based on 
[Jonathan Olmsted's](mailto:jpolmsted@gmail.com) RcppTN 0.1-8 at
https://github.com/olmjo/RcppTN,
[Christopher Jackson's](chris.jackson@mrc-bsu.cam.ac.uk) R codes in msm package,
and Robert's paper (1995, Statistics & Computing). 

## Reference
* Heathcote, A., Lin, Y.-S., Reynolds, A., Strickland. L. Gretton, M., & Matzke, D. (2018). Dynamic models of choice, _Behavior Research Methods_. https://doi.org/10.3758/s13428-018-1067-y
* Lin, Y.-S. and Strickland, L. (2020). [Evidence accumulation models with R: A practical guide to hierarchical Bayesian methods]([https://www.tascl.org/uploads/4/9/3/3/49339445/lin_strickland_2019.pdf](https://www.tqmp.org/RegularArticles/vol16-2/p133/p133.pdf)). _The Quantitative Methods for Psychology_.

