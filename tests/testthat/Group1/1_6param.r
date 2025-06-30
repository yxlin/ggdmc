# q(save = "no")
cat("\n\n-------------------- Run 6 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

wkdir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(wkdir, "data/lba_data1.rda")
load(fn)


model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = "M", sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "lba"
)

dmis <- ggdmcModel::BuildDMI(dat, model)
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

priors <- set_priors(p_prior = p_prior)
# nchain <- model@npar * 3
theta_input <- setThetaInput(nmc = 500, pnames = model@pnames, thin = 1)
samples0 <- ggdmc::initialise_theta(theta_input, priors, dmis[[1]],
    seed = 846671,
    verbose = TRUE
)

de_input <- ggdmc::setDEInput(
    sub_migration_prob = 0.05,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)
configs <- ggdmc::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]


fit0 <- run_subject(config_r = cfg, dmi = dmis[[1]], samples = samples0)


theta_input <- setThetaInput(nmc = 500, pnames = model@pnames, thin = 1)
de_input <- ggdmc::setDEInput(
    sub_migration_prob = 0.00,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)
configs <- ggdmcDE::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]
fit1 <- run_subject(cfg, dmis[[1]], fit0)

nmc <- 1000
theta_input <- setThetaInput(nmc = nmc, pnames = model@pnames, thin = 4)
de_input <- ggdmcDE::setDEInput(
    sub_migration_prob = 0.00,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)

configs <- ggdmc::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]
fit2 <- run_subject(cfg, dmis[[1]], fit1)


fit <- fit2
plot(fit)
plot(fit, den = T, pll = F)

options(digits = 2)
est0 <- compare(fit, ps = p_vector, verbose = TRUE)
#  A      B mean_v.false mean_v.true sd_v.true      t0
# True           0.7500 1.2500       1.5000      2.5000    0.1000  0.1500
# 2.5% Estimate  0.5679 0.9448       1.3031      2.1957    0.0749  0.0050
# 50% Estimate   1.1548 2.1866       3.1034      3.8705    0.1842  0.0895
# 97.5% Estimate 2.1641 3.8581       5.6637      6.4482    0.4239  0.2774
# Median-True    0.4048 0.9366       1.6034      1.3705    0.0842 -0.0605
