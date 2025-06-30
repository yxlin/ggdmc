# q(save = "no")
cat("\n\n--------------------Testing Drift Rate Model--------------------")
rm(list = ls())
pkg <- c("ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

wkdir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(wkdir, "data/lba_data2.rda")
load(fn)


model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = c("POLITICAL_VIEW", "M"), sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2"), POLITICAL_VIEW = c("liberal", "conservative")),
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
nmc <- 500
nchain <- model@npar * 3
thin <- 1
theta_input <- setThetaInput(nmc = nmc, nchain = nchain, pnames = model@pnames, thin = thin)
samples0 <- initialise_theta(theta_input, priors, dmis[[1]], seed = 929726)

samples0@pnames
samples0@npar

# Burn-in ----------------------
de_input <- ggdmc::setDEInput(
    sub_migration_prob = 0.05,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)
configs <- ggdmc::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]

fit0 <- run_subject(config_r = cfg, dmi = dmis[[1]], samples = samples0)


# Sampling ----------------------
de_input <- ggdmc::setDEInput(
    sub_migration_prob = 0.00,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)

configs <- ggdmcDE::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]
fit1 <- run_subject(cfg, dmis[[1]], fit0)
fit <- fit1

options(digits = 2)
est0 <- compare(fit, start = 1, ps = p_vector[model@pnames], verbose = TRUE)
#                     A      B mean_v.conservative.false mean_v.conservative.true
# True           0.7500 1.2500                    1.5000                   2.1500
# 2.5% Estimate  0.5850 0.9095                    0.7691                   1.6877
# 50% Estimate   0.8307 1.6002                    1.7131                   2.4229
# 97.5% Estimate 1.0847 2.2294                    2.4786                   3.1535
# Median-True    0.0807 0.3502                    0.2131                   0.2729
#                mean_v.liberal.false mean_v.liberal.true sd_v.true      t0
# True                         1.5000              2.5000    0.1000  0.1500
# 2.5% Estimate                0.7816              1.9778    0.0713  0.0032
# 50% Estimate                 1.7272              2.7725    0.1030  0.0705
# 97.5% Estimate               2.5679              3.5635    0.1439  0.2130
# Median-True                  0.2272              0.2725    0.0030 -0.0795

#                     A      B mean_v.conservative.false mean_v.conservative.true
# True           0.7500 1.2500                    1.5000                   2.1500
# 2.5% Estimate  0.5426 0.8538                    0.7223                   1.6010
# 50% Estimate   0.8056 1.5409                    1.6352                   2.3501
# 97.5% Estimate 1.0500 2.1288                    2.3897                   3.0325
# Median-True    0.0556 0.2909                    0.1352                   0.2001
#                mean_v.liberal.false mean_v.liberal.true sd_v.true      t0
# True                         1.5000              2.5000    0.1000  0.1500
# 2.5% Estimate                0.6389              1.8439    0.0710  0.0032
# 50% Estimate                 1.6453              2.6911    0.1009  0.0827
# 97.5% Estimate               2.4415              3.4350    0.1370  0.2150
# Median-True                  0.1453              0.1911    0.0009 -0.0673
