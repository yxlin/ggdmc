# q(save = "no")
cat("\n\n--------------------Testing B Model--------------------")
rm(list = ls())
pkg <- c("ggdmc")

suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

wkdir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(wkdir, "data/lba_data3.rda")
load(fn)

model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = c("S", "政黨傾向"), mean_v = "M", sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(S = c("紅", "黃", "藍", "綠"), 政黨傾向 = c("自由派", "保守派")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
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

priors <- ggdmcPrior::set_priors(p_prior = p_prior)
nmc <- 500
nchain <- model@npar * 3
thin <- 4
theta_input <- setThetaInput(nmc = nmc, nchain = nchain, pnames = model@pnames, thin = thin)
samples0 <- initialise_theta(theta_input, priors, dmis[[1]], seed = 929726)


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
configs <- ggdmc::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]
fit1 <- run_subject(cfg, dmis[[1]], fit0)
fit <- fit1

options(digits = 2)
est0 <- compare(fit, start = 1, ps = p_vector[model@pnames], verbose = TRUE)
#                      A    B.紅.保守派 B.紅.自由派   B.綠.保守派 B.綠.自由派
# True            0.2500      2.4000      1.2000      4.2000      2.1000
# 2.5% Estimate   0.0047      1.6774      0.6351      3.1318      1.4086
# 50% Estimate    0.1233      2.4144      1.2296      4.1593      2.0982
# 97.5% Estimate  0.6184      3.2029      1.7997      5.1885      2.8129
# Median-True    -0.1267      0.0144      0.0296     -0.0407     -0.0018
#                  B.藍.保守派   B.藍.自由派  B.黃.保守派  B.黃.自由派 mean_v.false
# True                3.6000      1.8000      3.0000      1.5000       1.1500
# 2.5% Estimate       2.6822      1.1314      2.2000      0.9273       0.6324
# 50% Estimate        3.5775      1.8189      3.0532      1.5497       1.1933
# 97.5% Estimate      4.5049      2.4884      3.9173      2.1907       1.7097
# Median-True        -0.0225      0.0189      0.0532      0.0497       0.0433
#                mean_v.true sd_v.true      t0
# True                2.8000    0.8000  0.1000
# 2.5% Estimate       2.3742    0.6575  0.0072
# 50% Estimate        2.8530    0.7684  0.0968
# 97.5% Estimate      3.3391    0.8635  0.2022
# Median-True         0.0530   -0.0316 -0.0032


# expected_values <- "tests/testthat/Group0/data/expected_model3.rda"
expected_values <- "data/expected_model3.rda"
# samples <- samples0
# est <- est0
# save(samples, est, file = expected_values)

load(expected_values)
testthat::expect_true(all(samples0@theta[, , 1] == samples@theta[, , 1]))
testthat::expect_true(all(est0 == est))
