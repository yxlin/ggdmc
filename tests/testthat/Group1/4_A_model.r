# q(save = "no")
cat("\n\n--------------------Testing A Model--------------------")
rm(list = ls())
pkg <- c("ggdmcDE", "ggdmcModel", "ggdmcPrior", "ggdmcPhi", "ggdmcLikelihood", "ggplot2")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
fn <- "tests/testthat/Group1/data/gdmc_data4.rda"
# fn <- "data/gdmc_A_model.RData"
load(fn)
model <- ggdmcModel::BuildModel(
    p_map = list(
        B = "1", A = c("S", "政黨傾向"), mean_v = "M", sd_v = "M", st0 = "1",
        t0 = "1"
    ),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(
        S = c("紅", "黃", "藍", "綠"),
        政黨傾向 = c("自由派", "保守派")
    ),
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
priors <- set_priors(p_prior = p_prior)

nmc <- 500
nchain <- model@npar * 3
thin <- 4
theta_input <- setThetaInput(nmc = nmc, nchain = nchain, pnames = model@pnames, thin = thin)
samples0 <- initialise_theta(theta_input, priors, dmis[[1]], seed = 929726)

# Burn-in ----------------------
de_input <- ggdmcDE::setDEInput(
    migration_prob_Hu = 0.05,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)
configs <- ggdmcDE::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]
fit0 <- run_subject(cfg, dmis[[1]], samples0, seed = 123)

# gdmc:::plot(fit0, start = 100)
# pdf("tests/testthat/Group1/Rplot.pdf")
# p0 <- gdmc:::plot(fit0, pll = F, start = 300)
# p1 <- gdmc:::plot(fit0, pll = F, den = TRUE, start = 300)
# dev.off()

# Sampling ----------------------
de_input <- ggdmcDE::setDEInput(
    migration_prob_Hu = 0.00,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)

configs <- ggdmcDE::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]
fit1 <- run_subject(cfg, dmis[[1]], fit0, seed = 123)

est0 <- gdmc:::summary(fit1, start = 1, recovery = TRUE, ps = p_vector[model@pnames], verbose = TRUE)
#                  A.紅.保守派 A.紅.自由派   A.綠.保守派   A.綠.自由派  A.藍.保守派
# True                1.2300      0.2500      1.8900      0.5500      1.6700
# 2.5% Estimate       1.0242      0.0606      1.6298      0.5560      1.5725
# 50% Estimate        1.3836      0.2586      2.1201      0.8195      2.0566
# 97.5% Estimate      1.6898      0.4568      2.5045      1.0700      2.5768
# Median-True         0.1536      0.0086      0.2301      0.2695      0.3866
#                A.藍.自由派    A.黃.保守派  A.黃.自由派      B mean_v.false
# True                0.4500      1.4500      0.3500 1.2500       1.1500
# 2.5% Estimate       0.2998      1.4688      0.0920 1.0097       0.9491
# 50% Estimate        0.5585      1.9324      0.3049 1.3333       1.3917
# 97.5% Estimate      0.8063      2.3755      0.4881 1.7094       1.7631
# Median-True         0.1085      0.4824     -0.0451 0.0833       0.2417
#                mean_v.true sd_v.true      t0
# True                2.8000    0.8000  0.1000
# 2.5% Estimate       2.6221    0.7126  0.0285
# 50% Estimate        3.0387    0.8100  0.0940
# 97.5% Estimate      3.4017    0.9124  0.1545
# Median-True         0.2387    0.0100 -0.0060


# expected_values <- "tests/testthat/Group0/data/expected_model4.rda"
expected_values <- "data/expected_model4.rda"
# samples <- samples0
# est <- est0
# save(samples, est, file = expected_values)

load(expected_values)
testthat::expect_true(all(samples0@theta[, , 1] == samples@theta[, , 1]))
testthat::expect_true(all(est0 == est))
