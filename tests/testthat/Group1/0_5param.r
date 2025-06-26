# q(save = "no")

cat("\n\n-------------------- 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data0.rda"
load(fn)

# New procedure ------------------------------------------
# theta_input <- setThetaInput(nmc = nmc, nchain = nchain, pnames = model@pnames, thin = 1)
de_input <- ggdmc::setDEInput(
    sub_migration_prob = 0.00,
    nparameter = as.integer(sub_theta_input@nparameter),
    nchain = as.integer(sub_theta_input@nchain)
)
configs <- ggdmc::set_configs(prior = sub_priors, theta_input = sub_theta_input, de_input = de_input)

theta_input <- ggdmc::setThetaInput(nmc = 2, nchain = 3, pnames = model@pnames, thin = 1)
de_input <- ggdmc::setDEInput(
    sub_migration_prob = 0.00,
    nparameter = as.integer(theta_input@nparameter),
    nchain = as.integer(theta_input@nchain)
)
configs <- ggdmc::set_configs(prior = sub_priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]

# pop_priors <- ggdmcPrior::set_priors(p_prior = model_likelihood, l_prior = l_prior, s_prior = s_prior)
# pop_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = pop_priors@pnames)

pop_samples <- ggdmc::initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671, verbose = FALSE)


# samples <- ggdmc::initialise_phi(theta_input, pop_priors, sub_dmis[[1]],
#     seed = 123
# )

# hyper_model <- ggdmcModel::BuildModel(
#     p_map = list(A = "1", B = "1", mean_v = "M", sd_v = "1", st0 = "1", t0 = "1"),
#     match_map = list(M = list(s1 = "r1", s2 = "r2")),
#     factors = list(S = c("s1", "s2")),
#     constants = c(sd_v = 1, st0 = 0),
#     accumulators = c("r1", "r2"),
#     type = "hyper"
# )
# hyper_dmi


# pop_samples$phi
fit0 <- run_hyper(cfg, hyper_dmi, pop_samples$phi)
head(run_hyper)
# fit0 <- run_subject(cfg, sub_dmis[[1]], samples, debug = T)


# de_input <- ggdmcDE::setDEInput(
#     sub_migration_prob = 0.00,
#     nparameter = as.integer(sub_theta_input@nparameter),
#     nchain = as.integer(sub_theta_input@nchain)
# )

# configs <- ggdmcDE::set_configs(prior = sub_priors, theta_input = sub_theta_input, de_input = de_input)
# cfg <- configs[[1]]
# fit1 <- run_subject(cfg, sub_dmis[[1]], fit0, debug = FALSE)

# fit <- fit0
# start <- fit@nmc * 0.5
# options(digits = 2)
# est_theta <- ggdmcPhi::compare(fit, ps = p_vector)

#                    A      B mean_v.false mean_v.true     t0
# True            0.750  1.250         1.50        2.50 0.1500
# 5% Estimate     0.049  0.625         0.49        1.70 0.0358
# 50% Estimate    0.460  1.186         1.25        2.29 0.1568
# 97.5% Estimate  1.269  1.996         2.13        3.15 0.2790
# Median-True    -0.290 -0.064        -0.25       -0.21 0.0068

# gdmc::gelman(fit)
