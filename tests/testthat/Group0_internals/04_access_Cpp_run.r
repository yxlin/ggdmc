#!/usr/bin/env Rscript
# q(save = "no")
cat("\n\n-------------------- Access Cpp run --------------------")
rm(list = ls())

pkg <- c("lbaModel", "ggdmcPrior", "ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)

wkdir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/Group1/"
cat("\nWorking directory: ", getwd(), "\n")
# helper_path <- paste0(wkdir, "helpers.r")
# save_path <- paste0(wkdir, "data/lba_data0.rda")
# source(helper_path)


model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    accumulators = c("r1", "r2"),
    type = "lba"
)

pop_mean <- c(A = .4, B = .5, mean_v.false = 0.15, mean_v.true = 2.5, t0 = 0.3)
pop_scale <- c(A = .1, B = .1, mean_v.false = .2, mean_v.true = .2, t0 = 0.05)
pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = c(0, 0, 0, 0, 0),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

# ---------------------------------------
sub_model <- setLBA(model)
pop_model <- setLBA(model, population_distribution = pop_dist)

p_vector <- c(A = .75, B = 1.25, mean_v.false = 1.5, mean_v.true = 2.5, t0 = .15)
dat <- simulate(sub_model, nsim = 256, parameter_vector = p_vector, n_subject = 1)
hdat <- simulate(pop_model, nsim = 128, n_subject = 32)

sub_dmis <- ggdmcModel::BuildDMI(dat, model)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)

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
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames, thin = 2)
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]], seed = 846671, verbose = FALSE)


de_input <- ggdmc::setDEInput(
    sub_migration_prob = .05,
    gamma_precursor = 2.38, rp = 0.001,
    is_pblocked = FALSE,
    nparameter = as.integer(model@npar),
    nchain = as.integer(model@npar * 3),
    sub_debug = FALSE
)

config_list <- ggdmc::set_configs(
    prior = sub_priors, theta_input = sub_theta_input,
    de_input = de_input,
    ncore = 1, seed = 123
)

fit0 <- ggdmc::run_subject(
    config_r = config_list[[1]], dmi = sub_dmis[[1]],
    samples = sub_samples
)

fit1 <- ggdmc::run_subject(
    config_r = config_list[[1]], dmi = sub_dmis[[1]],
    samples = fit0
)

fit <- fit1
hat <- ggdmc::gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)

pdf("Rplot.pdf")
p0 <- ggdmc:::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
dev.off()
