# q(save = "no")
cat("\n\n-------------------- Generate DDM a model --------------------")
rm(list = ls())
pkg <- c("ddModel", "ggdmcPrior")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/data/"
save_path <- paste0(wkdir, "ddm_data3.rda")


hyper_model <- ggdmcModel::BuildModel(
    p_map = list(
        a = c("S", "政黨傾向"), v = "1", z = "1", d = "1", sz = "1", sv = "1",
        t0 = "1", st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2"), 政黨傾向 = c("自由派", "保守派")),
    constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
    accumulators = c("r1", "r2"),
    type = "hyper",
    verbose = FALSE
)

model <- ggdmcModel::BuildModel(
    p_map = list(
        a = c("S", "政黨傾向"), v = "1", z = "1", d = "1", sz = "1", sv = "1",
        t0 = "1", st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2"), 政黨傾向 = c("自由派", "保守派")),
    constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
    accumulators = c("r1", "r2"),
    type = "fastdm"
)

pop_mean <- c(
    a.s1.保守派 = 1,
    a.s1.自由派 = 2.5,
    a.s2.保守派 = 1.5,
    a.s2.自由派 = 3.5,
    sz = 0.25, t0 = 0.15, v = 2.2, z = 0.38
)
pop_scale <- c(
    a.s1.保守派 = 0.05,
    a.s1.自由派 = 0.06,
    a.s2.保守派 = 0.05,
    a.s2.自由派 = 0.08,
    sz = 0.02, t0 = 0.02, v = 0.5, z = 0.03
)

pop_dist <- ggdmcPrior::BuildPrior(
    p0    = pop_mean,
    p1    = pop_scale,
    lower = c(0, 0, 0, 0, 0, 0, -10, 0),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)


# ---------------------------------------
sub_model <- setDDM(model)
pop_model <- setDDM(model, population_distribution = pop_dist)

p_vector <- c(
    a.s1.保守派 = 1,
    a.s1.自由派 = 2.5,
    a.s2.保守派 = 1.5,
    a.s2.自由派 = 3.5,
    sz = 0.25, t0 = 0.15, v = 2.2, z = 0.38
)


dat <- simulate(sub_model, nsim = 64, parameter_vector = p_vector, n_subject = 1, debug = FALSE)
hdat <- simulate(pop_model,
    nsim = 256, n_subject = 32
)

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
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]], seed = 846671, verbose = F)

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

pop_samples <- ggdmc::initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671, verbose = FALSE)

save(hyper_model, model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi,
    sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input,
    file = save_path
)
