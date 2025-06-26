# q(save = "no")
cat("\n\n-------------------- Generate model 1 --------------------")
rm(list = ls())
pkg <- c("lbaModel", "ggdmcPrior", "ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
save_path <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data1.rda"
helper_path <- "~/Documents/ggdmc/tests/testthat/Group1/data/helpers.r"
source(helper_path)

hyper_model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "M", st0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "hyper", verbose = FALSE
)


model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "M", st0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "lba"
)
pop_mean <- c(A = .4, B = .5, mean_v.false = .15, mean_v.true = 2.5, sd_v.true = .1, t0 = .3)
pop_scale <- c(A = .1, B = .1, mean_v.false = .2, mean_v.true = .2, sd_v.true = .1, t0 = .05)
pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(FALSE, model@npar)
)


# ---------------------------------------
sub_model <- setLBA(model)
pop_model <- setLBA(model, population_distribution = pop_dist)

p_vector <- c(
    A = .75, B = 1.25, mean_v.false = 1.5, mean_v.true = 2.5,
    sd_v.true = 0.1, t0 = 0.15
)
dat <- simulate(sub_model, nsim = 256, parameter_vector = p_vector, n_subject = 1)
hdat <- simulate(pop_model, nsim = 256, n_subject = 32)

sub_dmis <- ggdmcModel::BuildDMI(dat, model)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)
hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)

res <- simple_get_accuracy(dat)
res <- simple_get_accuracy(hdat)
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

nmc <- 1000
nchain <- model@npar * 3
thin <- 2

sub_theta_input <- setThetaInput(nmc = nmc, nchain = nchain, pnames = model@pnames, thin = thin, report_length = 200)

sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]], seed = 846671, verbose = FALSE)


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

pop_nchain <- pop_priors@nparameter * 3

pop_theta_input <- ggdmc::setThetaInput(
    nmc = nmc, nchain = pop_nchain,
    pnames = pop_priors@pnames, thin = thin, report_length = 200
)

pop_samples <- ggdmcPhi::initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671, verbose = FALSE)

save(hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi, sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input,
    file = save_path
)
