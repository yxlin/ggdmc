# q(save = "no")
cat("\n\n-------------------- Generate BV data --------------------")
rm(list = ls())
pkg <- c("lbaModel", "ggdmcPrior", "ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ggdmc/tests/testthat/Group1/data/"
helper_path <- paste0(wkdir, "helpers.r")
save_path <- paste0(wkdir, "lba_data6.rda")
source(helper_path)


hyper_model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = c("S", "COLOR"), t0 = "1", mean_v = c("NOISE", "M"), sd_v = "M", st0 = "1"),
    match_map = list(M = list(left = "z_key", right = "x_key")),
    factors = list(
        S = c("left", "right"),
        COLOR = c("red", "blue"),
        NOISE = c("high", "moderate", "low")
    ),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("z_key", "x_key"),
    type = "hyper", verbose = FALSE
)


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

res <- zx_get_accuracy(dat)
res <- zx_get_accuracy(hdat)
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
thin <- 2
sub_theta_input <- setThetaInput(nmc = nmc, pnames = model@pnames, thin = thin, report_length = 200)

sub_samples <- initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]], seed = 483128, verbose = FALSE)

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

pop_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = pop_priors@pnames, thin = thin, report_length = 200)

pop_samples <- initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671, verbose = FALSE)


save(hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi,
    sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input,
    file = save_path
)
