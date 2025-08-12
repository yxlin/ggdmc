q(save = "no")
cat("\n\n-------------------- Generate DDM a, v , z, t0 model --------------------")
rm(list = ls())
pkg <- c("ddModel", "ggdmcPrior")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/data/"
save_path <- paste0(wkdir, "ddm_data8.rda")


hyper_model <- ggdmcModel::BuildModel(
    p_map = list(
        a = "COLOUR", v = c("NOISE", "S"), z = "P", d = "1", sz = "1", sv = "1",
        t0 = "MOTOR_SKILL", st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(left = "z_key", right = "x_key")),
    factors = list(
        S = c("left", "right"),
        NOISE = c("high", "moderate", "low"),
        COLOUR = c("red", "blue"),
        MOTOR_SKILL = c("excellent", "good", "poor"),
        P = c("20", "80")
    ),
    constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
    accumulators = c("z_key", "x_key"),
    type = "hyper",
    verbose = FALSE
)
model <- ggdmcModel::BuildModel(
    p_map = list(
        a = "COLOUR", v = c("NOISE", "S"), z = "P", d = "1", sz = "1", sv = "1",
        t0 = "MOTOR_SKILL", st0 = "1", s = "1", precision = "1"
    ),
    match_map = list(M = list(left = "z_key", right = "x_key")),
    factors = list(
        S = c("left", "right"),
        NOISE = c("high", "moderate", "low"),
        COLOUR = c("red", "blue"),
        MOTOR_SKILL = c("excellent", "good", "poor"),
        P = c("20", "80")
    ),
    constants = c(d = 0, s = 1, st0 = 0, sv = 0, precision = 3),
    accumulators = c("z_key", "x_key"),
    type = "fastdm"
)


pop_mean <- c(
    a.blue = 2.0, # Slightly lower boundary for blue stimuli
    a.red = 2.2, # Slightly higher boundary for red stimuli
    sz = 0.25, # Start point variability (same as before)
    t0.excellent = 0.12, # Fastest non-decision time
    t0.good = 0.15, # Intermediate non-decision time
    t0.poor = 0.18, # Slowest non-decision time
    v.left.high = 1.0, # High drift rate for left choices
    v.left.low = 2.5, # Low drift rate for left choices
    v.left.moderate = 1.8, # Moderate drift rate for left choices
    v.right.high = 1.1, # High drift rate for right choices
    v.right.low = 2.4, # Low drift rate for right choices
    v.right.moderate = 1.7, # Moderate drift rate for right choices
    z.20 = 0.2, # Bias toward lower boundary (20%)
    z.80 = 0.8 # Bias toward upper boundary (80%)
)

pop_scale <- c(
    a.blue = 0.1,
    a.red = 0.12,
    sz = 0.02,
    t0.excellent = 0.015,
    t0.good = 0.02,
    t0.poor = 0.025,
    v.left.high = 0.5,
    v.left.low = 0.3,
    v.left.moderate = 0.4,
    v.right.high = 0.5,
    v.right.low = 0.3,
    v.right.moderate = 0.4,
    z.20 = 0.05,
    z.80 = 0.05
)

pop_dist <- ggdmcPrior::BuildPrior(
    p0    = pop_mean,
    p1    = pop_scale,
    lower = c(0, 0, 0, 0, 0, 0, -10, -10, -10, -10, -10, -10, 0, 0),
    upper = c(NA, NA, NA, 1, 1, 1, rep(NA, 8)), # t0.* bounded at 1
    dists = rep("tnorm", length(pop_mean)),
    log_p = rep(FALSE, length(pop_mean))
)

# ---------------------------------------
sub_model <- setDDM(model)
pop_model <- setDDM(model, population_distribution = pop_dist)

p_vector <- c(
    a.blue = 2.0, # Slightly lower boundary for blue stimuli
    a.red = 2.2, # Slightly higher boundary for red stimuli
    sz = 0.25, # Start point variability (same as before)
    t0.excellent = 0.12, # Fastest non-decision time
    t0.good = 0.15, # Intermediate non-decision time
    t0.poor = 0.18, # Slowest non-decision time
    v.left.high = 1.0, # High drift rate for left choices
    v.left.low = 2.5, # Low drift rate for left choices
    v.left.moderate = 1.8, # Moderate drift rate for left choices
    v.right.high = 1.1, # High drift rate for right choices
    v.right.low = 2.4, # Low drift rate for right choices
    v.right.moderate = 1.7, # Moderate drift rate for right choices
    z.20 = 0.2, # Bias toward lower boundary (20%)
    z.80 = 0.8 # Bias toward upper boundary (80%)
)

dat <- simulate(sub_model,
    nsim = 1728,
    parameter_vector = p_vector, n_subject = 1
)

hdat <- simulate(pop_model,
    nsim = 1728, n_subject = 32
)


sub_dmis <- ggdmcModel::BuildDMI(dat, model)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)
hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)

options(digits = 3)
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
