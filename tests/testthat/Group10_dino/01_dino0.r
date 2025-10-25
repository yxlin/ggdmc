#!/usr/bin/env Rscript
# Baseline CDM Model: DINO Rule
# Subject-level parameter estimation with N=3000
cat("\n\n--------- DINO Baseline Model (Subject-level) ----------\n")
rm(list = ls())
# q(save = "no")

pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group10_dino/data")
fig_dir <- file.path(home_dir, "Group10_dino/figs")
save_path <- paste0(home_dir, "Group10_dino/data/dino0.rda")

# -------------------- Q-Matrix Setup --------------------
# 3 items, 2 skills: Item 1 (Algebra only), Item 2 (Geometry only), Item 3 (both)
K <- 3
Q <- matrix(c(
    1, 1, 1,
    0, 0, 0,
    0, 1, 0,
    1, 0, 0,
    0, 0, 1,
    1, 0, 1,
    0, 1, 1,
    1, 1, 0
), ncol = K, byrow = TRUE)
J <- nrow(Q)
colnames(Q) <- c("Algebra", "Geometry", "Science")
rownames(Q) <- paste0("Item", seq_len(J))

# Calculate skill probabilities (no correlation, sigma = 0)
sigma <- .2
means <- c(.1, .15, .45)
Sigma <- cdModel::build_correlation_matrix(K, sigma)
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)
print(skill_probs)

# -------------------- Model Definition --------------------
simulation_model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        guess4 = "1", guess5 = "1", guess6 = "1",
        guess7 = "1", guess8 = "1",
        mean1 = "1", mean2 = "1", mean3 = "1", simga = "1",
        slip1 = "1", slip2 = "1", slip3 = "1",
        slip4 = "1", slip5 = "1", slip6 = "1",
        slip7 = "1", slip8 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

# Set DINA rule without multivariate normal
use_mvn <- TRUE
sub_model <- setCDM(simulation_model,
    q_matrix = simulation_model@cdm_info$q_matrix,
    rule = "DINO",
    use_mvn = use_mvn
)

# -------------------- Data Simulation --------------------
# True parameter values
sim_p_vector <- c(
    guess1 = .1, guess2 = .01, guess3 = .03,
    guess4 = .03, guess5 = .03, guess6 = .04,
    guess7 = .04, guess8 = .05,
    mean1 = means[1], mean2 = means[2], mean3 = means[3],
    sigma = sigma,
    slip1 = .05, slip2 = .02, slip3 = .02,
    slip4 = .02, slip5 = .02, slip6 = .03,
    slip7 = .03, slip8 = .05
)

true_p_vector <- c(
    guess1 = .1, guess2 = .01, guess3 = .03,
    guess4 = .03, guess5 = .03, guess6 = .04,
    guess7 = .04, guess8 = .05,
    pi_000 = skill_probs$probability[1],
    pi_100 = skill_probs$probability[2],
    pi_010 = skill_probs$probability[3],
    pi_110 = skill_probs$probability[4],
    pi_001 = skill_probs$probability[5],
    pi_101 = skill_probs$probability[6],
    pi_011 = skill_probs$probability[7],
    slip1 = .05, slip2 = .02, slip3 = .02,
    slip4 = .02, slip5 = .02, slip6 = .03,
    slip7 = .03, slip8 = .05
)

N <- 3000
dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = sim_p_vector,
    nschool = 1,
    debug = FALSE,
    seed = 123
)


pop_mean <- c(
    guess1 = .1, guess2 = .2, guess3 = .3,
    pi1 = 0.1, pi2 = 0.3, pi3 = 0.5,
    slip1 = .02,
    slip2 = .05, slip3 = .07
)

pop_scale <- c(
    guess1 = .01, guess2 = .01, guess3 = .01,
    pi1 = .01, pi2 = .03, pi3 = .05,
    slip1 = .02, slip2 = .05, slip3 = .07
)


pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINO"
)
pop_model <- setCDM(model,
    population_distribution = pop_dist,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINO"
)

p_vector <- c(
    guess1 = .1, guess2 = .2, guess3 = .3,
    pi1 = 0.1, pi2 = 0.3, pi3 = 0.5,
    slip1 = .02,
    slip2 = .05, slip3 = .07
)

N_total_student <- 1000
dat <- simulate(sub_model,
    nsim = N_total_student, parameter_vector = p_vector,
    nschool = 1,
    seed = 123,
    debug = FALSE
)

hdat <- simulate(pop_model,
    nsim = N_total_student,
    nschool = 32,
    seed = 123
)
ps <- attr(hdat, "parameters")

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = sub_model@q_matrix,
    profile_probability = sub_model@profile_probability,
    rule = "DINO"
)
pop_dmis <- BuildDMI(hdat$responses, model,
    q_matrix = sub_model@q_matrix,
    profile_probability = sub_model@profile_probability,
    rule = "DINO"
)


true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)

# Fit subject/school theta starting samples -------------------------
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(1.1, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

sub_priors <- set_priors(p_prior = p_prior)

nmc <- 500
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames)
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]],
    verbose = T
)

# Generate pop samples -------------------------
p0 <- runif(model@npar)
names(p0) <- model@pnames
model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    # Setting population model likelihod to 10 will result in
    # very difficult to convert; TODO: test if increase the sample size
    # will result this problem
    # p1 = rep(10, model@npar),

    # We know that the probability is within 0 and 1, bcz this
    # is the definition of the guess and slip parameters. That is,
    # they are probability values.
    p1 = rep(1.1, model@npar),
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

pop_priors <- ggdmcPrior::set_priors(
    p_prior = model_likelihood, l_prior = l_prior,
    s_prior = s_prior
)
pop_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = pop_priors@pnames)

pop_samples <- ggdmc::initialise_phi(pop_theta_input, pop_priors, pop_dmis,
    seed = 846671, verbose = FALSE
)

save(model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis,
    sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input, Q,
    file = save_path
)
