q(save = "no")
cat("\n\n-------------------- Generate RRUM model 0 --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
save_path <- paste0(home_dir, "Group9_gen_cdm/data/rrum0.rda")
# At the moment, the profile_probability is fixed via the setCDM, and BuildDMI functions.
# The remaining questions are:
# Q1. Should we put it in BuildModel to estimate?
# Q2. If not where should we put it?
# Q3. How should we adjust pop_dist if we want to estimate profile_probability, too?
# Q4. Perhaps as a second step, we should not allow it to vary with other factors.

Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("Algebra", "Geometry")
rownames(Q) <- c("Item1", "Item2", "Item3")

model <- BuildModel(
    p_map = list(
        pi1 = "1", pi2 = "1", pi3 = "1",
        pi_item1 = "1", pi_item2 = "1", pi_item3 = "1",
        r11 = "1", r22 = "1",
        r31 = "1", r32 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL, # Must enter NULL
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)
pop_mean <- c(
    pi1 = .1, pi2 = .2, pi3 = .3,
    pi_item1 = .1, pi_item2 = .2, pi_item3 = .3,
    r11 = .1,
    r22 = .4,
    r31 = .3, r32 = .6
)
pop_scale <- c(
    pi1 = .01, pi2 = .02, pi3 = .03,
    pi_item1 = .01, pi_item2 = .02, pi_item3 = .03,
    r11 = .01,
    r22 = .04,
    r31 = .03, r32 = .06
)


pop_dist <- BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "RRUM"
)


pop_model <- setCDM(model,
    population_distribution = pop_dist,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "RRUM"
)

p_vector <- c(
    pi1 = .1, pi2 = .2, pi3 = .3,
    pi_item1 = .1, pi_item2 = .2, pi_item3 = .3,
    r11 = .1,
    r22 = .4,
    r31 = .3, r32 = .6
)


p0 <- rep(0, model@npar)
names(p0) <- model@pnames

p_prior <- BuildPrior(
    p0 = p0,
    p1 = rep(1, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

sub_priors <- set_priors(p_prior = p_prior)


N_total_student <- 1000
dat <- simulate(sub_model,
    nsim = N_total_student, parameter_vector = p_vector,
    nschool = 1,
    seed = 123,
    debug = F
)

hdat <- simulate(pop_model,
    nsim = N_total_student,
    nschool = 32,
    seed = 123
)


options(digits = 3)
ps <- attr(hdat, "parameters")



sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = sub_model@q_matrix,
    profile_probability = sub_model@profile_probability, rule = "RRUM"
)
pop_dmis <- BuildDMI(hdat$responses, model,
    q_matrix = pop_model@q_matrix,
    profile_probability = pop_model@profile_probability, rule = "RRUM"
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
    p1 = rep(1.05, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

sub_priors <- set_priors(p_prior = p_prior)


nmc <- 500
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames)
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]],
    verbose = F
)


# Generate pop samples -------------------------
p0 <- runif(model@npar)
names(p0) <- model@pnames
model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    # 1. Setting p1 = 10 tends to hinder convergence.
    # 2. TODO: Test whether increasing the sample size mitigates this.
    # 3. Guess and slip are probabilities; constrain them to [0, 1].
    # p1 <- rep(10, model@npar)
    p1 = rep(1.05, model@npar),
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
    # 1) Weakly informative hyperpriors on population location/scale are acceptable.
    # 2) In a hierarchical model, p_prior represents the model likelihood
    #    at the population level. It denotes p(theta | phi). This presumes
    #    the optimiser knows the data (which, in this case, is theta), given
    #    the phi parameter (for which the optimiser makes a new guess at each
    #    iteration).
    # 3) The subject-level likelihood, denoted p(y | theta), is given by the
    #    DDM, LBA, CDM, etc., likelihood functions. This presumes that
    #    the optimiser knows the data (which, in this case, is y), given
    #    the theta parameter (for which the optimiser makes a new guess at each
    #    iteration).
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
