# q(save = "no")
cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
save_path <- paste0(home_dir, "Group9_gen_cdm/data/dina9.rda")

Q <- as.matrix(CDM::data.ecpe$q.matrix[1:5, ])
rownames(Q) <- c("Item1", "Item2", "Item3", "Item4", "Item5")
str(Q)

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        guess4 = "1", guess5 = "1",
        pi1 = "1", pi2 = "1", pi3 = "1", pi4 = "1", pi5 = "1",
        pi6 = "1", pi7 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1",
        slip4 = "1", slip5 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

pop_mean <- c(
    guess1 = 0.71, guess2 = 0.72, guess3 = 0.44,
    guess4 = 0.48, guess5 = 0.76,
    pi1 = 0.31, pi2 = 0.01, pi3 = 0.04, pi4 = 0.05,
    pi5 = 0.01, pi6 = 0.03, pi7 = 0.10,
    slip1 = 0.085, slip2 = 0.11, slip3 = 0.27, slip4 = 0.16,
    slip5 = 0.04
)

pop_scale <- c(
    guess1 = 0.01, guess2 = 0.026, guess3 = 0.01,
    guess4 = 0.02, guess5 = 0.01,
    pi1 = 0.017, pi2 = 0.009, pi3 = 0.012, pi4 = 0.010,
    pi5 = 0.008, pi6 = 0.007, pi7 = 0.011,
    slip1 = 0.009, slip2 = 0.008, slip3 = 0.01, slip4 = 0.01,
    slip5 = 0.005
)



pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)


# 1. The probability of occurrence in the population, written as P(alpha_l)
# 2. profile probability, gamma (𝛾) = [P(alpha_1), P(alpha_2), ..., P(alpha_{2^K})].
sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINA"
)
pop_model <- setCDM(model,
    population_distribution = pop_dist,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINA"
)


p_vector <- c(
    guess1 = 0.71, guess2 = 0.72, guess3 = 0.44,
    guess4 = 0.48, guess5 = 0.76,
    pi1 = 0.31, pi2 = 0.01, pi3 = 0.04, pi4 = 0.05,
    pi5 = 0.01, pi6 = 0.03, pi7 = 0.10,
    slip1 = 0.085, slip2 = 0.11, slip3 = 0.27, slip4 = 0.16,
    slip5 = 0.04
)



N_total_student <- 2000
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
# str(hdat)
# head(hdat$response)
# table(hdat$response$s)
# table(dat$response$s)

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = sub_model@q_matrix,
    profile_probability = sub_model@profile_probability, rule = "DINA"
)



pop_dmis <- BuildDMI(hdat$responses, model,
    q_matrix = sub_model@q_matrix,
    profile_probability = sub_model@profile_probability, rule = "DINA"
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
    verbose = F
)


# Generate pop samples -------------------------
p0 <- runif(model@npar)
names(p0) <- model@pnames
model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    # 1. Setting p1 to 10 in the hyper lileihood function, rendering
    # the model hard to converge;
    # 2. TODO: test if increase the sample size will resolve this problem
    # p1 = rep(10, model@npar),
    # 3. We know that the probability is within 0 and 1, bcz this
    # is what CDM defines it item parameters. They are probabilities.
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
save_path

save(model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis,
    sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input, Q,
    file = save_path
)
cat(paste0("save file to ", save_path), "\n")
