# q(save = "no")
cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
save_path <- paste0(home_dir, "Group9_gen_cdm/data/dino1.rda")

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = NULL,
    type = "cdm",
    verbose = TRUE
)

hyper_model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = NULL,
    type = "hyper",
    verbose = FALSE
)


pop_mean <- c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4,
    slip1 = .01, slip2 = .02, slip3 = .03, slip4 = .04
)
pop_scale <- c(
    guess1 = .01, guess2 = .01, guess3 = .01, guess4 = .01,
    slip1 = .05,
    slip2 = .03, slip3 = .02, slip4 = .01
)


pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1,
    1, 0
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")

n_item <- nrow(Q)
n_skill <- ncol(Q)
n_profile <- 2^(n_skill)
pi_uniform <- rep(1 / n_profile, n_profile)

# Investigate why using the default setting
sub_model <- setCDM(model, q_matrix = Q, prior_pi = pi_uniform, rule = "DINO")
pop_model <- setCDM(model,
    population_distribution = pop_dist,
    q_matrix = Q, prior_pi = pi_uniform, rule = "DINO"
)


p_vector <- c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4,
    slip1 = .01, slip2 = .02, slip3 = .03, slip4 = .04
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
    prior_pi = sub_model@prior_pi,
    rule = "DINO"
)
pop_dmis <- BuildDMI(hdat$responses, model,
    q_matrix = pop_model@q_matrix,
    prior_pi = pop_model@prior_pi,
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
    verbose = F
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
    pop_priors, pop_samples, pop_theta_input,
    file = save_path
)
