#!/usr/bin/env Rscript
cat("\n\n-------------------- Generate model 0 --------------------")
# q(save = "no")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
save_path <- paste0(home_dir, "Group9_gen_cdm/data/dina11.rda")

Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1
), ncol = 2, byrow = TRUE)

colnames(Q) <- c("Algebra", "Geometry")
rownames(Q) <- c("Item1", "Item2", "Item3")

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        pi1 = "1", pi2 = "1", pi3 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)


# pop_mean <- c(
#     guess1 = .1, guess2 = .2, guess3 = .3,
#     pi1 = 0.1, pi2 = 0.3, pi3 = 0.5,
#     slip1 = .02, slip2 = .05, slip3 = .07
# )
pop_mean <- c(
    guess1 = .1, guess2 = .2, guess3 = .3,
    pi1 = 2, pi2 = 2, pi3 = 2,
    slip1 = .02, slip2 = .05, slip3 = .07
)

# pop_scale <- c(
#     guess1 = .01, guess2 = .01, guess3 = .01,
#     pi1 = .01, pi2 = .03, pi3 = .05,
#     slip1 = .02, slip2 = .05, slip3 = .07
# )

pop_scale <- c(
    guess1 = .01, guess2 = .01, guess3 = .01,
    pi1 = NA, pi2 = NA, pi3 = NA,
    slip1 = .02, slip2 = .05, slip3 = .07
)

pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = c(
        "tnorm", "tnorm", "tnorm", "dirichlet", "dirichlet", "dirichlet",
        "tnorm", "tnorm", "tnorm"
    ),
    log_p = rep(F, model@npar)
)
print_prior(pop_dist)
plot_prior(pop_dist)



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
# sub_model <- setCDM(model)
# pop_model <- setCDM(model, population_distribution = pop_dist)
# p_vector <- c(
#     guess1 = .1, guess2 = .2, guess3 = .3,
#     pi1 = 2, pi2 = 2, pi3 = 2,
#     slip1 = .02,
#     slip2 = .05, slip3 = .07
# )
set.seed(123)
p_vector <- rprior(pop_dist)[, 1]
cat("p_vector:\n")
cat("Sum of first 3 components:", sum(p_vector[4:6]), "\n")
cat("Implied 4th component (1 - sum):", 1 - sum(p_vector[4:6]), "\n")
names(p_vector) <- names(pop_mean)

cat("True p_vector\n")
options(digits = 2)
print(p_vector)

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

# true_vector

# Fit subject/school theta starting samples -------------------------
p0 <- c(0, 0, 0, 10, 10, 10, 0, 0, 0)
names(p0) <- model@pnames

# pop_dist <- ggdmcPrior::BuildPrior(
#     p0 = pop_mean,
#     p1 = pop_scale,
#     lower = rep(0, model@npar),
#     upper = rep(NA, model@npar),
#     dists = c(
#         "tnorm", "tnorm", "tnorm", "dirichlet", "dirichlet", "dirichlet",
#         "tnorm", "tnorm", "tnorm"
#     ),
#     log_p = rep(F, model@npar)
# )

p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = c(1, 1, 1, NA, NA, NA, 1, 1, 1),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = c(
        "unif", "unif", "unif", "dirichlet", "dirichlet", "dirichlet",
        "unif", "unif", "unif"
    ),
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
    p1 = c(1, 1, 1, NA, NA, NA, 1, 1, 1),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = c(
        "unif", "unif", "unif", "dirichlet", "dirichlet", "dirichlet",
        "unif", "unif", "unif"
    ),
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

cat(paste0("save file to ", save_path), "\n")
