#!/usr/bin/env Rscript
# Generate CDM Data with 3 Skills - Simple Test Case
q(save = "no")
cat("\n\n--------- Generate CDM Data (DINA + 3 Skills) ----------\n")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
data_dir <- paste0(home_dir, "Group9_gen_cdm/data/")

# Simple Q-matrix: 3 skills, 5 items
# - Items 1-3: Single skills (identity matrix)
# - Item 4: Skills 1+2
# - Item 5: All skills
Q <- matrix(c(
    1, 0, 0, # Item1: Skill1 only
    0, 1, 0, # Item2: Skill2 only
    0, 0, 1, # Item3: Skill3 only
    1, 1, 0, # Item4: Skills 1+2
    1, 1, 1 # Item5: All skills
), ncol = 3, byrow = TRUE)

colnames(Q) <- c("Skill1", "Skill2", "Skill3")
rownames(Q) <- paste0("Item", 1:5)
cat("\nQ-matrix:\n")
print(Q)

# === STEP 1: Generate data using MVN profile probabilities ===
model_gen <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1", guess5 = "1",
        mean1 = "1", mean2 = "1", mean3 = "1", sigma = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

use_mvn <- TRUE
sub_model_gen <- setCDM(model_gen,
    q_matrix = model_gen@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = use_mvn
)

# Calculate profile probabilities from MVN parameters
K <- ncol(Q)
sigma <- 0 # Independent skills
means <- rep(0, K) # Equal mastery probability
Sigma <- cdModel::build_correlation_matrix(K, sigma)
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)

cat("\nProfile probabilities (from MVN with mean=0, sigma=0):\n")
print(skill_probs)

# True parameter vector for data generation
sim_p_vector <- c(
    guess1 = .10, guess2 = .15, guess3 = .20, guess4 = .25, guess5 = .30,
    mean1 = means[1], mean2 = means[2], mean3 = means[3], sigma = sigma,
    slip1 = .08, slip2 = .10, slip3 = .12, slip4 = .15, slip5 = .18
)

cat("\nTrue parameters:\n")
print(sim_p_vector)

# Generate data
N <- 100
set.seed(123)
dat <- simulate(sub_model_gen,
    nsim = N,
    parameter_vector = sim_p_vector,
    nschool = 1,
    debug = FALSE,
    seed = 123
)

cat("\nGenerated", N, "subjects\n")
cat("Response summary:\n")
print(summary(dat$responses))


# === STEP 2: Fit model with discrete profile probabilities ===
# 3 skills → 2^3 = 8 profiles → estimate 7 profile probs
model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1", guess5 = "1",
        pi_000 = "1", pi_100 = "1", pi_010 = "1", pi_110 = "1",
        pi_001 = "1", pi_101 = "1", pi_011 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

cat("\nFitting model - Parameter names:\n")
print(model@pnames)
cat("Total parameters:", model@npar, "\n")

# Setup prior with Dirichlet for profile probabilities
p0 <- rep(0, model@npar)
names(p0) <- model@pnames

# Identify profile probability parameters
pi_params <- grepl("^pi_", model@pnames)
cat("\nProfile probability parameters:\n")
print(model@pnames[pi_params])

# Set p0: alpha MUST be > 0 for Dirichlet
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p0[pi_params] <- 1.0 # ✓ FIX: Set alpha = 1 (or 2, or any positive value)

# Set p1: NA for Dirichlet
p1 <- rep(1.0, model@npar)
names(p1) <- model@pnames
p1[pi_params] <- NA # ✓ This part was correct

# Build prior
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = p1,
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dists = ifelse(pi_params, "dirichlet", "unif"), # Note: 'dists' not 'dist'!
    log_p = rep(TRUE, model@npar)
)

# IMPORTANT: For Dirichlet, p0 = alpha (concentration parameter) must be > 0
# Common choices:
#   alpha = 1.0: Uniform over simplex (flat prior)
#   alpha = 2.0: Symmetric, slightly concentrated
#   alpha = 10.0: Strongly concentrated around equal probabilities
# Using alpha = 1.0 for uniform prior
# p0 <- c(
#     0, 0, 0, 0, 0,              # guess parameters (not Dirichlet)
#     1.0, 1.0, 1.0, 1.0,         # pi_000, pi_100, pi_010, pi_110 (Dirichlet alpha)
#     1.0, 1.0, 1.0,              # pi_001, pi_101, pi_011 (Dirichlet alpha)
#     0, 0, 0, 0, 0               # slip parameters (not Dirichlet)
# )
# p1 <- c(
#     1, 1, 1, 1, 1,
#     NA, NA, NA, NA, NA,
#     NA, NA,
#     1, 1, 1, 1, 1
# )
# # p0 <- rep(0, model@npar)
# names(p0) <- model@pnames
# names(p1) <- model@pnames

# p_prior <- ggdmcPrior::BuildPrior(
#     p0 = p0,
#     p1 = rep(1.0, model@npar),
#     lower = rep(NA, model@npar),
#     upper = rep(NA, model@npar),
#     dist = rep("unif", model@npar),
#     log_p = rep(TRUE, model@npar)
# )
# sub_priors <- set_priors(p_prior = p_prior)
# print_prior(p_prior)

cat("\nPrior distributions:\n")
print(data.frame(
    param = model@pnames,
    dist = ifelse(pi_params, "dirichlet", "unif")
))

sub_priors <- set_priors(p_prior = p_prior)

# Build DMI
sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = FALSE
)

cat("\nStarting MCMC sampling...\n")
head(StartSampling_subject)
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    nchain = 30, rp = 0.5,
    sub_migration_prob = 0.15, thin = 4, is_pblocked = TRUE,
    seed = 9032
)
fits0[[1]]@nchain

cat("\nRestarting with migration_prob = 0...\n")
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.05, thin = 4, is_pblocked = TRUE
)

# 3 * model@npar
# Plot chains
cat("\nPlotting chains from start...\n")
p0 <- ggdmc::plot(fits0[[1]], start = 1, hide_legend = FALSE)
p0 <- ggdmc::plot(fits0[[2]], start = 1)
p0 <- ggdmc::plot(fits0[[3]], start = 1)

cat("\nPlotting chains after burn-in...\n")
p0 <- ggdmc::plot(fits0[[1]], start = 200)
p0 <- ggdmc::plot(fits0[[2]], start = 200)
p0 <- ggdmc::plot(fits0[[3]], start = 200)

cat("\nPlotting chains after burn-in...\n")
p0 <- ggdmc::plot(fits1[[1]], start = 200)
p0 <- ggdmc::plot(fits1[[2]], start = 200)
p0 <- ggdmc::plot(fits1[[3]], start = 200, hide_legend = FALSE)

# Check convergence
fits <- fits1
fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("\n=== Convergence Diagnostics ===\n")
cat("Overall mpsrf =", hat$mpsrf, "\n")

ncore <- 3
for (i in seq_len(ncore)) {
    hat <- gelman(fits[[i]])
    cat("Chain", i, "mpsrf =", hat$mpsrf, "\n")
}

# Compare with true values
options(digits = 3)
true_vector <- c(
    guess1 = .10, guess2 = .15, guess3 = .20, guess4 = .25, guess5 = .30,
    pi_000 = skill_probs$probability[1],
    pi_100 = skill_probs$probability[2],
    pi_010 = skill_probs$probability[3],
    pi_110 = skill_probs$probability[4],
    pi_001 = skill_probs$probability[5],
    pi_101 = skill_probs$probability[6],
    pi_011 = skill_probs$probability[7],
    slip1 = .08, slip2 = .10, slip3 = .12, slip4 = .15, slip5 = .18
)

cat("\n=== Parameter Estimates vs True Values ===\n")
est_theta <- ggdmc::compare(fit, ps = true_vector)

cat("\n=== Summary Statistics ===\n")
est_summary <- ggdmc::summary(fit)
print(est_summary$quantiles)

cat("\n=== Done ===\n")
