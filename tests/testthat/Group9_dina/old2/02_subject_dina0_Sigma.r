#!/usr/bin/env Rscript
# Baseline CDM Model: DINA Rule without MVN
# Subject-level parameter estimation with N=2000
cat("\n\n--------- CDM DINA Baseline Model (Subject-level) ----------\n")
rm(list = ls())
# q(save = "no")
# Load packages
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat/"
data_dir <- paste0(home_dir, "Group9_gen_cdm/data/")

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
colnames(Q) <- c("Algebra", "Geometry", "Science")
rownames(Q) <- paste0("Item", 1:nrow(Q))


# Calculate skill probabilities (no correlation, sigma = 0)
sigma <- 0
means <- rep(0, K)
Sigma <- cdModel::build_correlation_matrix(K, sigma)
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)
skill_probs


# -------------------- Model Definition --------------------
model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        guess4 = "1", guess5 = "1", guess6 = "1",
        guess7 = "1", guess8 = "1",
        pi_000 = "1", pi_100 = "1", pi_010 = "1",
        pi_110 = "1", pi_001 = "1", pi_101 = "1",
        pi_011 = "1",
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
use_mvn <- FALSE
sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = use_mvn
)

# -------------------- Data Simulation --------------------
# True parameter values
sim_p_vector <- c(
    guess1 = .3, guess2 = .001, guess3 = .1,
    guess4 = .1, guess5 = .1, guess6 = .2,
    guess7 = .22, guess8 = .24,
    pi_000 = skill_probs$probability[1],
    pi_100 = skill_probs$probability[2],
    pi_010 = skill_probs$probability[3],
    pi_110 = skill_probs$probability[4],
    pi_001 = skill_probs$probability[5],
    pi_101 = skill_probs$probability[6],
    pi_011 = skill_probs$probability[7],
    slip1 = .1, slip2 = .001, slip3 = .03,
    slip4 = .01, slip5 = .02, slip6 = .03,
    slip7 = .03, slip8 = .03
)
# length(sim_p_vector)
N <- 2000
dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = sim_p_vector,
    nschool = 1,
    debug = FALSE,
    seed = 123
)

# -------------------- Prior Setup --------------------
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(1.0, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
sub_priors <- set_priors(p_prior = p_prior)
# plot_prior(p_prior)

# Build data-model-info objects
sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = use_mvn
)

# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration plus blocking at the
# subject level (exploration phase)
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 2, is_pblocked = TRUE,
    seed = 9032
)
# Processing time (in R): 783.57 secs.

# Stage 1: First restart without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
)

# Stage 2: Final sampling
fits2 <- ggdmc:::RestartSampling_subject(fits1,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
)

# -------------------- Diagnostics (Optional) --------------------
# Check Stage 0: Burn-in chains
p0 <- ggdmc::plot(fits0[[1]], start = 1)
p0 <- ggdmc::plot(fits0[[2]], start = 1)
p0 <- ggdmc::plot(fits0[[3]], start = 1)

p0 <- ggdmc::plot(fits0[[1]], start = 200)
p0 <- ggdmc::plot(fits0[[2]], start = 200)
p0 <- ggdmc::plot(fits0[[3]], start = 200)

# Check Stage 1: First restart
p0 <- ggdmc::plot(fits1[[1]])
p0 <- ggdmc::plot(fits1[[2]])
p0 <- ggdmc::plot(fits1[[3]])

# Check Stage 2: Final chains
p0 <- ggdmc::plot(fits2[[1]])
p0 <- ggdmc::plot(fits2[[2]])
p0 <- ggdmc::plot(fits2[[3]])

# -------------------- Convergence Check --------------------
fits <- fits2
fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("Overall mpsrf = ", hat$mpsrf, "\n")

# Check individual chains
ncore <- 3
for (i in seq_len(ncore)) {
    hat <- gelman(fits[[i]])
    cat("Chain", i, "mpsrf = ", hat$mpsrf, "\n")
}

# -------------------- Parameter Recovery --------------------
options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = sim_p_vector)

# # Example output (N=2000):
#               guess1 guess2 guess3 guess4 guess5 guess6 guess7 guess8  pi_000
# True            0.30  0.001  0.100   0.10   0.10   0.20   0.22   0.24  0.1250
# 5 Estimate      0.30  0.051  0.078   0.10   0.11   0.18   0.21   0.23  0.1045
# 50 Estimate     0.90  0.503  0.987   0.98   0.99   0.99   0.97   0.98  0.1188
# 97.5 Estimate   0.96  0.976  1.000   1.00   1.00   1.00   1.00   1.00  0.1376
# Median-True     0.60  0.502  0.887   0.88   0.89   0.79   0.75   0.74 -0.0062
#                 pi_001   pi_010  pi_011   pi_100  pi_101 pi_110 slip1   slip2
# True           0.12500  0.12500  0.1250  0.12500  0.1250  0.125 0.100 0.00100
# 5 Estimate     0.00028  0.00024  0.0034  0.00019  0.0075  0.004 0.051 0.00068
# 50 Estimate    0.00449  0.00350  0.0117  0.00274  0.0173  0.013 0.695 0.00182
# 97.5 Estimate  0.14758  0.14842  0.1334  0.13378  0.1427  0.129 0.723 0.00435
# Median-True   -0.12051 -0.12150 -0.1133 -0.12226 -0.1077 -0.112 0.595 0.00082
#                slip3  slip4  slip5  slip6 slip7  slip8
# True          0.0300 0.0100 0.0200 0.0300 0.030 0.0300
# 5 Estimate    0.0021 0.0052 0.0044 0.0069 0.015 0.0032
# 50 Estimate   0.5187 0.4914 0.5070 0.6978 0.688 0.6567
# 97.5 Estimate 0.5488 0.5230 0.5385 0.7287 0.717 0.6860
# Median-True   0.4887 0.4814 0.4870 0.6678 0.658 0.6267
