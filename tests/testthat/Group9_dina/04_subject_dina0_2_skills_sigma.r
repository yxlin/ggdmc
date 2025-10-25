#!/usr/bin/env Rscript
# Baseline CDM Model: DINA Rule without MVN
# Subject-level parameter estimation with N=2000
cat("\n\n--------- CDM DINA Baseline Model (Subject-level) ----------\n")
rm(list = ls())
# q(save = "no")

# Load packages
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group9_gen_cdm/data")
fig_dir <- file.path(home_dir, "Group9_gen_cdm/figs")


# -------------------- Q-Matrix Setup --------------------
# 3 items, 2 skills: Item 1 (Algebra only), Item 2 (Geometry only), Item 3 (both)
Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("Algebra", "Geometry")
rownames(Q) <- c("Item1", "Item2", "Item3")

# Calculate skill probabilities (no correlation, sigma = 0)
K <- ncol(Q)
sigma <- .2
means <- rep(0, K)
Sigma <- cdModel::build_correlation_matrix(K, sigma)
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)


# -------------------- Model Definition --------------------
simulation_model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        mean1 = "1", mean2 = "1", simga = "1",
        slip1 = "1", slip2 = "1", slip3 = "1"
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
    rule = "DINA",
    use_mvn = use_mvn
)

# -------------------- Data Simulation --------------------
# True parameter values
sim_p_vector <- c(
    guess1 = .10, guess2 = .20, guess3 = .30,
    mean1 = means[1], mean2 = means[2],
    sigma = sigma,
    slip1 = .01, slip2 = .03, slip3 = .05
)
true_p_vector <- c(
    guess1 = .10, guess2 = .20, guess3 = .30,
    pi_00 = skill_probs$probability[1],
    pi_10 = skill_probs$probability[2],
    pi_01 = skill_probs$probability[3],
    slip1 = .01, slip2 = .03, slip3 = .05
)
N <- 2000
dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = sim_p_vector,
    nschool = 1,
    debug = FALSE,
    seed = 123
)

# -------------------- Prior Setup --------------------
fit_model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        pi_00 = "1", pi_10 = "1", pi_01 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

sub_priors <- cdModel::setup_cdm_prior(fit_model)
print_prior(sub_priors@p_prior)


# Build data-model-info objects
sub_dmis <- BuildDMI(dat$responses, fit_model,
    q_matrix = fit_model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = FALSE
)

# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration plus blocking at the
# subject level (exploration phase)
save_path <- file.path(data_dir, "04_subject_dina0_2_skills_sigma.rda")


fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 2, is_pblocked = TRUE,
    seed = 9032
)

# Stage 1: First restart without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
)

save(fits0, fits1, file = save_path)

fits <- fits1
fit <- RebuildPosterior(fits)

# -------------------- Diagnostics (Optional) --------------------
# Check Stage 0: Burn-in chains

figure_name <- file.path(fig_dir, "04_subject_dina0_2_skills_sigma.pdf")
pdf(figure_name)


p0 <- ggdmc::plot(fits0[[1]], start = 1)
p0 <- ggdmc::plot(fits0[[2]], start = 1)
p0 <- ggdmc::plot(fits0[[3]], start = 1)

p0 <- ggdmc::plot(fits0[[1]], start = 200)
p0 <- ggdmc::plot(fits0[[2]], start = 200)
p0 <- ggdmc::plot(fits0[[3]], start = 200)

p0 <- ggdmc::plot(fits1[[1]])
p0 <- ggdmc::plot(fits1[[2]])
p0 <- ggdmc::plot(fits1[[3]])

p1 <- ggdmc::plot(fit)
p1 <- ggdmc::plot(fit, den = TRUE, pll = FALSE, hide_legend = FALSE)
dev.off()

# -------------------- Convergence Check --------------------
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
est_theta <- ggdmc::compare(fit, ps = true_p_vector)
# Example output (N=2000):
#                guess1 guess2 guess3  pi_00  pi_01 pi_10  slip1 slip2    slip3
# True           0.1000  0.200  0.300  0.282  0.218 0.218 0.0100 0.030  0.05000
# 5 Estimate     0.0089  0.012  0.239  0.181  0.152 0.177 0.0056 0.012  0.00454
# 50 Estimate    0.0962  0.111  0.277  0.225  0.204 0.235 0.0411 0.053  0.04959
# 97.5 Estimate  0.2298  0.258  0.314  0.274  0.273 0.304 0.0998 0.107  0.13698
# Median-True   -0.0038 -0.089 -0.023 -0.057 -0.014 0.017 0.0311 0.023 -0.00041

#                 guess1  guess2 guess3  pi_00 pi_01  pi_10  slip1 slip2   slip3
# True           1.0e-01  0.2000  0.300  0.282 0.218  0.218 0.0100 0.030  0.0500
# 5 Estimate     9.5e-03  0.0099  0.238  0.181 0.178  0.153 0.0056 0.011  0.0044
# 50 Estimate    1.0e-01  0.1150  0.277  0.225 0.235  0.205 0.0412 0.053  0.0470
# 97.5 Estimate  2.3e-01  0.2564  0.314  0.273 0.304  0.271 0.1006 0.107  0.1337
# Median-True   -1.9e-07 -0.0850 -0.023 -0.057 0.017 -0.013 0.0312 0.023 -0.0030
