#!/usr/bin/env Rscript
# Baseline CDM Model: DINO Rule without MVN
# Subject-level parameter estimation with N=2000
cat("\n\n--------- CDM DINO Baseline Model (Subject-level) ----------\n")
rm(list = ls())
# q(save = "no")

# Load packages
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group10_dino/data")
fig_dir <- file.path(home_dir, "Group10_dino/figs")
save_path <- file.path(data_dir, "03_dino0_2_skills_mean_sigma.rda")

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
means <- c(.1, .5)
Sigma <- cdModel::build_correlation_matrix(K, sigma)
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)
skill_probs

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
    rule = "DINO",
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

sub_priors <- cdModel::setup_cdm_prior(fit_model, verbose = TRUE)
print_prior(sub_priors@p_prior)

# Build data-model-info objects
sub_dmis <- BuildDMI(dat$responses, fit_model,
    q_matrix = fit_model@cdm_info$q_matrix,
    rule = "DINO",
    use_mvn = FALSE
)

# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration plus blocking at the
# subject level (exploration phase)
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 2, is_pblocked = TRUE,
    seed = 9032
)
save(fits0, file = save_path)

# Stage 1: Restart without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
)
save(fits0, fits1, true_p_vector, file = save_path)

fits <- fits1
fit <- RebuildPosterior(fits)
# -------------------- Diagnostics (Optional) --------------------
# Check Stage 0: Burn-in chains
figure_name <- file.path(fig_dir, "03_dino0_2_skills_mean_sigma.pdf")
pdf(figure_name)
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
#               guess1 guess2  guess3  pi_00  pi_01 pi_10 slip1 slip2  slip3
# True           0.100  0.200  0.3000  0.170  0.290 0.138 0.010 0.030  0.050
# 5 Estimate     0.032  0.171  0.0099  0.110  0.141 0.090 0.011 0.006  0.026
# 50 Estimate    0.081  0.234  0.1061  0.132  0.227 0.154 0.109 0.059  0.039
# 97.5 Estimate  0.146  0.309  0.2906  0.171  0.313 0.222 0.234 0.139  0.055
# Median-True   -0.019  0.034 -0.1939 -0.038 -0.062 0.015 0.099 0.029 -0.011

est_theta2 <- ggdmc::summary(fit)
est_theta2$statistics
est_theta2$quantile

#            5%   50% 97.5%
# guess1 0.0322 0.081 0.146
# guess2 0.1715 0.234 0.309
# guess3 0.0099 0.106 0.291
# pi_00  0.1104 0.132 0.171
# pi_10  0.0902 0.154 0.222
# pi_01  0.1405 0.227 0.313
# slip1  0.0112 0.109 0.234
# slip2  0.0060 0.059 0.139
# slip3  0.0262 0.039 0.055

true_p_vector
# guess1 guess2 guess3  pi_00  pi_10  pi_01  slip1  slip2  slip3
#   0.10   0.20   0.30   0.17   0.14   0.29   0.01   0.03   0.05
