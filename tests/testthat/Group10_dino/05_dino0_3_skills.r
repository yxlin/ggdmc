#!/usr/bin/env Rscript
# Baseline CDM Model: DINO Rule without MVN
# Subject-level parameter estimation with N=3000
cat("\n\n--------- CDM DINO Baseline Model (Subject-level) ----------\n")
rm(list = ls())
# q(save = "no")
# Load packages
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))


home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group10_dino/data")
fig_dir <- file.path(home_dir, "Group10_dino/figs")
save_path <- file.path(data_dir, "05_dino0_3_skills.rda")
figure_name <- file.path(fig_dir, "05_dino0_3_skills.pdf")

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

# -------------------- Prior Setup --------------------
p_map <- list(
    guess1 = "1", guess2 = "1", guess3 = "1",
    guess4 = "1", guess5 = "1", guess6 = "1",
    guess7 = "1", guess8 = "1",
    pi_000 = "1",
    pi_100 = "1",
    pi_010 = "1",
    pi_110 = "1",
    pi_001 = "1",
    pi_101 = "1",
    pi_011 = "1",
    slip1 = "1", slip2 = "1", slip3 = "1",
    slip4 = "1", slip5 = "1", slip6 = "1",
    slip7 = "1", slip8 = "1"
)

fit_model <- BuildModel(
    p_map = p_map,
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
    rule = "DINO",
    use_mvn = FALSE
)

# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration plus blocking at the
# subject level (exploration phase)
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 8, is_pblocked = FALSE,
    seed = 9032
)
save(fits0, file = save_path)

# Stage 1: First restart without migration
# load(save_path)
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
)

save(fits0, fits1, true_p_vector, file = save_path)

fits <- fits1
fit <- RebuildPosterior(fits)
# -------------------- Diagnostics (Optional) --------------------
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
#               guess1  guess2  guess3 guess4 guess5 guess6 guess7 guess8 pi_000
# True           0.100 0.01000  0.0300 0.0300 0.0300 0.0400  0.040  0.050  0.104
# 5 Estimate     0.050 0.00807  0.0158 0.0214 0.0213 0.0308  0.012  0.013  0.084
# 50 Estimate    0.074 0.01083  0.0265 0.0341 0.0349 0.0494  0.024  0.025  0.093
# 97.5 Estimate  0.108 0.01481  0.0427 0.0512 0.0551 0.0765  0.046  0.043  0.104
# Median-True   -0.026 0.00083 -0.0035 0.0041 0.0049 0.0094 -0.016 -0.025 -0.011
#               pi_001 pi_010  pi_011 pi_100  pi_101 pi_110  slip1 slip2  slip3
# True          0.1305  0.075  0.1506 0.0687  0.1374 0.0785  0.050 0.020 0.0200
# 5 Estimate    0.1234  0.065  0.1295 0.0667  0.1239 0.0749  0.040 0.048 0.0111
# 50 Estimate   0.1340  0.072  0.1419 0.0747  0.1359 0.0844  0.046 0.487 0.0229
# 97.5 Estimate 0.1475  0.082  0.1573 0.0849  0.1501 0.0968  0.055 0.972 0.0404
# Median-True   0.0035 -0.003 -0.0088 0.0059 -0.0015 0.0059 -0.004 0.467 0.0029
#                slip4   slip5   slip6  slip7   slip8
# True          0.0200  0.0200  0.0300  0.030  0.0500
# 5 Estimate    0.0116  0.0066  0.0204  0.021  0.0403
# 50 Estimate   0.0241  0.0140  0.0257  0.027  0.0482
# 97.5 Estimate 0.0427  0.0258  0.0330  0.034  0.0585
# Median-True   0.0041 -0.0060 -0.0043 -0.003 -0.0018
