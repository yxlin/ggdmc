#!/usr/bin/env Rscript
# Baseline CDM Model: DINO Rule 3 skills using Dirichlet
# Subject-level parameter estimation with N=3000
cat("\n\n----- CDM DINO0 using 3 skill Dirichlet Prior (Subject-level) -----\n")
rm(list = ls())
# q(save = "no")
# Load packages
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group10_dino/data")
fig_dir <- file.path(home_dir, "Group10_dino/figs")
save_path <- file.path(data_dir, "04_dino0_3_skills_dirichlet.rda")
figure_name <- file.path(fig_dir, "04_dino0_3_skills_dirichlet.pdf")

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
    rule = "DINO",
    use_mvn = use_mvn
)

# -------------------- Data Simulation --------------------
# True parameter values
sim_p_vector <- c(
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
sub_priors <- cdModel::setup_cdm_prior(model)
print_prior(sub_priors@p_prior)
# sub_priors <- set_priors(p_prior = p_prior)

# Build data-model-info objects
sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINO",
    use_mvn = use_mvn
)

# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration plus blocking at the
# subject level (exploration phase)
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 4, is_pblocked = FALSE,
    seed = 9032
)
save(fits0, file = save_path)

# Stage 1: Restart without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 4, is_pblocked = FALSE
)
save(fits0, fits1, file = save_path)

fits <- fits1
fit <- RebuildPosterior(fits)

# -------------------- Diagnostics (Optional) --------------------
# Check Stage 0: Burn-in chains
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
est_theta <- ggdmc::compare(fit, ps = sim_p_vector)
#               guess1 guess2 guess3  guess4  guess5 guess6  guess7  guess8
# True           0.100 0.0100  0.030  0.0300  0.0300 0.0400 4.0e-02  0.0500
# 5 Estimate     0.059 0.0075  0.010  0.0174  0.0136 0.0297 2.6e-02  0.0272
# 50 Estimate    0.082 0.0102  0.019  0.0271  0.0221 0.0448 4.0e-02  0.0411
# 97.5 Estimate  0.115 0.0141  0.033  0.0414  0.0350 0.0669 6.1e-02  0.0611
# Median-True   -0.018 0.0002 -0.011 -0.0029 -0.0079 0.0048 3.2e-05 -0.0089
#               pi_000  pi_001 pi_010 pi_011  pi_100 pi_101 pi_110   slip1 slip2
# True           0.125  0.1250 0.1250  0.125 0.12500 0.1250  0.125  0.0500 0.020
# 5 Estimate     0.104  0.1133 0.1172  0.104 0.11564 0.1217  0.127  0.0371 0.051
# 50 Estimate    0.113  0.1233 0.1274  0.114 0.12593 0.1325  0.138  0.0434 0.495
# 97.5 Estimate  0.125  0.1362 0.1403  0.127 0.13864 0.1460  0.152  0.0515 0.976
# Median-True   -0.012 -0.0017 0.0024 -0.011 0.00093 0.0075  0.013 -0.0066 0.475
#                slip3 slip4   slip5   slip6  slip7 slip8
# True          0.0200 0.020  0.0200  0.0300 0.0300 0.050
# 5 Estimate    0.0152 0.018  0.0060  0.0186 0.0264 0.047
# 50 Estimate   0.0259 0.029  0.0138  0.0242 0.0328 0.055
# 97.5 Estimate 0.0419 0.045  0.0279  0.0319 0.0415 0.066
# Median-True   0.0059 0.009 -0.0062 -0.0058 0.0028 0.005
