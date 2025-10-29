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
data_dir <- file.path(home_dir, "Group9_dina/data")
fig_dir <- file.path(home_dir, "Group9_dina/figs")
save_path <- file.path(data_dir, "01_subject_dina0.rda")
figure_name <- file.path(fig_dir, "01_subject_dina0.pdf")

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
sigma <- 0
means <- rep(0, K)
Sigma <- cdModel::build_correlation_matrix(K, sigma)
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)
# skill_probs
#   profile probability
# 1   (0,0)        0.25
# 2   (1,0)        0.25
# 3   (0,1)        0.25
# 4   (1,1)        0.25
# -------------------- Model Definition --------------------
model <- BuildModel(
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
# print_prior(p_prior)


# Build data-model-info objects
sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = use_mvn
)
# sub_dmis[[1]]@use_mvn
# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration no blocking at the
# e subject level (exploration phase)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 2, is_pblocked = TRUE,
    seed = 9032
)
save(fits0, file = save_path)
# Stage 1: Sampling without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
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
p1 <- ggdmc::plot(fit, den = TRUE, pll = FALSE)
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

# Example output (N=2000):
#                guess1 guess2 guess3  pi_00  pi_01 pi_10 slip1  slip2  slip3
# True           0.1000  0.200  0.300  0.250  0.250 0.250 0.010 0.0300 0.0500
# 2.5%           0.0101  0.011  0.234  0.145  0.188 0.211 0.006 0.0092 0.0054
# 50%            0.0952  0.102  0.276  0.187  0.239 0.265 0.047 0.0564 0.0599
# 97.5%          0.2273  0.246  0.312  0.235  0.311 0.339 0.124 0.1269 0.1579
# Median-True   -0.0048 -0.098 -0.024 -0.063 -0.011 0.015 0.037 0.0264 0.0099

#                guess1 guess2 guess3 pi_00 pi_01 pi_10  slip1 slip2  slip3
# True           0.1000  0.200  0.300  0.25 0.250  0.25 0.0100 0.030 0.0500
# 5 Estimate     0.0097  0.012  0.231  0.15 0.212  0.19 0.0061 0.010 0.0050
# 50 Estimate    0.0975  0.110  0.276  0.19 0.264  0.24 0.0467 0.057 0.0574
# 97.5 Estimate  0.2253  0.248  0.311  0.24 0.335  0.31 0.1217 0.131 0.1532
# Median-True   -0.0025 -0.090 -0.024 -0.06 0.014 -0.01 0.0367 0.027 0.0074

#                guess1 guess2 guess3 pi_00 pi_01   pi_10  slip1 slip2  slip3
# True           0.1000  0.200  0.300  0.25 0.250  0.2500 0.0100 0.030 0.0500
# 5 Estimate     0.0098  0.011  0.232  0.15 0.211  0.1889 0.0064 0.010 0.0055
# 50 Estimate    0.0978  0.110  0.275  0.19 0.264  0.2401 0.0486 0.058 0.0561
# 97.5 Estimate  0.2275  0.243  0.311  0.24 0.336  0.3060 0.1227 0.129 0.1551
# Median-True   -0.0022 -0.090 -0.025 -0.06 0.014 -0.0099 0.0386 0.028 0.0061

#     guess1 guess2 guess3  pi_00 pi_01  pi_10 slip1 slip2 slip3
# True          0.1000  0.200  0.300  0.250 0.250  0.250 0.010 0.030 0.050
# 5 Estimate    0.0094  0.011  0.233  0.145 0.210  0.185 0.005 0.010 0.006
# 50 Estimate   0.1069  0.108  0.277  0.192 0.267  0.238 0.047 0.061 0.065
# 97.5 Estimate 0.2412  0.258  0.314  0.242 0.342  0.315 0.124 0.133 0.174
# Median-True   0.0069 -0.092 -0.023 -0.058 0.017 -0.012 0.037 0.031 0.015
