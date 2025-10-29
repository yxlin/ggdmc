#!/usr/bin/env Rscript
# Baseline CDM Model: DINA Rule using Dirichlet
# Subject-level parameter estimation with N=2000
cat("\n\n----- CDM DINA0 using Dirichlet Prior (Subject-level) -----\n")
rm(list = ls())
# q(save = "no")
# Load packages
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group9_dina/data")
fig_dir <- file.path(home_dir, "Group9_dina/figs")

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

p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = p1,
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dists = ifelse(pi_params, "dirichlet", "unif"), # Note: 'dists' not 'dist'!
    log_p = rep(TRUE, model@npar)
)
print_prior(p_prior)
sub_priors <- set_priors(p_prior = p_prior)

# --------------------------------------------------------------------
# | Parameter |   p0   |   p1   | Lower  | Upper  | Log_P | Dist  |
# --------------------------------------------------------------------
# | guess1 |   0.00 |   1.00 |   0.00 |   1.00 |  True | uniform |
# | guess2 |   0.00 |   1.00 |   0.00 |   1.00 |  True | uniform |
# | guess3 |   0.00 |   1.00 |   0.00 |   1.00 |  True | uniform |
# |  pi_00 |   1.00 |    nan |   0.00 |   1.00 |  True | dirichlet |
# |  pi_01 |   1.00 |    nan |   0.00 |   1.00 |  True | dirichlet |
# |  pi_10 |   1.00 |    nan |   0.00 |   1.00 |  True | dirichlet |
# |  slip1 |   0.00 |   1.00 |   0.00 |   1.00 |  True | uniform |
# |  slip2 |   0.00 |   1.00 |   0.00 |   1.00 |  True | uniform |
# |  slip3 |   0.00 |   1.00 |   0.00 |   1.00 |  True | uniform |
# --------------------------------------------------------------------
# [1] "guess1" "guess2" "guess3" "pi_00"  "pi_01"  "pi_10"  "slip1"  "slip2"
# [9] "slip3"

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

# Stage 1: First restart without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
)

save_path <- file.path(data_dir, "02_subject_dina0_dirichlet.rda")
save(fits0, fits1, file = save_path)
# load(save_path)

fits <- fits1
fit <- RebuildPosterior(fits)
# -------------------- Diagnostics (Optional) --------------------
# Check Stage 0: Burn-in chains
figure_name <- file.path(fig_dir, "02_subject_dina0_dirichlet.pdf")
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
# Uniform Prior
#                guess1 guess2 guess3  pi_00  pi_01 pi_10 slip1  slip2  slip3
# True           0.1000  0.200  0.300  0.250  0.250 0.250 0.010 0.0300 0.0500
# 2.5%           0.0101  0.011  0.234  0.145  0.188 0.211 0.006 0.0092 0.0054
# 50%            0.0952  0.102  0.276  0.187  0.239 0.265 0.047 0.0564 0.0599
# 97.5%          0.2273  0.246  0.312  0.235  0.311 0.339 0.124 0.1269 0.1579
# Median-True   -0.0048 -0.098 -0.024 -0.063 -0.011 0.015 0.037 0.0264 0.0099
# Dirichlet Prior
#                guess1 guess2 guess3  pi_00  pi_01 pi_10 slip1  slip2  slip3
# True           0.1000  0.200  0.300  0.250  0.250 0.250 0.010 0.0300 0.0500
# 5 Estimate     0.0099  0.010  0.233  0.144  0.190 0.212 0.006 0.0097 0.0064
# 50 Estimate    0.0966  0.107  0.276  0.189  0.239 0.264 0.048 0.0574 0.0591
# 97.5 Estimate  0.2231  0.242  0.312  0.236  0.308 0.335 0.122 0.1276 0.1580
# Median-True   -0.0034 -0.093 -0.024 -0.061 -0.011 0.014 0.038 0.0274 0.0091

#                guess1 guess2 guess3  pi_00 pi_01   pi_10  slip1  slip2  slip3
# True           0.1000  0.200  0.300  0.250 0.250  0.2500 0.0100 0.0300 0.0500
# 5 Estimate     0.0090  0.011  0.233  0.143 0.210  0.1906 0.0064 0.0095 0.0054
# 50 Estimate    0.0945  0.107  0.277  0.189 0.263  0.2408 0.0476 0.0559 0.0603
# 97.5 Estimate  0.2258  0.248  0.313  0.237 0.331  0.3099 0.1238 0.1264 0.1552
# Median-True   -0.0055 -0.093 -0.023 -0.061 0.013 -0.0092 0.0376 0.0259 0.0103
