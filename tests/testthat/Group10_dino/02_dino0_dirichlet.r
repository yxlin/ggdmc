#!/usr/bin/env Rscript
# Baseline CDM Model: DINO Rule Dirchlet
# Subject-level parameter estimation with N=2000
cat("\n\n--------- DINO Baseline Model (Subject-level) ----------\n")
rm(list = ls())
# q(save = "no")

pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group10_dino/data")
fig_dir <- file.path(home_dir, "Group10_dino/figs")
save_path <- paste0(home_dir, "Group10_dino/data/dino0_dirichlet.rda")

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


# -------------------- Model Definition --------------------
simulation_model <- BuildModel(
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
sub_model <- setCDM(simulation_model,
    q_matrix = simulation_model@cdm_info$q_matrix,
    rule = "DINO",
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
pi_params <- grepl("^pi_", simulation_model@pnames)
cat("\nProfile probability parameters:\n")
print(simulation_model@pnames[pi_params])

# Set p0: alpha MUST be > 0 for Dirichlet
p0 <- rep(0, simulation_model@npar)
names(p0) <- simulation_model@pnames
p0[pi_params] <- 1.0 # ✓ FIX: Set alpha = 1 (or 2, or any positive value)

# Set p1: NA for Dirichlet
p1 <- rep(1.0, simulation_model@npar)
names(p1) <- simulation_model@pnames
p1[pi_params] <- NA # ✓ This part was correct

p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = p1,
    lower = rep(NA, simulation_model@npar),
    upper = rep(NA, simulation_model@npar),
    dists = ifelse(pi_params, "dirichlet", "unif"), # Note: 'dists' not 'dist'!
    log_p = rep(TRUE, simulation_model@npar)
)
print_prior(p_prior)
sub_priors <- set_priors(p_prior = p_prior)


# Build data-model-info objects
sub_dmis <- BuildDMI(dat$responses, simulation_model,
    q_matrix = simulation_model@cdm_info$q_matrix,
    rule = "DINO",
    use_mvn = use_mvn
)


# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration plus blocking at the
# subject level (exploration phase)
save_path <- file.path(data_dir, "01_dino0_dirichlet.rda")

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 2, is_pblocked = TRUE,
    seed = 9032
)
save(fits0, file = save_path)

# Stage 1: Restart without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
)

save(fits0, fits1, sim_p_vector, file = save_path)

fits <- fits1
fit <- RebuildPosterior(fits)
# -------------------- Diagnostics (Optional) --------------------
# Check Stage 0: Burn-in chains
figure_name <- file.path(fig_dir, "01_dino0_dirichlet.pdf")
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

#               guess1 guess2 guess3  pi_00  pi_01 pi_10  slip1  slip2   slip3
# True           0.100  0.200  0.300  0.250  0.250 0.250 0.0100 0.0300  0.0500
# 5 Estimate     0.014  0.168  0.011  0.152  0.169 0.220 0.0078 0.0081  0.0242
# 50 Estimate    0.065  0.218  0.113  0.182  0.225 0.288 0.0789 0.0819  0.0427
# 97.5 Estimate  0.134  0.282  0.282  0.235  0.294 0.363 0.1880 0.1945  0.0638
# Median-True   -0.035  0.018 -0.187 -0.068 -0.025 0.038 0.0689 0.0519 -0.0073

#               guess1 guess2 guess3 pi_00  pi_01 pi_10  slip1  slip2   slip3
# True           0.100  0.200  0.300  0.25  0.250 0.250 0.0100 0.0300  0.0500
# 5 Estimate     0.014  0.169  0.011  0.15  0.171 0.217 0.0082 0.0096  0.0239
# 50 Estimate    0.063  0.217  0.107  0.18  0.227 0.285 0.0785 0.0879  0.0431
# 97.5 Estimate  0.135  0.279  0.289  0.23  0.295 0.363 0.1863 0.1951  0.0635
# Median-True   -0.037  0.017 -0.193 -0.07 -0.023 0.035 0.0685 0.0579 -0.0069
