#!/usr/bin/env Rscript
# Baseline CDM Model: DINA Rule 3 skills using Dirichlet
# Subject-level parameter estimation with N=3000
cat("\n\n----- CDM DINA0 using 3 skill Dirichlet Prior (Subject-level) -----\n")
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
    rule = "DINA",
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
    rule = "DINA",
    use_mvn = use_mvn
)

# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration plus blocking at the
# subject level (exploration phase)
save_path <- file.path(data_dir, "03_subject_dina0_3_skills_dirichlet.rda")

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 4, is_pblocked = FALSE,
    seed = 9032
)
save(fits0, file = save_path)

# Stage 1: First restart without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 4, is_pblocked = FALSE
)
save(fits0, fits1, file = save_path)

fits <- fits1
fit <- RebuildPosterior(fits)
# -------------------- Diagnostics (Optional) --------------------
# Check Stage 0: Burn-in chains
figure_name <- file.path(fig_dir, "03_subject_dina0_3_skills_dirichlet.pdf")
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



est_theta <- ggdmc::compare(fits[[1]], ps = sim_p_vector)
#                guess1 guess2  guess3 guess4 guess5  guess6 guess7  guess8
# True          0.10000  0.010  0.0300 0.0300  0.030 0.04000 0.0400  0.0500
# 5 Estimate    0.09185  0.046  0.0118 0.0210  0.017 0.03342 0.0349  0.0363
# 50 Estimate   0.10095  0.485  0.0222 0.0326  0.027 0.04071 0.0419  0.0434
# 97.5 Estimate 0.11281  0.976  0.0390 0.0492  0.041 0.05011 0.0513  0.0525
# Median-True   0.00095  0.475 -0.0078 0.0026 -0.003 0.00071 0.0019 -0.0066
#                pi_000 pi_001 pi_010  pi_011  pi_100 pi_101 pi_110  slip1
# True           0.1250 0.1250 0.1250  0.1250  0.1250   0.12 0.1250  0.050
# 5 Estimate     0.1053 0.1168 0.1154  0.1072  0.1127   0.12 0.1231  0.018
# 50 Estimate    0.1159 0.1274 0.1261  0.1171  0.1232   0.14 0.1333  0.032
# 97.5 Estimate  0.1292 0.1407 0.1393  0.1294  0.1365   0.15 0.1461  0.055
# Median-True   -0.0091 0.0024 0.0011 -0.0079 -0.0018   0.01 0.0083 -0.018
#                 slip2   slip3   slip4  slip5   slip6 slip7  slip8
# True           0.0200  0.0200  0.0200 0.0200  0.0300 0.030 0.0500
# 5 Estimate     0.0133  0.0083  0.0092 0.0170  0.0112 0.029 0.0409
# 50 Estimate    0.0168  0.0156  0.0171 0.0271  0.0203 0.044 0.0568
# 97.5 Estimate  0.0217  0.0276  0.0297 0.0414  0.0362 0.065 0.0795
# Median-True   -0.0032 -0.0044 -0.0029 0.0071 -0.0097 0.014 0.0068

est_theta <- ggdmc::compare(fits[[2]], ps = sim_p_vector)
#    guess1 guess2  guess3 guess4 guess5  guess6 guess7  guess8 pi_000
# True          0.1000   0.01  0.0300 0.0300  0.030 0.04000  0.040  0.0500  0.125
# 5 Estimate    0.0916   0.05  0.0120 0.0208  0.017 0.03317  0.035  0.0361  0.105
# 50 Estimate   0.1011   0.50  0.0226 0.0325  0.027 0.04052  0.042  0.0431  0.116
# 97.5 Estimate 0.1129   0.98  0.0390 0.0480  0.042 0.05010  0.052  0.0526  0.130
# Median-True   0.0011   0.49 -0.0074 0.0025 -0.003 0.00052  0.002 -0.0069 -0.009
#               pi_001 pi_010  pi_011  pi_100 pi_101 pi_110  slip1   slip2
# True           0.125 0.1250  0.1250  0.1250   0.12 0.1250  0.050  0.0200
# 5 Estimate     0.116 0.1160  0.1069  0.1125   0.12 0.1229  0.018  0.0133
# 50 Estimate    0.127 0.1265  0.1171  0.1234   0.14 0.1333  0.032  0.0169
# 97.5 Estimate  0.141 0.1403  0.1297  0.1368   0.15 0.1461  0.056  0.0222
# Median-True    0.002 0.0015 -0.0079 -0.0016   0.01 0.0083 -0.018 -0.0031
#                 slip3   slip4  slip5   slip6 slip7  slip8
# True           0.0200  0.0200 0.0200  0.0300 0.030 0.0500
# 5 Estimate     0.0081  0.0092 0.0165  0.0113 0.029 0.0404
# 50 Estimate    0.0155  0.0171 0.0271  0.0206 0.044 0.0569
# 97.5 Estimate  0.0277  0.0291 0.0423  0.0351 0.067 0.0792
# Median-True   -0.0045 -0.0029 0.0071 -0.0094 0.014 0.0069

est_theta <- ggdmc::compare(fits[[3]], ps = sim_p_vector)
