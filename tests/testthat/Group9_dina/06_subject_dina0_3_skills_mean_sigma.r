#!/usr/bin/env Rscript
# Baseline CDM Model: DINA Rule without MVN
# Subject-level parameter estimation with N=3000
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
    rule = "DINA",
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
    rule = "DINA",
    use_mvn = FALSE
)

# -------------------- MCMC Sampling --------------------
# Stage 0: Burn-in with migration plus blocking at the
# subject level (exploration phase)
save_path <- file.path(data_dir, "06_subject_dina0_3_skills_mean_sigam.rda")

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
# Check Stage 0: Burn-in chains
figure_name <- file.path(fig_dir, "06_subject_dina0_3_skills_mean_sigma.pdf")
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
#               guess1 guess2 guess3 guess4  guess5 guess6 guess7  guess8  pi_000
# True          0.1000  0.010 0.0300 0.0300  0.0300 0.0400 0.0400  0.0500  0.1038
# 5 Estimate    0.0897  0.053 0.0181 0.0211  0.0075 0.0374 0.0313  0.0334  0.0839
# 50 Estimate   0.1029  0.505 0.0318 0.0355  0.0230 0.0483 0.0419  0.0428  0.0961
# 97.5 Estimate 0.9662  0.972 0.9985 0.9985  0.9991 0.9808 0.9855  0.9591  0.2583
# Median-True   0.0029  0.495 0.0018 0.0055 -0.0070 0.0083 0.0019 -0.0072 -0.0077
#               pi_001  pi_010   pi_011   pi_100   pi_101   pi_110  slip1   slip2
# True          0.1304  0.0753  0.15063  0.06874  0.13737  0.07850 0.0500  0.0200
# 5 Estimate    0.0016  0.0015  0.00017  0.00120  0.00027  0.00022 0.0385  0.0131
# 50 Estimate   0.1321  0.0678  0.14041  0.06865  0.13380  0.07756 0.0542  0.0166
# 97.5 Estimate 0.1490  0.0807  0.15732  0.08137  0.14989  0.09091 0.9116  0.0214
# Median-True   0.0016 -0.0075 -0.01022 -0.00009 -0.00357 -0.00094 0.0042 -0.0034
#                slip3   slip4   slip5  slip6   slip7  slip8
# True          0.0200  0.0200  0.0200 0.0300  0.0300 0.0500
# 5 Estimate    0.0153  0.0103  0.0108 0.0215  0.0163 0.0441
# 50 Estimate   0.0252  0.0184  0.0187 0.0329  0.0267 0.0598
# 97.5 Estimate 0.6024  0.6047  0.4440 0.7919  0.7870 0.8691
# Median-True   0.0052 -0.0016 -0.0013 0.0029 -0.0033 0.0098

#               guess1 guess2 guess3 guess4  guess5 guess6 guess7  guess8 pi_000
# True          0.1000  0.010 0.0300 0.0300  0.0300 0.0400 0.0400  0.0500  0.104
# 5 Estimate    0.0898  0.053 0.0182 0.0214  0.0077 0.0376 0.0318  0.0333  0.084
# 50 Estimate   0.1029  0.508 0.0319 0.0361  0.0232 0.0483 0.0423  0.0427  0.096
# 97.5 Estimate 0.9665  0.975 0.9987 0.9988  0.9989 0.9813 0.9852  0.9588  0.259
# Median-True   0.0029  0.498 0.0019 0.0061 -0.0068 0.0083 0.0023 -0.0073 -0.008
#               pi_001  pi_010   pi_011   pi_100   pi_101   pi_110  slip1   slip2
# True          0.1305  0.0753  0.15064  0.06875  0.13745  0.07853 0.0500  0.0200
# 5 Estimate    0.0016  0.0016  0.00018  0.00120  0.00026  0.00023 0.0388  0.0131
# 50 Estimate   0.1325  0.0678  0.14065  0.06856  0.13383  0.07763 0.0542  0.0166
# 97.5 Estimate 0.1493  0.0806  0.15700  0.08196  0.14981  0.09012 0.9107  0.0217
# Median-True   0.0021 -0.0075 -0.00999 -0.00019 -0.00361 -0.00090 0.0042 -0.0034
#                slip3   slip4   slip5  slip6   slip7  slip8
# True          0.0200  0.0200  0.0200 0.0300  0.0300 0.0500
# 5 Estimate    0.0152  0.0101  0.0107 0.0212  0.0164 0.0434
# 50 Estimate   0.0247  0.0185  0.0185 0.0328  0.0267 0.0598
# 97.5 Estimate 0.6015  0.6047  0.4443 0.7914  0.7870 0.8706
# Median-True   0.0047 -0.0015 -0.0015 0.0028 -0.0033 0.0098

est_theta <- ggdmc::compare(fits[[2]], ps = true_p_vector)
#                guess1 guess2  guess3  guess4  guess5 guess6   guess7 guess8
# True           0.1000  0.010  0.0300 0.03000  0.0300 0.0400  0.04000  0.050
# 5 Estimate     0.0885  0.046  0.0172 0.02002  0.0066 0.0367  0.03101  0.033
# 50 Estimate    0.0982  0.505  0.0272 0.03094  0.0171 0.0446  0.03909  0.040
# 97.5 Estimate  0.1113  0.979  0.0428 0.04650  0.0356 0.0556  0.05025  0.049
# Median-True   -0.0018  0.495 -0.0028 0.00094 -0.0129 0.0046 -0.00091 -0.010
#               pi_000 pi_001  pi_010  pi_011 pi_100  pi_101 pi_110    slip1
# True           0.104 0.1305  0.0753  0.1506 0.0688 0.13745 0.0785  0.05000
# 5 Estimate     0.083 0.1262  0.0630  0.1347 0.0636 0.12807 0.0730  0.03764
# 50 Estimate    0.092 0.1369  0.0712  0.1451 0.0720 0.13837 0.0812  0.04903
# 97.5 Estimate  0.103 0.1504  0.0816  0.1582 0.0830 0.15186 0.0913  0.06522
# Median-True   -0.012 0.0064 -0.0041 -0.0055 0.0033 0.00093 0.0027 -0.00097
#                 slip2  slip3   slip4   slip5   slip6   slip7  slip8
# True           0.0200 0.0200  0.0200  0.0200  0.0300  0.0300 0.0500
# 5 Estimate     0.0132 0.0144  0.0096  0.0100  0.0204  0.0157 0.0414
# 50 Estimate    0.0168 0.0214  0.0156  0.0158  0.0289  0.0232 0.0542
# 97.5 Estimate  0.0216 0.0313  0.0245  0.0242  0.0426  0.0342 0.0706
# Median-True   -0.0032 0.0014 -0.0044 -0.0042 -0.0011 -0.0068 0.0042

est_theta <- ggdmc::compare(fits[[3]], ps = true_p_vector)
#                guess1 guess2 guess3 guess4 guess5 guess6 guess7 guess8 pi_000
# True            0.10  0.010   0.03   0.03   0.03   0.04   0.04   0.05   0.10
# 5 Estimate      0.94  0.064   0.98   0.98   0.98   0.96   0.97   0.93   0.23
# 50 Estimate     0.96  0.529   0.99   0.99   0.99   0.97   0.98   0.95   0.25
# 97.5 Estimate   0.97  0.974   1.00   1.00   1.00   0.98   0.99   0.96   0.26
# Median-True     0.86  0.519   0.96   0.96   0.96   0.93   0.94   0.90   0.14
#                 pi_001   pi_010   pi_011   pi_100   pi_101   pi_110 slip1
# True           0.13046  0.07532  0.15064  0.06875  1.4e-01  0.07853  0.05
# 5 Estimate     0.00058  0.00058  0.00006  0.00039  9.9e-05  0.00007  0.89
# 50 Estimate    0.00465  0.00458  0.00071  0.00401  1.0e-03  0.00090  0.90
# 97.5 Estimate  0.01195  0.01157  0.00363  0.01018  4.2e-03  0.00365  0.91
# Median-True   -0.12581 -0.07074 -0.14993 -0.06475 -1.4e-01 -0.07763  0.85
#                 slip2 slip3 slip4 slip5 slip6 slip7 slip8
# True           0.0200  0.02  0.02  0.02  0.03  0.03  0.05
# 5 Estimate     0.0131  0.57  0.57  0.41  0.77  0.76  0.85
# 50 Estimate    0.0166  0.59  0.59  0.43  0.78  0.77  0.86
# 97.5 Estimate  0.0218  0.61  0.61  0.45  0.80  0.79  0.87
# Median-True   -0.0034  0.57  0.57  0.41  0.75  0.74  0.81
