#!/usr/bin/env Rscript
# Baseline CDM Model: DINA Rule without MVN
# Subject-level parameter estimation ecpe
cat("\n\n--------- CDM DINA Baseline Model (Subject-level) ----------\n")
rm(list = ls())
# q(save = "no")
# Load packages
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group9_gen_cdm/data")
fig_dir <- file.path(home_dir, "Group9_gen_cdm/figs")
save_path <- file.path(data_dir, "07_subject_dina0_ecpe.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))

load(save_path)
# -------------------- Q-Matrix Setup --------------------
# 3 items, 2 skills: Item 1 (Algebra only), Item 2 (Geometry only), Item 3 (both)
K <- 3
ecpe_dat <- CDM::data.ecpe$data[, c("id", paste0("E", 1:10))]
Q <- CDM::data.ecpe$q.matrix[1:10, ]

ecpe <- CDM::din(data = ecpe_dat[, -1], q.matrix = Q)
param <- CDM::IRT.se(ecpe, extended = TRUE)
p <- split(param, param$partype)


cdm_result <- tibble::tibble(
    guess = p$guess$est,
    slip = p$slip$est,
    guess005 = p$guess$est - 1.96 * p$guess$se,
    guess975 = p$guess$est + 1.96 * p$guess$se,
    slip005 = p$slip$est - 1.96 * p$slip$se,
    slip975 = p$slip$est + 1.96 * p$slip$se
)

# Rearrange to match ggdmc format
nitem <- length(p$guess$est)
cdm_formatted <- data.frame(
    `5%` = c(
        p$guess$est - 1.96 * p$guess$se, # guess lower bounds
        p$slip$est - 1.96 * p$slip$se # slip lower bounds
    ),
    `50%` = c(p$guess$est, p$slip$est), # medians (estimates)
    `97.5%` = c(
        p$guess$est + 1.96 * p$guess$se, # guess upper bounds
        p$slip$est + 1.96 * p$slip$se # slip upper bounds
    ),
    check.names = FALSE
)
rownames(cdm_formatted) <- c(
    paste0("guess", 1:nitem),
    paste0("slip", 1:nitem)
)

cat("\n=== CDM Package Results ===\n")
print(round(cdm_formatted, 3))

cat("\n=== Original cdm_result ===\n")
print(cdm_result)


Qmat <- as.matrix(Q)
colnames(Qmat) <- c("Algebra", "Geometry", "Science")
rownames(Qmat) <- paste0("Item", 1:nrow(Qmat))

# -------------------- Prior Setup --------------------
p_map <- list(
    guess1 = "1", guess2 = "1", guess3 = "1",
    guess4 = "1", guess5 = "1", guess6 = "1",
    guess7 = "1", guess8 = "1",
    guess9 = "1", guess10 = "1",
    pi_000 = "1",
    pi_100 = "1",
    pi_010 = "1",
    pi_110 = "1",
    pi_001 = "1",
    pi_101 = "1",
    pi_011 = "1",
    slip1 = "1", slip2 = "1", slip3 = "1",
    slip4 = "1", slip5 = "1", slip6 = "1",
    slip7 = "1", slip8 = "1",
    slip9 = "1", slip10 = "1"
)

fit_model <- BuildModel(
    p_map = p_map,
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Qmat,
    type = "cdm",
    verbose = TRUE
)

sub_priors <- cdModel::setup_cdm_prior(fit_model)
print_prior(sub_priors@p_prior)



dat_long <- wide2long_base(ecpe_dat)
dat_long$s <- 1

dat_responses <- dat_long[, 2:5]

# Build data-model-info objects
sub_dmis <- BuildDMI(dat_responses, fit_model,
    q_matrix = fit_model@cdm_info$q_matrix,
    rule = "DINA",
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


# Stage 1: First restart without migration
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE
)

save(fits0, fits1, file = save_path)

fits <- fits1
fit <- RebuildPosterior(fits)

# -------------------- Diagnostics (Optional) --------------------
# Check Stage 0: Burn-in chains
figure_name <- file.path(fig_dir, "07_subject_dina0_ecpe.pdf")
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

# -------------------- Parameter summary --------------------
options(digits = 2)
est_theta <- ggdmc::summary(fit)
print(est_theta$quantile)

cdm_formatted[1:10, ] - est_theta$quantile[1:10, ]
#               5%     50%   97.5%
# guess1  -0.01014 -0.0089 -0.0103
# guess2  -0.03336 -0.0431 -0.0431
# guess3  -0.00827 -0.0057 -0.0095
# guess4   0.00198 -0.0049 -0.0141
# guess5   0.00016 -0.0031 -0.0072
# guess6   0.00202 -0.0024 -0.0080
# guess7  -0.00568 -0.0089 -0.0153
# guess8  -0.03162 -0.0429 -0.0440
# guess9  -0.00278 -0.0029 -0.0066
# guess10  0.00097 -0.0184 -0.0293

cdm_formatted[11:20, ] - est_theta$quantile[18:27, ]
#              5%     50%    97.5%
# slip1   0.01135 0.00816  0.00200
# slip2   0.01125 0.00829  0.00296
# slip3   0.00551 0.00513 -0.00041
# slip4   0.00287 0.00135 -0.00339
# slip5  -0.00125 0.00017 -0.00111
# slip6  -0.00069 0.00056 -0.00128
# slip7   0.00678 0.00480 -0.00260
# slip8   0.00907 0.00657  0.00221
# slip9   0.00084 0.00131 -0.00148
# slip10  0.00702 0.00436 -0.00154

#    [,1]  [,2]  [,3]
# 000  0.1286 0.175 0.220
# 100 -0.0126 0.021 0.055
# 010  0.0360 0.078 0.119
# 001  0.0304 0.065 0.099
# 110 -0.0050 0.023 0.051
# 101 -0.0036 0.018 0.039
# 011  0.0675 0.104 0.141
# 111  0.4856 0.517 0.548


cdm_est_profile_probs <- cbind(
    p$probs$est - 1.96 * p$probs$se, p$probs$est,
    p$probs$est + 1.96 * p$probs$se
)
rownames(cdm_est_profile_probs) <- rownames(ecpe$attribute.patt)
cdm_est_profile_probs

ggdmc_est <- est_theta$quantile[11:17, ]
ggdmc_est


# Full comparison with output
result <- compare_profile_probs(ggdmc_est, cdm_est_profile_probs)
# === Profile Probability Comparison ===

# Matrix 1 (ggdmc_est):
#         5%   50% 97.5%
# 000 0.1579 0.220 0.292
# 100 0.0012 0.015 0.059
# 010 0.0045 0.040 0.109
# 110 0.0022 0.018 0.052
# 001 0.0460 0.092 0.160
# 101 0.0053 0.032 0.084
# 011 0.0326 0.086 0.147

# Matrix 2 (cdm_est_profile_probs):
#        [,1]  [,2]  [,3]
# 000  0.1286 0.175 0.220
# 100 -0.0126 0.021 0.055
# 010  0.0360 0.078 0.119
# 110 -0.0050 0.023 0.051
# 001  0.0304 0.065 0.099
# 101 -0.0036 0.018 0.039
# 011  0.0675 0.104 0.141

# Differences (Matrix 1 - Matrix 2):
#          5%     50%   97.5%
# 000  0.0293  0.0457  0.0716
# 100  0.0137 -0.0064  0.0049
# 010 -0.0315 -0.0373 -0.0107
# 110  0.0072 -0.0048  0.0007
# 001  0.0156  0.0271  0.0606
# 101  0.0089  0.0143  0.0454
# 011 -0.0349 -0.0183  0.0059

# Absolute differences:
#         5%    50%  97.5%
# 000 0.0293 0.0457 0.0716
# 100 0.0137 0.0064 0.0049
# 010 0.0315 0.0373 0.0107
# 110 0.0072 0.0048 0.0007
# 001 0.0156 0.0271 0.0606
# 101 0.0089 0.0143 0.0454
# 011 0.0349 0.0183 0.0059

# Max absolute difference: 0.072
# Mean absolute difference: 0.024
# Or just get differences
diff <- quick_compare(ggdmc_est, cdm_est_profile_probs)
diff
#          5%     50%    97.5%
# 000  0.0293  0.0457  0.07160
# 100  0.0137 -0.0064  0.00488
# 010 -0.0315 -0.0373 -0.01073
# 110  0.0072 -0.0048  0.00069
# 001  0.0156  0.0271  0.06063
# 101  0.0089  0.0143  0.04544
# 011 -0.0349 -0.0183  0.00594
