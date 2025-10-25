#!/usr/bin/env Rscript
# DINA Model fit ECPE data
# Subject-level parameter estimation
cat("\n\n--------- CDM DINA ECPE Model (Subject-level) ----------\n")
rm(list = ls())
# q(save = "no")
# Load packages
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_dir <- file.path(home_dir, "Group9_gen_cdm/data")
fig_dir <- file.path(home_dir, "Group9_gen_cdm/figs")
save_path <- file.path(data_dir, "08_subject_dina0_ecpe.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))

load(save_path)
# -------------------- Q-Matrix Setup --------------------
n_student <- dim(CDM::data.ecpe$data)[1]
n_item <- dim(CDM::data.ecpe$data)[2] - 1
K <- dim(CDM::data.ecpe$q.matrix)[2]

ecpe_dat <- CDM::data.ecpe$data
Q <- CDM::data.ecpe$q.matrix

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
    paste0("guess", 1:n_item),
    paste0("slip", 1:n_item)
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
    guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1", guess5 = "1",
    guess6 = "1", guess7 = "1", guess8 = "1", guess9 = "1", guess10 = "1",
    guess11 = "1", guess12 = "1", guess13 = "1", guess14 = "1", guess15 = "1",
    guess16 = "1", guess17 = "1", guess18 = "1", guess19 = "1", guess20 = "1",
    guess21 = "1", guess22 = "1", guess23 = "1", guess24 = "1", guess25 = "1",
    guess26 = "1", guess27 = "1", guess28 = "1",
    pi_000 = "1",
    pi_100 = "1",
    pi_010 = "1",
    pi_110 = "1",
    pi_001 = "1",
    pi_101 = "1",
    pi_011 = "1",
    slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1",
    slip6 = "1", slip7 = "1", slip8 = "1", slip9 = "1", slip10 = "1",
    slip11 = "1", slip12 = "1", slip13 = "1", slip14 = "1", slip15 = "1",
    slip16 = "1", slip17 = "1", slip18 = "1", slip19 = "1", slip20 = "1",
    slip21 = "1", slip22 = "1", slip23 = "1", slip24 = "1", slip25 = "1",
    slip26 = "1", slip27 = "1", slip28 = "1"
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
figure_name <- file.path(fig_dir, "08_subject_dina0_ecpe_all.pdf")
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
# print(est_theta$quantile)

cdm_formatted[1:28, ] - est_theta$quantile[1:28, ]
#               5%      50%    97.5%
# guess1  -0.00462 -0.00080  3.5e-04
# guess2  -0.01650 -0.01221 -9.2e-03
# guess3  -0.00442  0.00065  5.0e-04
# guess4  -0.00316  0.00125  7.4e-04
# guess5  -0.00159  0.00244  2.5e-03
# guess6  -0.00255  0.00202  3.2e-03
# guess7  -0.00378  0.00054  2.5e-05
# guess8  -0.01502 -0.01228 -1.0e-02
# guess9  -0.00400  0.00073  1.3e-03
# guess10 -0.00290  0.00026  5.1e-04
# guess11 -0.00361  0.00015  5.6e-04
# guess12 -0.00363 -0.00062 -1.3e-03
# guess13 -0.00247  0.00076 -2.0e-04
# guess14 -0.00299  0.00064  7.5e-04
# guess15 -0.00087  0.00193  1.9e-03
# guess16 -0.00453  0.00058  1.2e-03
# guess17 -0.00673 -0.00292 -1.3e-03
# guess18 -0.00156  0.00230  2.1e-03
# guess19 -0.00212  0.00335  3.5e-03
# guess20 -0.00407 -0.00012 -1.0e-03
# guess21 -0.00429  0.00023  7.3e-04
# guess22 -0.00045  0.00330  3.8e-04
# guess23 -0.02417 -0.02161 -2.2e-02
# guess24 -0.02597 -0.02335 -2.4e-02
# guess25 -0.00456  0.00013  5.7e-04
# guess26 -0.00357  0.00164  2.3e-03
# guess27 -0.00329  0.00042 -2.5e-04
# guess28 -0.00188  0.00230  2.7e-03 
est_theta$quantiles 
cdm_formatted[29:56, ] - est_theta$quantile[36:63, ]
#             5%      50%    97.5%
# slip1   0.00312  0.00514  0.00336
# slip2   0.00301  0.00475  0.00263
# slip3  -0.00496 -0.00050 -0.00044
# slip4  -0.00444 -0.00148 -0.00228
# slip5  -0.00283 -0.00039 -0.00125
# slip6  -0.00357 -0.00089 -0.00096
# slip7  -0.00385 -0.00066 -0.00121
# slip8   0.00154  0.00310  0.00157
# slip9  -0.00527 -0.00169 -0.00157
# slip10 -0.00447 -0.00063 -0.00177
# slip11 -0.00414 -0.00099 -0.00115
# slip12 -0.00485 -0.00078 -0.00087
# slip13 -0.00413 -0.00030 -0.00097
# slip14 -0.00490 -0.00037 -0.00082
# slip15 -0.00322 -0.00088 -0.00134

cdm_est_profile_probs <- cbind(
    p$probs$est - 1.96 * p$probs$se, p$probs$est,
    p$probs$est + 1.96 * p$probs$se
)
rownames(cdm_est_profile_probs) <- rownames(ecpe$attribute.patt)
cdm_est_profile_probs

ggdmc_est <- est_theta$quantile[29:35, ]
ggdmc_est

# Full comparison with output
result <- compare_profile_probs(ggdmc_est, cdm_est_profile_probs)
diff <- quick_compare(ggdmc_est, cdm_est_profile_probs)
diff
#          5%      50%    97.5%
# 000  0.0194  0.01684  0.01812
# 100  0.0115  0.00095  0.00088
# 010 -0.0137 -0.02459 -0.02019
# 110  0.0051 -0.00074 -0.00117
# 001  0.0113  0.01291  0.02149
# 101  0.0097  0.01473  0.02748
# 011 -0.0087 -0.00971 -0.00451