#  q(save = "no")
cat("\n\n-------------------- DDM z model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
fn <- "~/Documents/ggdmc/tests/testthat/Group6/data/ddm_data7.rda"
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/"
save_path <- paste0(wkdir, "fit_data/7_z_model.rda")

load(fn)
load(save_path)
# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.05,
#     thin = 4L, pop_debug = F, seed = 9032
# )
# save(fits0, file = save_path)

# fits1 <- RestartSampling(fits0,
#     pop_migration_prob = 0.02,
#     sub_migration_prob = 0.00,
#     thin = 4L, seed = 9032
# )

# save(fits0, fits1, file = save_path)

# fits2 <- RestartSampling(fits1,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.00,
#     thin = 2L, seed = 9032
# )
# save(fits0, fits1, fits2, file = save_path)

fits <- fits2
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)
est_phi <- compare(phi, ps = true_vector)
#  loc_a loc_sv loc_sz loc_t0 loc_v.L loc_v.R loc_z.20 loc_z.80
# True           1.500  0.100  0.250 0.1500   2.350    2.25     0.38    0.680
# 5% Estimate    1.475  0.014  0.029 0.1439   2.112    1.73     0.39    0.679
# 50% Estimate   1.533  0.140  0.228 0.1518   2.285    2.12     0.42    0.731
# 97.5% Estimate 1.601  0.392  0.503 0.1610   2.538    2.46     0.46    0.780
# Median-True    0.033  0.040 -0.022 0.0018  -0.065   -0.13     0.04    0.051
#                sca_a sca_sv sca_sz   sca_t0 sca_v.L sca_v.R sca_z.20 sca_z.80
# True           0.050   0.01   0.01  0.02000    0.50    0.50    0.020    0.010
# 5% Estimate    0.081   0.14   0.13  0.01481    0.34    0.64    0.041    0.049
# 50% Estimate   0.116   0.32   0.26  0.01987    0.49    0.84    0.065    0.078
# 97.5% Estimate 0.200   0.67   0.49  0.02848    0.82    1.55    0.107    0.135
# Median-True    0.066   0.31   0.25 -0.00013   -0.01    0.34    0.045    0.068
rhat <- gelman(phi)

plot(phi, pll = F, den = T)

est_theta <- compare_many(thetas, ps = ps)
#    a      sv     sz      t0  v.L   v.R   z.20   z.80
# Mean  1.530  0.2793  0.329  0.1527 2.29  2.17  0.418  0.730
# True  1.484  0.0997  0.252  0.1517 2.30  2.29  0.376  0.684
# Diff -0.047 -0.1796 -0.077 -0.0010 0.01  0.12 -0.042 -0.046
# Sd    0.083  0.0693  0.099  0.0164 0.39  0.70  0.045  0.043
# True  0.045  0.0094  0.012  0.0185 0.52  0.40  0.021  0.011
# Diff -0.038 -0.0599 -0.087  0.0021 0.12 -0.30 -0.024 -0.033
# est_theta

# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     sub_migration_prob = 0.02,
#     thin = 2, seed = 9032
# )
# fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 8, seed = 9032)

# fits <- fits1
# fit <- RebuildPosterior(fits)
# options(digits = 2)
# est_theta <- ggdmc::compare(fit, ps = p_vector)
#                  a   sv    sz     t0   v.L  v.R z.20 z.80
# True           1.50 0.10 0.250 0.1500 2.350 2.25 0.38 0.68
# 5% Estimate    1.46 0.13 0.035 0.1310 1.779 1.81 0.35 0.73
# 50% Estimate   1.67 0.97 0.331 0.1563 2.389 2.51 0.45 0.88
# 97.5% Estimate 2.44 3.01 0.946 0.1733 4.288 4.60 0.61 1.25
# Median-True    0.17 0.87 0.081 0.0063 0.039 0.26 0.07 0.20
# plot(fit, den = T, pll = F)
