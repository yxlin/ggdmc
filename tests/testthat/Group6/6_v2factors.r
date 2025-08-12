#  q(save = "no")
cat("\n\n-------------------- DDM v x (S & NOISE) model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
fn <- "~/Documents/ggdmc/tests/testthat/Group6/data/ddm_data6.rda"
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/"
save_path <- paste0(wkdir, "fit_data/6_v2factors.rda")
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

# fits3 <- RestartSampling(fits2,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.01,
#     nmc = 1000L, report_length = 200L,
#     thin = 4L, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, file = save_path)

fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)
est_phi <- compare(phi, ps = true_vector)
est_theta <- compare_many(thetas, ps = ps)
#           a      sz        t0 v.s1.high v.s1.low v.s1.moderate v.s2.high
# Mean  1.0083  0.1490  0.15033     1.281    2.501       1.98483     1.690
# True  1.0057  0.2504  0.15262     1.348    2.579       1.98534     1.727
# Diff -0.0026  0.1014  0.00229     0.067    0.078       0.00051     0.037
# Sd    0.0486  0.0295  0.02178     0.437    0.195       0.08929     0.191
# True  0.0591  0.0094  0.02175     0.475    0.276       0.18426     0.369
# Diff  0.0105 -0.0202 -0.00003     0.038    0.081       0.09496     0.178
#      v.s2.low v.s2.moderate       z
# Mean   2.6593         2.083  0.3812
# True   2.6603         2.188  0.3811
# Diff   0.0011         0.105 -0.0001
# Sd     0.1747         0.076  0.0113
# True   0.0893         0.160  0.0087
# Diff  -0.0854         0.084 -0.0026
#                loc_a loc_sz loc_t0 loc_v.s1.high loc_v.s1.low loc_v.s1.moderate
# True           1.000  0.250 0.1500          1.50        2.500            2.0000
# 5% Estimate    0.985  0.014 0.1433          0.47        1.758            0.9187
# 50% Estimate   1.012  0.139 0.1519          1.24        2.516            2.0027
# 97.5% Estimate 1.050  0.338 0.1626          1.57        3.043            3.0355
# Median-True    0.012 -0.111 0.0019         -0.26        0.016            0.0027
#                loc_v.s2.high loc_v.s2.low loc_v.s2.moderate loc_z sca_a sca_sz
# True                  1.7000        2.700              2.20 0.380 0.050  0.010
# 5% Estimate           1.4854        2.454              1.84 0.373 0.044  0.069
# 50% Estimate          1.6943        2.723              2.06 0.386 0.061  0.171
# 97.5% Estimate        1.9224        9.530              2.26 0.403 0.095  0.326
# Median-True          -0.0057        0.023             -0.14 0.006 0.011  0.161
#                sca_t0 sca_v.s1.high sca_v.s1.low sca_v.s1.moderate
# True           0.0200          0.50         0.25              0.20
# 5% Estimate    0.0179          0.39         0.23              0.16
# 50% Estimate   0.0228          0.63         0.45              0.36
# 97.5% Estimate 0.0317          1.40         4.56              9.68
# Median-True    0.0028          0.13         0.20              0.16
#                sca_v.s2.high sca_v.s2.low sca_v.s2.moderate sca_z
# True                   0.400         0.10              0.15 0.010
# 5% Estimate            0.209         0.24              0.12 0.012
# 50% Estimate           0.369         0.45              0.27 0.022
# 97.5% Estimate         0.749         8.55              0.50 0.042
# Median-True           -0.031         0.35              0.12 0.012
rhat <- gelman(phi)

# rhat <- gelman(fits[[1]]$phi)
# rhat <- gelman(fits[[2]]$phi)
# rhat <- gelman(fits[[3]]$phi)

# plot(fits[[1]]$phi)
# plot(fits[[2]]$phi)
# plot(fits[[3]]$phi)

# plot(fits[[1]]$phi, den = T, pll = F)
# plot(fits[[2]]$phi, den = T, pll = F)
# plot(fits[[3]]$phi, den = T, pll = F)

plot(phi, den = T, pll = F)
