#  q(save = "no")
cat("\n\n-------------------- t0 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
fn <- "~/Documents/ggdmc/tests/testthat/Group6/data/ddm_data5.rda"
load(fn)

wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/"
save_path <- paste0(wkdir, "fit_data/5_t0_model.rda")
load(save_path)

# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.05,
#     thin = 8L, pop_debug = F, seed = 9032
# )
# save(fits0, file = save_path)


# fits1 <- RestartSampling(fits0,
#     pop_migration_prob = 0.02,
#     sub_migration_prob = 0.02,
#     thin = 2L, seed = 9032
# )
# save(fits0, fits1, file = save_path)

# fits2 <- RestartSampling(fits1,
#     pop_migration_prob = 0.01,
#     sub_migration_prob = 0.00,
#     thin = 2L, seed = 9032
# )
# save(fits0, fits1, fits2, file = save_path)

# fits3 <- RestartSampling(fits2,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.00,
#     thin = 2L, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, file = save_path)

fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)
est_phi <- compare(phi, ps = true_vector)
#                  loc_a loc_sz loc_t0.excellent loc_t0.good loc_t0.poor  loc_v
# True            1.0000  0.250           0.0500      0.1000      0.2000  2.500
# 5% Estimate     0.9727  0.023           0.0484      0.0891      0.1901  2.265
# 50% Estimate    0.9973  0.176           0.0523      0.0977      0.2017  2.427
# 97.5% Estimate  1.0273  0.344           0.0585      0.1068      0.2147  2.662
# Median-True    -0.0027 -0.074           0.0023     -0.0023      0.0017 -0.073
#                 loc_z sca_a sca_sz sca_t0.excellent sca_t0.good sca_t0.poor
# True           0.3800 0.050  0.010           0.0100       0.020      0.0300
# 5% Estimate    0.3695 0.053  0.029           0.0060       0.018      0.0279
# 50% Estimate   0.3813 0.069  0.107           0.0087       0.023      0.0344
# 97.5% Estimate 0.3962 0.098  0.257           0.0131       0.031      0.0453
# Median-True    0.0013 0.019  0.097          -0.0013       0.003      0.0044
#                sca_v  sca_z
# True           0.500 0.0100
# 5% Estimate    0.433 0.0051
# 50% Estimate   0.551 0.0162
# 97.5% Estimate 0.758 0.0325
# Median-True    0.051 0.0062
rhat <- gelman(phi)
rhat

est_theta <- compare_many(thetas, ps = ps)
#           a     sz t0.excellent  t0.good t0.poor      v       z
# Mean  0.9973  0.211       0.0527 9.8e-02  0.2021  2.437  0.3817
# True  0.9936  0.248       0.0539 9.9e-02  0.2032  2.526  0.3806
# Diff -0.0037  0.037       0.0012 1.1e-03  0.0012  0.088 -0.0011
# Sd    0.0588  0.039       0.0066 2.1e-02  0.0326  0.489  0.0076
# True  0.0535  0.009       0.0095 2.1e-02  0.0312  0.374  0.0091
# Diff -0.0053 -0.030       0.0029 5.2e-05 -0.0015 -0.115  0.0015
plot(phi, pll = F, den = T)
