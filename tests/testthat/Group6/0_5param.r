# q(save = "no")
cat("\n\n-------------------- 5 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
fn <- "~/Documents/ggdmc/tests/testthat/Group6/data/ddm_data0.rda"
load(fn)
gplot_path <- "~/Documents/ggdmc/tests/testthat/Group5/"
helper_path <- paste0(gplot_path, "gplot.r")
source(helper_path)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.02,
    thin = 2, seed = 9032
)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 1, seed = 9032)


fits <- fits1
fit <- RebuildPosterior(fits)
options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
#                    a    sz    t0     v     z
# True           1.000 0.250 0.150  2.50 0.380
# 5% Estimate    0.958 0.044 0.147  2.06 0.363
# 50% Estimate   1.025 0.394 0.161  2.37 0.415
# 97.5% Estimate 1.125 0.701 0.170  2.78 0.487
# Median-True    0.025 0.144 0.011 -0.13 0.035

# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F)
# p1 <- ggdmc::plot(fit, facet_chains = T, start = 20)

# p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = T, start = fit@nmc * 0.5)
# p0 <- ggdmc::plot(fit, pll = F, den = T)
# p0 <- ggdmc::plot(fit, pll = F, den = F)


# p1 <- ggdmc::plot(fits[[1]], facet_chains = F)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F)


# fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors, sub_migration_prob = 0.05, thin = 2)
# fits1 <- RestartSampling_hyper(fits0)

# fits <- fits1
# fit <- RebuildPosterior(fits)

# est_phi <- ggdmc::compare(fit, ps = true_vector)
#                  loc_a loc_sz loc_t0  loc_v    loc_z   sca_a sca_sz  sca_t0
# True            1.0000  0.250  0.150  2.500  0.38000  0.0500 0.0100  0.0200
# 5% Estimate     0.9837  0.246  0.147  2.344  0.37663  0.0395 0.0085  0.0155
# 50% Estimate    0.9976  0.249  0.152  2.483  0.37914  0.0482 0.0103  0.0187
# 97.5% Estimate  1.0143  0.253  0.159  2.670  0.38211  0.0637 0.0137  0.0247
# Median-True    -0.0024 -0.001  0.002 -0.017 -0.00086 -0.0018 0.0003 -0.0013
#                  sca_v   sca_z
# True            0.5000  0.0100
# 5% Estimate     0.4074  0.0071
# 50% Estimate    0.4967  0.0087
# 97.5% Estimate  0.6728  0.0113
# Median-True    -0.0033 -0.0013


# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.05,
#     thin = 4L, pop_debug = F, seed = 9032
# )

# # save(fits0, file = save_path)
# fits1 <- RestartSampling(fits0,
#     pop_migration_prob = 0.02,
#     sub_migration_prob = 0.00,
#     thin = 2L, seed = 9032
# )

# fits2 <- RestartSampling(fits1,
#     pop_migration_prob = 0.01,
#     sub_migration_prob = 0.00,
#     thin = 2L, seed = 9032
# )

# fits3 <- RestartSampling(fits2,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.00,
#     thin = 2L, seed = 9032
# )


# fits <- fits2
# phi <- RebuildHyper(fits)
# thetas <- RebuildPosteriors(fits)

# est_phi <- compare(phi, ps = true_vector)
#   loc_a loc_sz loc_t0  loc_v   loc_z  sca_a sca_sz sca_t0 sca_v
# True            1.0000  0.250 0.1500  2.500  0.3800  0.050  0.010 0.0200 0.500
# 5% Estimate     0.9779  0.014 0.1478  2.307  0.3633  0.013  0.063 0.0173 0.432
# 50% Estimate    0.9963  0.126 0.1544  2.486  0.3754  0.033  0.141 0.0219 0.552
# 97.5% Estimate  1.0175  0.359 0.1633  2.723  0.3886  0.061  0.290 0.0303 0.730
# Median-True    -0.0037 -0.124 0.0044 -0.014 -0.0046 -0.017  0.131 0.0019 0.052
#                 sca_z
# True           0.0100
# 5% Estimate    0.0060
# 50% Estimate   0.0136
# 97.5% Estimate 0.0263
# Median-True    0.0036


# DT <- ggdmcPhi::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
# p1 <- plot_thetas(DT)
# p2 <- gplot_thetas(DT)

# p3 <- plot(fits[[1]]$phi, facet_chains = F)
# p3 <- plot(fits[[2]]$phi, facet_chains = F)
# p3 <- plot(fits[[3]]$phi, facet_chains = F)


# p3 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
# p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
