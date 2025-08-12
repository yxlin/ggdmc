# q(save = "no")
cat("\n\n-------------------- 6 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/"
fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ddm_data1.rda"
load(fn)
save_path <- paste0(wkdir, "data/6_ddm_param.rda")
load(save_path)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.04,
    thin = 4, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.05, thin = 2, seed = 9032)
fits2 <- RestartSampling_subject(fits1, sub_migration_prob = 0.05, thin = 16, seed = 9032)
save(fits0, fits1, file = save_path)

fits <- fits2
fit <- RebuildPosterior(fits)
options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
#                     a    sv    sz       t0     v      z
# True            1.000 0.100 0.250  0.15000 2.500  0.380
# 5% Estimate     0.909 0.033 0.035  0.14058 2.308  0.338
# 50% Estimate    0.959 0.316 0.269  0.14919 2.597  0.377
# 97.5% Estimate  1.038 1.140 0.548  0.15846 3.048  0.432
# Median-True    -0.041 0.216 0.019 -0.00081 0.097 -0.003


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

# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, start = fits[[1]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F, start = fits[[2]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F, start = fits[[3]]@nmc * 0.5)


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
# gplot_path <- "~/Documents/ggdmc/tests/testthat/Group5/"
# helper_path <- paste0(gplot_path, "gplot.r")
# source(helper_path)

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
