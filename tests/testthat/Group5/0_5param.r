#    q(save = "no")
cat("\n\n-------------------- 5 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc")

suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

wkdir <- "~/Documents/ggdmc/tests/testthat/Group5/"
data_dir <- "~/Documents/ggdmc/tests/testthat/Group5/data/"
data_path <- paste0(data_dir, "ddm_data0.rda")
save_path <- paste0(wkdir, "fit_data/5param.rda")

options(digits = 2)
load(data_path)
load(save_path)

# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.00,
#     nmc = 1000,
#     is_pblocked = TRUE,
#     thin = 2L, seed = 9032
# )
# save(fits0, file = save_path)

# fits1 <- RestartSampling(fits0,
#     pop_migration_prob = 0.01,
#     sub_migration_prob = 0.01,
#     thin = 2L, seed = 9032
# )
# save(fits0, fits1, file = save_path)

# fits2 <- RestartSampling(fits1,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.01,
#     nmc = 500,
#     thin = 2L, seed = 9032
# )

# save(fits0, fits1, fits2, file = save_path)


fits3 <- RestartSampling(fits2,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    is_hblocked = TRUE,
    thin = 2L, seed = 9032
)

save(fits0, fits1, fits2, fits3, file = save_path)
fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

hat <- gelman(phi)
cat("mpsrf = ", hat$mpsrf, "\n")
hats <- lapply(thetas, gelman)
sort(unlist(lapply(hats, function(x) x$mpsrf)), decreasing = TRUE)

# DT1 <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.85)
# DT2 <- ggdmc::prepare_thetas_data(fits[[2]]$subject_theta, start = fits[[2]]$phi@nmc * 0.85)
# DT3 <- ggdmc::prepare_thetas_data(fits[[3]]$subject_theta, start = fits[[3]]$phi@nmc * 0.85)
# p1 <- plot_thetas(DT1)
# p1 <- plot_thetas(DT2)
# p1 <- plot_thetas(DT3)
# p2 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
# p2 <- plot(fits[[2]]$phi, facet_chains = F, start = fits[[2]]$phi@nmc * 0.5)
# p2 <- plot(fits[[3]]$phi, facet_chains = F, start = fits[[3]]$phi@nmc * 0.5)

# p2 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
# p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)


# est_phi <- summary(phi)
est_phi <- compare(phi, ps = true_vector)
#                 loc_a loc_sz  loc_t0 loc_v loc_z sca_a sca_sz  sca_t0 sca_v
# True           1.0000  0.250  0.1500  2.50 0.380 0.050   0.01  0.0200 0.500
# 5% Estimate    0.9761  0.023  0.1382  2.45 0.386 0.053   0.05  0.0137 0.454
# 50% Estimate   1.0043  0.216  0.1453  2.65 0.403 0.072   0.16  0.0175 0.593
# 97.5% Estimate 1.0390  0.418  0.1536  2.89 0.424 0.110   0.34  0.0264 0.862
# Median-True    0.0043 -0.034 -0.0047  0.15 0.023 0.022   0.15 -0.0025 0.093
#                 sca_z
# True           0.0100
# 5% Estimate    0.0041
# 50% Estimate   0.0159
# 97.5% Estimate 0.0419
# Median-True    0.0059

# post_summary <- summary(phi, start = 501) # Discard first 500 as burn-in
# Custom quantiles
# detailed_summary <- summary(phi, probability = seq(0.1, 0.9, by = 0.1))
# detailed_summary$quantiles


# Extract specific elements
# posterior_means <- post_summary$statistics[, "Mean"]
# credible_intervals <- post_summary$quantiles[, c("5%", "97.5%")]
# result <- summary_many(thetas)
# result <- summary_many(thetas, verbose = T)


# options(digits = 2)
# result <- compare_many(thetas, ps = ps)
# result <- compare_many(thetas, ps = ps, verbose = TRUE)

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

#                 loc_a loc_sz  loc_t0 loc_v loc_z sca_a sca_sz  sca_t0 sca_v
# True           1.0000  0.250  0.1500  2.50 0.380 0.050  0.010  0.0200  0.50
# 5% Estimate    0.9614  0.027  0.1340  0.23 0.227 0.054  0.079  0.0138  0.46
# 50% Estimate   1.0049  0.223  0.1453  2.60 0.405 0.076  0.203  0.0181  0.64
# 97.5% Estimate 1.2946  0.617  0.1711  2.99 0.464 1.017  1.160  0.1174  4.06
# Median-True    0.0049 -0.027 -0.0047  0.10 0.025 0.026  0.193 -0.0019  0.14
#                 sca_z
# True           0.0100
# 5% Estimate    0.0083
# 50% Estimate   0.0208
# 97.5% Estimate 0.6235
# Median-True    0.0108




# DT <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
# DT <- ggdmc::prepare_thetas_data(thetas, start = 5000)

# p1 <- plot_thetas(DT)
# p1 <- plot_thetas(DT, start = 300, end = 400)
# p1 <- plot_thetas(DT, start = 300, end = 400, subjects = 5)
# p1 <- plot_thetas(DT, start = 300, end = 400, subjects = as.character(1:10))
# p1 <- plot_thetas(DT, start = 300, end = 400, max_subjects = 8)

# p2 <- gplot_thetas(DT)

# p3 <- plot(fits[[1]]$phi, facet_chains = F)
# p3 <- plot(fits[[2]]$phi, facet_chains = F)
# p3 <- plot(fits[[3]]$phi, facet_chains = F)

# p3 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
# p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)

# Fit one subject -------------------------
# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     sub_migration_prob = 0.05,
#     thin = 2, seed = 9032
# )
# fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 2, seed = 9032)


# fits <- fits1
# fit <- RebuildPosterior(fits)

# hat <- gelman(fit)
# cat("mpsrf = ", hat$mpsrf, "\n")


# est_theta <- ggdmc::compare(fit, ps = p_vector)
#                    a    sz    t0     v     z
# True           1.000 0.250 0.150  2.50 0.380
# 5% Estimate    0.958 0.044 0.147  2.06 0.363
# 50% Estimate   1.025 0.394 0.161  2.37 0.415
# 97.5% Estimate 1.125 0.701 0.170  2.78 0.487
# Median-True    0.025 0.144 0.011 -0.13 0.035

# Fit only hyper parameters -------------------------
# fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
#     sub_migration_prob = 0.05, thin = 2
# )
# fits1 <- RestartSampling_hyper(fits0)

# fits <- fits1
# fit <- RebuildPosterior(fits)

# est_phi <- ggdmc::compare(fit, ps = true_vector)
#                 loc_a loc_sz  loc_t0 loc_v  loc_z sca_a sca_sz   sca_t0 sca_v
# True           1.0000 0.2500  0.1500  2.50 0.3800 0.050 0.0100  0.02000  0.50
# 5% Estimate    0.9841 0.2479  0.1386  2.53 0.3817 0.049 0.0089  0.01571  0.45
# 50% Estimate   1.0023 0.2511  0.1443  2.69 0.3844 0.060 0.0110  0.01906  0.55
# 97.5% Estimate 1.0235 0.2549  0.1510  2.89 0.3877 0.080 0.0144  0.02503  0.73
# Median-True    0.0023 0.0011 -0.0057  0.19 0.0044 0.010 0.0010 -0.00094  0.05
#                   sca_z
# True            0.01000
# 5% Estimate     0.00773
# 50% Estimate    0.00946
# 97.5% Estimate  0.01239
# Median-True    -0.00054
