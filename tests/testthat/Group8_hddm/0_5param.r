# q(save = "no")
cat("\n\n-------------------- 5 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc")

suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group5_gen_ddm/data/ddm0.rda")
gplot_path <- file.path(home_dir, "Group6_ddm_subjects/gplot.r")
save_path <- file.path(home_dir, "Group8_hddm/data/5param.rda")
load(data_path)

fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.00,
    nmc = 1000,
    is_pblocked = TRUE,
    thin = 2L, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    pop_migration_prob = 0.01,
    sub_migration_prob = 0.01,
    thin = 2L, seed = 9032
)
save(fits0, fits1, file = save_path)

fits2 <- RestartSampling(fits1,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.01,
    nmc = 500,
    thin = 2L, seed = 9032
)

save(fits0, fits1, fits2, file = save_path)


fits3 <- RestartSampling(fits2,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    is_hblocked = TRUE,
    thin = 2L, seed = 9032
)

save(fits0, fits1, fits2, fits3, file = save_path)

# load(save_path)
fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

options(digits = 2)

hat <- gelman(phi)
cat("mpsrf = ", hat$mpsrf, "\n")
hats <- lapply(thetas, gelman)
sort(unlist(lapply(hats, function(x) x$mpsrf)), decreasing = TRUE)


# est_phi <- summary(phi)
est_phi <- compare(phi, ps = true_vector)
#    loc_a loc_sz  loc_t0  loc_v  loc_z  sca_a sca_sz   sca_t0
# True            1.0000  0.250  0.1500  2.500  0.380 0.0500  0.010  0.02000
# 5% Estimate     0.9654  0.012  0.1394  2.295  0.355 0.0359  0.031  0.01565
# 50% Estimate    0.9937  0.104  0.1462  2.486  0.369 0.0535  0.113  0.01951
# 97.5% Estimate  1.0232  0.308  0.1549  2.708  0.392 0.8863  0.302  0.04736
# Median-True    -0.0063 -0.146 -0.0038 -0.014 -0.011 0.0035  0.103 -0.00049
#                sca_v  sca_z
# True           0.500 0.0100
# 5% Estimate    0.425 0.0052
# 50% Estimate   0.548 0.0169
# 97.5% Estimate 0.785 0.4304
# Median-True    0.048 0.0069

#                 loc_a loc_sz  loc_t0 loc_v loc_z sca_a sca_sz  sca_t0 sca_v
# True           1.0000  0.250  0.1500 2.500 0.380 0.050  0.010 0.02000 0.500
#    5 Estimate  0.9617  0.015  0.1404 2.348 0.378 0.048  0.048 0.01604 0.406
#   50 Estimate  0.9923  0.143  0.1487 2.532 0.393 0.066  0.137 0.02021 0.528
# 97.5 Estimate  1.0269  0.337  0.1575 2.752 0.420 0.695  0.314 0.12006 0.784
# Median-True   -0.0077 -0.107 -0.0013 0.032 0.013 0.016  0.127 0.00021 0.028
#                sca_z
# True          0.0100
#    5 Estimate 0.0043
#   50 Estimate 0.0158
# 97.5 Estimate 0.4071
# Median-True   0.0058



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
