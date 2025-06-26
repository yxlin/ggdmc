q(save = "no")
cat("\n\n-------------------- 5 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/"
gplot_path <- "~/Documents/ggdmc/tests/testthat/Group5/"

fn_path <- paste0(wkdir, "data/ddm_data0.rda")
helper_path <- paste0(gplot_path, "gplot.r")
save_path <- paste0(wkdir, "fit_data/5param.rda")

load(fn_path)
source(helper_path)
load(save_path)
# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.05,
#     thin = 8L, pop_debug = F, seed = 9032
# )
# save(fits0, file = save_path)


# fits1 <- RestartSampling(fits0,
#     pop_migration_prob = 0.02,
#     sub_migration_prob = 0.00,
#     thin = 4L, seed = 9032
# )
# save(fits0, fits1, file = save_path)

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


fits <- fits1
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

est_phi <- summary(phi)
est_phi <- compare(phi, ps = true_vector)

post_summary <- summary(phi, start = 501) # Discard first 500 as burn-in
# Custom quantiles
detailed_summary <- summary(phi, probability = seq(0.1, 0.9, by = 0.1))
detailed_summary$quantiles


# Extract specific elements
posterior_means <- post_summary$statistics[, "Mean"]
credible_intervals <- post_summary$quantiles[, c("5%", "97.5%")]
result <- summary_many(thetas)
result <- summary_many(thetas, verbose = T)


options(digits = 2)
result <- compare_many(thetas, ps = ps)
result <- compare_many(thetas, ps = ps, verbose = TRUE)

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
hat <- gelman(phi)
hat$psrf


hat <- gelman(fits[[1]]$phi)
hat <- gelman(fits[[2]]$phi)
hat <- gelman(fits[[3]]$phi)
hat

DT <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
DT <- ggdmc::prepare_thetas_data(thetas, start = 5000)

p1 <- plot_thetas(DT)
p1 <- plot_thetas(DT, start = 300, end = 400)
p1 <- plot_thetas(DT, start = 300, end = 400, subjects = 5)
p1 <- plot_thetas(DT, start = 300, end = 400, subjects = as.character(1:10))
p1 <- plot_thetas(DT, start = 300, end = 400, max_subjects = 8)

p2 <- gplot_thetas(DT)

p3 <- plot(fits[[1]]$phi, facet_chains = F)
p3 <- plot(fits[[2]]$phi, facet_chains = F)
p3 <- plot(fits[[3]]$phi, facet_chains = F)

p3 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
