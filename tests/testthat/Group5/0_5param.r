q(save = "no")
cat("\n\n-------------------- 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ggdmc/tests/testthat/Group5/"
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/data/"

data_path <- paste0(data_dir, "lba_data0.rda")
save_path <- paste0(wkdir, "data/5param.rda")
options(digits = 2)

load(data_path)
load(save_path)

# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.05,
#     thin = 6L, pop_debug = F, seed = 9032
# )
# save(fits0, file = save_path)

# fits1 <- RestartSampling(fits0,
#     pop_migration_prob = 0.01,
#     sub_migration_prob = 0.00,
#     thin = 4L, seed = 9032
# )
# save(fits0, fits1, file = save_path)

# fits2 <- RestartSampling(fits1,
#     pop_migration_prob = 0.01,
#     sub_migration_prob = 0.00,
#     thin = 6L, seed = 9032
# )
# save(fits0, fits1, fits2, file = save_path)


fits <- fits2
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

est_phi <- compare(phi, ps = true_vector)
#                 loc_A loc_B loc_mean_v.false loc_mean_v.true loc_t0 sca_A sca_B
# True            0.400 0.500            0.150            2.50  0.300 0.100 0.100
# 5% Estimate     0.022 0.078            0.025            1.50  0.036 0.098 0.065
# 50% Estimate    0.269 0.664            0.389            2.78  0.267 0.275 0.241
# 97.5% Estimate  1.147 1.204            1.536            4.49  0.594 3.343 2.666
# Median-True    -0.131 0.164            0.239            0.28 -0.033 0.175 0.141
#                sca_mean_v.false sca_mean_v.true sca_t0
# True                       0.20            0.20  0.050
# 5% Estimate                0.22            0.21  0.045
# 50% Estimate               0.56            0.50  0.068
# 97.5% Estimate             3.34            7.69  1.666
# Median-True                0.36            0.30  0.018

rhat <- gelman(phi)
plot(phi, pll = F, den = T)
est_theta <- compare_many(thetas, ps = ps)
#    A      B mean_v.false mean_v.true     t0
# Mean  0.359  0.608        0.373       2.572 0.2976
# True  0.409  0.505        0.288       2.489 0.3163
# Diff  0.050 -0.103       -0.085      -0.082 0.0187
# Sd    0.174  0.030        0.074       0.099 0.0306
# True  0.103  0.105        0.156       0.220 0.0345
# Diff -0.071  0.075        0.082       0.121 0.0039

# hat <- gelman(fits[[1]]$phi, verbose = TRUE)
# hat <- gelman(fits[[2]]$phi, verbose = TRUE)
# hat <- gelman(fits[[3]]$phi, verbose = TRUE)

# fits2 <- RestartSampling(fits1, pop_migration_prob = 0.00, sub_migration_prob = 0.05, thin = 2, seed = 9032)
# save(fits0, fits1, fits2, file = save_path)

# fits3 <- RestartSampling(fits2, pop_migration_prob = 0.01, sub_migration_prob = 0.00, thin = 8, seed = 9032)
# save(fits0, fits1, fits2, fits3, file = save_path)



#         loc_A loc_B loc_mean_v.false loc_mean_v.true loc_t0 sca_A
# True            0.400  0.50            0.150            2.50  0.300 0.100
# 5% Estimate     0.022  0.47            0.075            2.48  0.236 0.058
# 50% Estimate    0.198  0.72            0.453            2.68  0.269 0.152
# 97.5% Estimate  0.570  0.94            0.887            3.22  0.324 1.486
# Median-True    -0.202  0.22            0.303            0.18 -0.031 0.052
#                  sca_B sca_mean_v.false sca_mean_v.true sca_t0
# True           0.10000             0.20           0.200 0.0500
# 5% Estimate    0.03992             0.16           0.174 0.0452
# 50% Estimate   0.10045             0.32           0.271 0.0585
# 97.5% Estimate 1.04802             1.65           2.643 0.3208
# Median-True    0.00045             0.12           0.071 0.0085


# DT <- ggdmcPhi::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
# p1 <- plot_thetas(DT)
# p2 <- gplot_thetas(DT)

# p3 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
# p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
