# q(save = "no")
cat("\n\n-------------------- 6 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group5_gen_ddm/data/ddm_data2.rda")
gplot_path <- file.path(home_dir, "Group6_ddm_subjects/gplot.r")
save_path <- file.path(home_dir, "Group8_hddm/data/v6param.rda")
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

load(save_path)
fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

options(digits = 3)
hat <- gelman(phi)
cat("mpsrf = ", hat$mpsrf, "\n")
hats <- lapply(thetas, gelman)
sort(unlist(lapply(hats, function(x) x$mpsrf)), decreasing = TRUE)

est_phi <- compare(phi, ps = true_vector)
#                 loc_a loc_sv loc_sz  loc_t0 loc_v  loc_z   sca_a sca_sv sca_sz
# True           1.0000  0.100  0.250  0.1500  2.50 0.3800  0.0500   0.01  0.010
# 5% Estimate    0.9725  0.028  0.016  0.1362  2.40 0.3642  0.0250   0.10  0.038
# 50% Estimate   1.0014  0.275  0.162  0.1457  2.61 0.3825  0.0452   0.33  0.130
# 97.5% Estimate 1.0610  0.827  0.410  0.1561  2.94 0.4228  0.5637   1.21  0.318
# Median-True    0.0014  0.175 -0.088 -0.0043  0.11 0.0025 -0.0048   0.32  0.120
#                sca_t0  sca_v sca_z
# True            0.020  0.500 0.010
# 5% Estimate     0.018  0.333 0.013
# 50% Estimate    0.023  0.459 0.028
# 97.5% Estimate  0.118  1.547 0.444
# Median-True     0.003 -0.041 0.018

#                 loc_a  loc_sz  loc_t0 loc_v.s1 loc_v.s2  loc_z     sca_a sca_sz
# True           1.0000  0.2500 0.15000   2.5000    2.100 0.3800  0.050000 0.0100
# 5% Estimate    0.9873  0.0225 0.14932   2.3378    2.075 0.3746  0.023062 0.0469
# 50% Estimate   1.0189  0.2157 0.15905   2.5699    2.228 0.3983  0.049063 0.1548
# 97.5% Estimate 1.0757  0.4808 0.17237   2.8555    2.419 0.4491  0.485693 0.3668
# Median-True    0.0189 -0.0343 0.00905   0.0699    0.128 0.0183 -0.000937 0.1448
#                 sca_t0 sca_v.s1 sca_v.s2   sca_z
# True           0.02000    0.500   0.2500 0.01000
# 5% Estimate    0.01726    0.352   0.0551 0.00798
# 50% Estimate   0.02205    0.542   0.1691 0.02468
# 97.5% Estimate 0.20416    0.918   0.6193 0.44796
# Median-True    0.00205    0.042  -0.0809 0.01468

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
# p4 <- plot(phi, den = F, pll = F, start = phi@nmc * 0.85, hide_legend = FALSE)

# p4 <- plot(fits[[2]]$phi, den = F, pll = F, start = fits[[2]]$phi@nmc * 0.5, hide_legend = TRUE)
