# q(save = "no")
cat("\n\n-------------------- 6 parameters --------------------")

rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
wkdir <- file.path(home_dir, "Group4_hlba/")
data_path <- file.path(home_dir, "Group0_gen_data/data/lba_data1.rda")
save_path <- file.path(home_dir, "Group4_hlba/data/6param.rda")

cat("\nWorking directory: ", wkdir, "\n")
load(data_path)


fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.02,
    thin = 10L, seed = 9032
)
# save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    pop_migration_prob = 0.01,
    sub_migration_prob = 0.00,
    is_hblocked = T,
    thin = 4L, seed = 9032
)
# save(fits0, fits1, file = save_path)


fits <- fits1
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

options(digits = 2)
est_phi <- compare(phi, ps = true_vector)

# Processing time: 211.11 secs.
#                loc_A loc_B loc_mean_v.false loc_mean_v.true loc_sd_v.true
# True            0.40 0.500            0.150            2.50         0.100
# 5% Estimate     0.36 0.417            0.020            2.45         0.014
# 50% Estimate    0.41 0.592            0.183            2.62         0.086
# 97.5% Estimate  0.50 0.718            0.577            3.21         0.162
# Median-True     0.01 0.092            0.033            0.12        -0.014
#                loc_t0 sca_A sca_B sca_mean_v.false sca_mean_v.true
# True            0.300 0.100 0.100            0.200          0.2000
# 5% Estimate     0.251 0.091 0.138            0.100          0.0988
# 50% Estimate    0.283 0.117 0.195            0.253          0.2074
# 97.5% Estimate  0.352 0.228 0.382            0.542          1.3448
# Median-True    -0.017 0.017 0.095            0.053          0.0074
#                sca_sd_v.true sca_t0
# True                   0.100  0.050
# 5% Estimate            0.076  0.036
# 50% Estimate           0.115  0.061
# 97.5% Estimate         0.202  0.121
# Median-True            0.015  0.011

# DT1 <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
# DT2 <- ggdmc::prepare_thetas_data(fits[[2]]$subject_theta, start = fits[[2]]$phi@nmc * 0.5)
# DT3 <- ggdmc::prepare_thetas_data(fits[[3]]$subject_theta, start = fits[[3]]$phi@nmc * 0.5)
# p1 <- plot_thetas(DT1)
# p1 <- plot_thetas(DT2)
# p1 <- plot_thetas(DT3)

# p2 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
# p2 <- plot(fits[[2]]$phi, facet_chains = F, start = fits[[2]]$phi@nmc * 0.5)
# p2 <- plot(fits[[3]]$phi, facet_chains = F, start = fits[[3]]$phi@nmc * 0.5)

# p2 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
# p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
