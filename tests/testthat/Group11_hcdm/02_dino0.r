# q(save = "no")

cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())
pkg <- c("ggdmc")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dino0.rda")
save_path <- file.path(home_dir, "Group11_hcdm/data/dino0.rda")
load(data_path)
# load(save_path)

# Fit random-effect models -------
fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.01, is_pblocked = TRUE,
    thin = 1L, pop_debug = F, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    sub_migration_prob = 0.00,
    is_pblocked = TRUE,
    thin = 1L, seed = 9032
)
save(fits0, fits1, file = save_path)

fits2 <- RestartSampling(fits1,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    is_pblocked = FALSE,
    is_hblocked = TRUE,
    thin = 1L, seed = 9032
)
save(fits0, fits1, fits2, file = save_path)


# fits3 <- RestartSampling(fits2,
#     pop_migration_prob = 0.01,
#     sub_migration_prob = 0.00,
#     is_pblocked = FALSE,
#     is_hblocked = TRUE,
#     thin = 1L, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, file = save_path)


fits <- fits2
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

options(digits = 2)
est_phi <- compare(phi, ps = true_vector)
#                 loc_guess1 loc_guess2 loc_guess3 loc_slip1 loc_slip2 loc_slip3
# True               0.1000      0.200      0.300    0.0100    0.0200    0.0300
# 5% Estimate        0.0055      0.017      0.023    0.0038    0.0026    0.0012
# 50% Estimate       0.0689      0.163      0.216    0.0487    0.0356    0.0162
# 97.5% Estimate     1.1982      2.083      0.814    2.0372    0.9558    1.3115
# Median-True       -0.0311     -0.037     -0.084    0.0387    0.0156   -0.0138
#                sca_guess1 sca_guess2 sca_guess3 sca_slip1 sca_slip2 sca_slip3
# True                 0.01       0.01       0.01      0.05      0.01      0.01
# 5% Estimate          0.40       0.30       0.36      0.38      0.28      0.11
# 50% Estimate         0.50       0.42       0.50      0.48      0.36      0.14
# 97.5% Estimate       2.45       5.50       2.30      3.76      2.67      2.65
# Median-True          0.49       0.41       0.49      0.43      0.35      0.13

#       loc_guess1 loc_guess2 loc_guess3 loc_slip1 loc_slip2 loc_slip3
# True               0.1000      0.200      0.300    0.0100    0.0200    0.0300
# 5% Estimate        0.0055      0.017      0.023    0.0038    0.0026    0.0012
# 50% Estimate       0.0689      0.163      0.216    0.0487    0.0356    0.0162
# 97.5% Estimate     1.1982      2.083      0.814    2.0372    0.9558    1.3115
# Median-True       -0.0311     -0.037     -0.084    0.0387    0.0156   -0.0138
#                sca_guess1 sca_guess2 sca_guess3 sca_slip1 sca_slip2 sca_slip3
# True                 0.01       0.01       0.01      0.05      0.01      0.01
# 5% Estimate          0.40       0.30       0.36      0.38      0.28      0.11
# 50% Estimate         0.50       0.42       0.50      0.48      0.36      0.14
# 97.5% Estimate       2.45       5.50       2.30      3.76      2.67      2.65
# Median-True          0.49       0.41       0.49      0.43      0.35      0.13


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
