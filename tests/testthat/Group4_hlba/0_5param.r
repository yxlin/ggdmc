# q(save = "no")
cat("\n\n---------------- 5 parameters HB-LBA ---------------")
rm(list = ls())
pkg <- c("ggdmc")

sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
wkdir <- file.path(home_dir, "Group4_hlba/")

data_path <- file.path(home_dir, "Group0_gen_data/data/lba_data0.rda")
save_path <- file.path(home_dir, "Group4_hlba/data/5param.rda")

cat("\nWorking directory: ", wkdir, "\n")
load(data_path)


# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.01, is_pblocked = TRUE,
#     thin = 1L, pop_debug = F, seed = 9032
# )

fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.01, is_pblocked = TRUE,
    thin = 1L, pop_debug = F, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    pop_migration_prob = 0.01,
    sub_migration_prob = 0.00,
    is_pblocked = FALSE,
    is_hblocked = TRUE,
    thin = 1L, seed = 9032
)
save(fits0, fits1, file = save_path)

fits2 <- RestartSampling(fits1,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    is_pblocked = FALSE,
    is_hblocked = FALSE,
    thin = 2L, seed = 9032
)

save(fits0, fits1, fits2, file = save_path)


fits <- fits2
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)
hat <- gelman(phi)
cat("mpsrf = ", hat$mpsrf, "\n")


options(digits = 2)
est_phi <- compare(phi, ps = true_vector, start = phi@nmc * 0.5)
#                 loc_A loc_B loc_mean_v.false loc_mean_v.true loc_t0 sca_A sca_B
# True            0.400  0.50            0.150            2.50  0.300  0.10 0.100
# 5% Estimate     0.018  0.49            0.066            2.50  0.196  0.13 0.054
# 50% Estimate    0.211  0.68            0.444            2.68  0.284  0.29 0.116
# 97.5% Estimate  1.006  1.01            0.818            3.01  3.495  4.01 1.541
# Median-True    -0.189  0.18            0.294            0.18 -0.016  0.19 0.016
#                sca_mean_v.false sca_mean_v.true sca_t0
# True                      0.200           0.200  0.050
# 5% Estimate               0.093           0.118  0.027
# 50% Estimate              0.263           0.229  0.038
# 97.5% Estimate            1.525           1.585  4.390
# Median-True               0.063           0.029 -0.012
# mpsrf =  1

#                loc_A loc_B loc_mean_v.false loc_mean_v.true loc_t0 sca_A sca_B
# True           0.400  0.50            0.150            2.50   0.30 0.100 0.100
#    5 Estimate  0.013  0.46            0.027            2.49   0.12 0.089 0.101
#   50 Estimate  0.150  0.66            0.264            2.64   0.28 0.215 0.165
# 97.5 Estimate  0.678  1.84            0.770            3.01   1.46 3.147 5.062
# Median-True   -0.250  0.16            0.114            0.14  -0.02 0.115 0.065
#               sca_mean_v.false sca_mean_v.true  sca_t0
# True                      0.20           0.200  0.0500
#    5 Estimate             0.26           0.103  0.0279
#   50 Estimate             0.47           0.223  0.0403
# 97.5 Estimate             1.85           1.669  3.3192
# Median-True               0.27           0.023 -0.0097




plot(phi, pll = F, den = T)
DT <- prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
p1 <- plot_thetas(DT)
p2 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
p2 <- plot(fits[[2]]$phi, facet_chains = F, start = fits[[2]]$phi@nmc * 0.5)
p2 <- plot(fits[[3]]$phi, facet_chains = F, start = fits[[3]]$phi@nmc * 0.5)
p2 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
