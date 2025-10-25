# q(save = "no")
cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())

pkg <- c("ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina0.rda")
save_path <- file.path(home_dir, "Group11_hcdm/data/dina0.rda")


cat("\nLoad data: ", data_path, "\n")
load(data_path)
load(save_path)

# Fit random-effect models -------
# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.02, is_pblocked = TRUE,
#     thin = 2L, pop_debug = F, seed = 9032
# )

# save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    sub_migration_prob = 0.00,
    is_pblocked = TRUE,
    thin = 4L, seed = 9032
)

save(fits0, fits1, file = save_path)

# fits2 <- RestartSampling(fits1,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.00,
#     is_pblocked = FALSE,
#     is_hblocked = TRUE,
#     thin = 1L, seed = 9032
# )
# save(fits0, fits1, fits2, file = save_path)

# fits3 <- RestartSampling(fits2,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.02,
#     is_pblocked = FALSE,
#     is_hblocked = TRUE,
#     thin = 2L, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, file = save_path)


fits <- fits1
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

options(digits = 2)
est_phi <- compare(phi, ps = true_vector)

#                loc_guess1 loc_guess2 loc_guess3 loc_slip1 loc_slip2 loc_slip3
# True               0.1000      0.200      0.300    0.0100    0.0200    0.0300
# 5% Estimate        0.0059      0.017      0.085    0.0043    0.0031    0.0054
# 50% Estimate       0.0782      0.175      0.410    0.0548    0.0413    0.0677
# 97.5% Estimate     1.6198      1.302      0.934    1.5039    1.2014    1.6167
# Median-True       -0.0218     -0.025      0.110    0.0448    0.0213    0.0377
#                sca_guess1 sca_guess2 sca_guess3 sca_slip1 sca_slip2 sca_slip3
# True                 0.01       0.01       0.01      0.05      0.01      0.01
# 5% Estimate          0.38       0.29       0.14      0.35      0.25      0.44
# 50% Estimate         0.50       0.44       0.19      0.46      0.36      0.59
# 97.5% Estimate       3.39       3.06       2.72      2.98      3.30      3.87
# Median-True          0.49       0.43       0.18      0.41      0.35      0.58
# 3200 * 2
# 4400 * 2
# 3800 * 2

#                loc_guess1 loc_guess2 loc_guess3 loc_pi1 loc_pi2 loc_pi3
# True                0.10        0.2       0.30   0.100  0.3000  0.5000
#    5 Estimate       0.32        0.3       0.39   0.077  0.0018  0.0037
#   50 Estimate       3.26        3.2       3.48   0.169  0.0215  0.0416
# 97.5 Estimate       9.22        9.0       9.10   2.398  2.1566  2.3427
# Median-True         3.16        3.0       3.18   0.069 -0.2785 -0.4584
#               loc_slip1 loc_slip2 loc_slip3 sca_guess1 sca_guess2 sca_guess3
# True              0.020      0.05     0.070       0.01       0.01       0.01
#    5 Estimate     0.063      0.24     0.056       3.62       3.96       4.04
#   50 Estimate     0.515      0.50     0.666       7.50       7.72       7.78
# 97.5 Estimate     1.728      1.93     1.272       9.89       9.91       9.90
# Median-True       0.495      0.45     0.596       7.49       7.71       7.77
#               sca_pi1 sca_pi2 sca_pi3 sca_slip1 sca_slip2 sca_slip3
# True            0.010   0.030   0.050      0.02     0.050      0.07
#    5 Estimate   0.010   0.060   0.065      0.12     0.091      0.15
#   50 Estimate   0.029   0.087   0.115      0.18     0.132      0.26
# 97.5 Estimate   4.091   3.797   4.776      3.96     4.337      3.93
# Median-True     0.019   0.057   0.065      0.16     0.082      0.19

# loc_guess1 loc_guess2 loc_guess3 loc_pi1 loc_pi2 loc_pi3
# True                0.10       0.20       0.30   0.100  0.3000  0.5000
#    5 Estimate       0.24       0.22       0.27   0.160  0.0017  0.0049
#   50 Estimate       1.99       1.94       2.22   0.168  0.0198  0.0470
# 97.5 Estimate       6.29       5.41       5.57   0.188  0.0681  0.1060
# Median-True         1.89       1.74       1.92   0.068 -0.2802 -0.4530
#               loc_slip1 loc_slip2 loc_slip3 sca_guess1 sca_guess2 sca_guess3
# True               0.02      0.05      0.07       0.01       0.01       0.01
#    5 Estimate      0.46      0.46      0.53       1.53       1.52       1.59
#   50 Estimate      0.53      0.50      0.68       3.53       3.65       3.67
# 97.5 Estimate      0.60      0.56      0.75       8.97       8.99       8.95
# Median-True        0.51      0.45      0.61       3.52       3.64       3.66
#               sca_pi1 sca_pi2 sca_pi3 sca_slip1 sca_slip2 sca_slip3
# True           0.0100   0.030   0.050     0.020     0.050      0.07
#    5 Estimate  0.0046   0.046   0.051     0.099     0.086      0.13
#   50 Estimate  0.0120   0.078   0.092     0.144     0.123      0.19
# 97.5 Estimate  0.9072   0.135   0.170     1.272     0.202      0.99
# Median-True    0.0020   0.048   0.042     0.124     0.073      0.12


DT <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta,
    start = fits[[1]]$phi@nmc * 0.5
)
# DT <- ggdmc::prepare_thetas_data(thetas, start = 5000)

# pdf("Rplot.pdf")
# p1 <- plot_thetas(DT)
# dev.off()

# pdf("Rplot.pdf")
# p3 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.7)
# dev.off()
# pdf("Rplot.pdf")
# p3 <- plot(fits[[2]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.7)
# dev.off()
# pdf("Rplot.pdf")
# p3 <- plot(fits[[3]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.7)
# dev.off()
# p3 <- plot(fits[[2]]$phi, facet_chains = F)
# p3 <- plot(fits[[3]]$phi, facet_chains = F)

# pdf("Rplot.pdf")
# p3 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
# p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
# dev.off()
