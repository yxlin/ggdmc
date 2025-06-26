#  q(save = "no")
cat("\n\n-------------------- Bv model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data6.rda"
save_fn <- "~/Documents/ggdmc/tests/testthat/Group3/fit_data/6_BV.rda"

load(fn)


fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors, sub_migration_prob = 0.01, thin = 8, ncore = 5)
fits1 <- RestartSampling_hyper(fits0, thin = 16)
fits2 <- RestartSampling_hyper(fits1, thin = 4)

fits <- fits2
fit <- RebuildPosterior(fits)

options(digits = 2)
est <- compare(fit, ps = true_vector)
save(fits0, fits1, fits2, est, file = save_fn)


#          loc_A loc_B.left.blue loc_B.left.red loc_B.right.blue
# True           0.250           2.550          2.330            2.350
# 5% Estimate    0.202           2.464          2.285            2.259
# 50% Estimate   0.297           2.535          2.352            2.319
# 97.5% Estimate 0.365           2.621          2.436            2.388
# Median-True    0.047          -0.015          0.022           -0.031
#                loc_B.right.red loc_mean_v.high.false loc_mean_v.high.true
# True                      2.13                1.8200                 2.10
# 5% Estimate               2.18                1.7573                 2.03
# 50% Estimate              2.24                1.8122                 2.11
# 97.5% Estimate            2.31                1.8769                 2.21
# Median-True               0.11               -0.0078                 0.01
#                loc_mean_v.low.false loc_mean_v.low.true
# True                          1.500                 4.1
# 5% Estimate                   1.419                 3.8
# 50% Estimate                  1.466                 3.9
# 97.5% Estimate                1.521                 4.1
# Median-True                  -0.034                -0.2
#                loc_mean_v.moderate.false loc_mean_v.moderate.true loc_sd_v.true
# True                               1.700                   3.1000        0.2000
# 5% Estimate                        1.631                   2.9720        0.0409
# 50% Estimate                       1.683                   3.0952        0.1909
# 97.5% Estimate                     1.747                   3.2406        0.2907
# Median-True                       -0.017                  -0.0048       -0.0091
#                 loc_t0  sca_A sca_B.left.blue sca_B.left.red sca_B.right.blue
# True            0.0500  0.200           0.255         0.2300            0.235
# 5% Estimate     0.0063  0.139           0.197         0.1898            0.164
# 50% Estimate    0.0428  0.181           0.241         0.2317            0.200
# 97.5% Estimate  0.0735  0.292           0.319         0.3065            0.267
# Median-True    -0.0072 -0.019          -0.014         0.0017           -0.035
#                sca_B.right.red sca_mean_v.high.false sca_mean_v.high.true
# True                     0.230                0.1820                0.250
# 5% Estimate              0.170                0.1509                0.223
# 50% Estimate             0.207                0.1848                0.273
# 97.5% Estimate           0.273                0.2433                0.362
# Median-True             -0.023                0.0028                0.023
#                sca_mean_v.low.false sca_mean_v.low.true
# True                         0.1500               0.450
# 5% Estimate                  0.1299               0.353
# 50% Estimate                 0.1588               0.432
# 97.5% Estimate               0.2091               0.566
# Median-True                  0.0088              -0.018
#                sca_mean_v.moderate.false sca_mean_v.moderate.true sca_sd_v.true
# True                              0.1700                    0.350         0.200
# 5% Estimate                       0.1441                    0.347         0.156
# 50% Estimate                      0.1752                    0.422         0.218
# 97.5% Estimate                    0.2309                    0.556         0.334
# Median-True                       0.0052                    0.072         0.018
#                sca_t0
# True            0.050
# 5% Estimate     0.045
# 50% Estimate    0.063
# 97.5% Estimate  0.094
# Median-True     0.013
