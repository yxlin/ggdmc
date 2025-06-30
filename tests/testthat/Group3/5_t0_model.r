# q(save = "no")
cat("\n\n-------------------- t0 model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data5.rda")
load(fn)

fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
    sub_migration_prob = 0.00, thin = 1, seed = 123, is_pblocked = TRUE
)

fits1 <- RestartSampling_hyper(fits0, sub_migration_prob = 0.00, thin = 1, seed = 123)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est <- ggdmc::compare(fit, ps = true_vector)
# mpsrf <- 1.054621
#                loc_A   loc_B loc_mean_v.紅.false loc_mean_v.紅.true
# True           0.250  1.3300               1.240              2.500
# 5% Estimate    0.181  1.2654               1.161              2.368
# 50% Estimate   0.297  1.3229               1.193              2.435
# 97.5% Estimate 0.374  1.3913               1.230              2.511
# Median-True    0.047 -0.0071              -0.047             -0.065
#                loc_mean_v.綠.false loc_mean_v.綠.true loc_mean_v.藍.false
# True                         1.720              3.100               1.440
# 5% Estimate                  1.666              3.058               1.383
# 50% Estimate                 1.699              3.169               1.414
# 97.5% Estimate               1.736              3.302               1.450
# Median-True                 -0.021              0.069              -0.026
#                loc_mean_v.藍.true loc_mean_v.黃.false loc_mean_v.黃.true
# True                        2.900               1.520             2.7000
# 5% Estimate                 2.812               1.521             2.6439
# 50% Estimate                2.869               1.558             2.7093
# 97.5% Estimate              2.931               1.605             2.7878
# Median-True                -0.031               0.038             0.0093
#                loc_sd_v.true loc_t0.left_hander.excellent
# True                   0.200                        0.110
# 5% Estimate            0.014                        0.043
# 50% Estimate           0.105                        0.172
# 97.5% Estimate         0.215                        0.242
# Median-True           -0.095                        0.062
#                loc_t0.left_hander.good loc_t0.left_hander.poor
# True                             0.120                  0.1300
# 5% Estimate                      0.034                  0.0065
# 50% Estimate                     0.165                  0.0718
# 97.5% Estimate                   0.260                  0.2011
# Median-True                      0.045                 -0.0582
#                loc_t0.right_hander.excellent loc_t0.right_hander.good
# True                                  0.0100                   0.0200
# 5% Estimate                           0.0035                   0.0062
# 50% Estimate                          0.0382                   0.0740
# 97.5% Estimate                        0.1568                   0.1922
# Median-True                           0.0282                   0.0540
#                loc_t0.right_hander.poor  sca_A  sca_B sca_mean_v.紅.false
# True                             0.0300  0.250  0.200              0.1000
# 5% Estimate                      0.0043  0.139  0.151              0.0779
# 50% Estimate                     0.0508  0.183  0.184              0.0955
# 97.5% Estimate                   0.1611  0.320  0.270              0.1700
# Median-True                      0.0208 -0.067 -0.016             -0.0045
#                sca_mean_v.紅.true sca_mean_v.綠.false sca_mean_v.綠.true
# True                        0.200              0.1000              0.300
# 5% Estimate                 0.165              0.0835              0.308
# 50% Estimate                0.203              0.1028              0.376
# 97.5% Estimate              0.278              0.2950              0.554
# Median-True                 0.003              0.0028              0.076
#                sca_mean_v.藍.false sca_mean_v.藍.true sca_mean_v.黃.false
# True                        0.1000              0.200               0.100
# 5% Estimate                 0.0751              0.143               0.097
# 50% Estimate                0.0929              0.175               0.118
# 97.5% Estimate              0.1443              0.241               0.180
# Median-True                -0.0071             -0.025               0.018
#                sca_mean_v.黃.true sca_sd_v.true sca_t0.left_hander.excellent
# True                       0.2000        0.2000                        0.200
# 5% Estimate                0.1677        0.1394                        0.115
# 50% Estimate               0.2055        0.1963                        0.161
# 97.5% Estimate             0.2922        0.3121                        0.290
# Median-True                0.0055       -0.0037                       -0.039
#                sca_t0.left_hander.good sca_t0.left_hander.poor
# True                           2.0e-01                   0.200
# 5% Estimate                    1.4e-01                   0.181
# 50% Estimate                   2.0e-01                   0.242
# 97.5% Estimate                 3.4e-01                   0.338
# Median-True                    2.3e-05                   0.042
#                sca_t0.right_hander.excellent sca_t0.right_hander.good
# True                                    0.20                   0.2000
# 5% Estimate                             0.20                   0.1525
# 50% Estimate                            0.25                   0.2077
# 97.5% Estimate                          0.34                   0.3446
# Median-True                             0.05                   0.0077
#                sca_t0.right_hander.poor
# True                              0.200
# 5% Estimate                       0.139
# 50% Estimate                      0.184
# 97.5% Estimate                    0.269
# Median-True                      -0.016
