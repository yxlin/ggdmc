# q(save = "no")
cat("\n\n--------------------t0 model (17 parameters)--------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data5.rda")
load(fn)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.01, thin = 4, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.02, thin = 4, seed = 9032)
fits2 <- RestartSampling_subject(fits1, sub_migration_prob = 0.00, thin = 4, seed = 9032)

fits <- fits2
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- compare(fit, ps = p_vector)
# mpsrf =  1.009562

#                    A     B mean_v.紅.false mean_v.紅.true mean_v.綠.false
# True           0.250  1.33            1.24           2.50            1.72
# 5% Estimate    0.191  0.82            0.59           1.90            1.02
# 50% Estimate   0.261  1.14            1.13           2.33            1.56
# 97.5% Estimate 0.336  1.46            1.63           2.76            2.07
# Median-True    0.011 -0.19           -0.11          -0.17           -0.16
#                mean_v.綠.true  mean_v.藍.false mean_v.藍.true mean_v.黃.false
# True                     3.10            1.44           2.90            1.52
# 5% Estimate              2.42            0.66           2.27            0.77
# 50% Estimate             2.89            1.22           2.73            1.30
# 97.5% Estimate           3.37            1.74           3.19            1.80
# Median-True             -0.21           -0.22          -0.17           -0.22
#                mean_v.黃.true sd_v.true t0.left_hander.excellent
# True                     2.70     0.200                    0.110
# 5% Estimate              2.08     0.160                    0.102
# 50% Estimate             2.52     0.189                    0.144
# 97.5% Estimate           2.97     0.231                    0.207
# Median-True             -0.18    -0.011                    0.034
#                t0.left_hander.good t0.left_hander.poor
# True                         0.120               0.130
# 5% Estimate                  0.115               0.123
# 50% Estimate                 0.157               0.165
# 97.5% Estimate               0.219               0.227
# Median-True                  0.037               0.035
#                t0.right_hander.excellent t0.right_hander.good
# True                              0.0100                0.020
# 5% Estimate                       0.0057                0.017
# 50% Estimate                      0.0475                0.059
# 97.5% Estimate                    0.1090                0.121
# Median-True                       0.0375                0.039
#                t0.right_hander.poor
# True                          0.030
# 5% Estimate                   0.024
# 50% Estimate                  0.066
# 97.5% Estimate                0.128
# Median-True                   0.036
