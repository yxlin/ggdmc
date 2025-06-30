# q(save = "no")
cat("\n\n-------------------- 6 parameters --------------------")

rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data1.rda")
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
# mpsrf =  1.080905
#                 loc_A  loc_B loc_mean_v.false loc_mean_v.true loc_sd_v.true
# True            0.400  0.500            0.150           2.500        0.1000
# 5% Estimate     0.343  0.459            0.014           2.456        0.0057
# 50% Estimate    0.369  0.498            0.124           2.515        0.0519
# 97.5% Estimate  0.410  0.541            0.272           2.593        0.1266
# Median-True    -0.031 -0.002           -0.026           0.015       -0.0481
#                loc_t0  sca_A   sca_B sca_mean_v.false sca_mean_v.true
# True           0.3000  0.100  0.1000            0.200           0.200
# 5% Estimate    0.2883  0.066  0.0790            0.172           0.151
# 50% Estimate   0.3055  0.081  0.0973            0.241           0.186
# 97.5% Estimate 2.3480  0.526  1.2933            0.391           0.747
# Median-True    0.0055 -0.019 -0.0027            0.041          -0.014
#                sca_sd_v.true sca_t0
# True                   0.100  0.050
# 5% Estimate            0.085  0.043
# 50% Estimate           0.116  0.053
# 97.5% Estimate         0.176  2.858
# Median-True            0.016  0.003
