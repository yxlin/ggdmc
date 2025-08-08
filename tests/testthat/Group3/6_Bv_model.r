#  q(save = "no")
cat("\n\n-------------------- Bv model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data6.rda")
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
# mpsrf =  1.056835
#                 loc_A loc_B.left.blue loc_B.left.red loc_B.right.blue
# True            0.250            2.55          2.330             2.35
# 5% Estimate     0.016            2.49          2.226             2.19
# 50% Estimate    0.130            2.58          2.284             2.27
# 97.5% Estimate  0.274            2.68          2.355             2.38
# Median-True    -0.120            0.03         -0.046            -0.08
#                loc_B.right.red loc_mean_v.high.false loc_mean_v.high.true
# True                      2.13                 1.820                2.100
# 5% Estimate               2.17                 1.704                1.938
# 50% Estimate              2.25                 1.755                2.018
# 97.5% Estimate            2.34                 1.809                2.110
# Median-True               0.12                -0.065               -0.082
#                loc_mean_v.low.false loc_mean_v.low.true
# True                          1.500              4.1000
# 5% Estimate                   1.416              3.9436
# 50% Estimate                  1.472              4.1031
# 97.5% Estimate                1.546              4.2797
# Median-True                  -0.028              0.0031
#                loc_mean_v.moderate.false loc_mean_v.moderate.true loc_sd_v.true
# True                              1.7000                    3.100         0.200
# 5% Estimate                       1.6550                    2.949         0.016
# 50% Estimate                      1.7048                    3.079         0.137
# 97.5% Estimate                    1.7635                    3.239         0.273
# Median-True                       0.0048                   -0.021        -0.063
#                 loc_t0 sca_A sca_B.left.blue sca_B.left.red sca_B.right.blue
# True            0.0500 0.200            0.26          0.230            0.235
# 5% Estimate     0.0032 0.185            0.24          0.153            0.224
# 50% Estimate    0.0280 0.258            0.29          0.187            0.274
# 97.5% Estimate  0.0644 0.385            0.41          0.259            0.383
# Median-True    -0.0220 0.058            0.04         -0.043            0.039
#                sca_B.right.red sca_mean_v.high.false sca_mean_v.high.true
# True                     0.230                 0.182                0.250
# 5% Estimate              0.200                 0.123                0.207
# 50% Estimate             0.243                 0.151                0.253
# 97.5% Estimate           0.336                 0.207                0.355
# Median-True              0.013                -0.031                0.003
#                sca_mean_v.low.false sca_mean_v.low.true
# True                           0.15               0.450
# 5% Estimate                    0.15               0.389
# 50% Estimate                   0.18               0.476
# 97.5% Estimate                 0.37               0.660
# Median-True                    0.03               0.026
#                sca_mean_v.moderate.false sca_mean_v.moderate.true sca_sd_v.true
# True                               0.170                     0.35         0.200
# 5% Estimate                        0.130                     0.37         0.176
# 50% Estimate                       0.159                     0.45         0.249
# 97.5% Estimate                     0.231                     0.62         0.378
# Median-True                       -0.011                     0.10         0.049
#                sca_t0
# True            0.050
# 5% Estimate     0.046
# 50% Estimate    0.062
# 97.5% Estimate  0.288
# Median-True     0.012
