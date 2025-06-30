# q(save = "no")
cat("\n\n--------------------BV Model--------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data6.rda")
load(fn)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.03, thin = 2, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 4, seed = 9032)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- compare(fit, ps = p_vector)
# mpsrf = 1.057539
#                     A B.left.blue B.left.red B.right.blue B.right.red
# True            0.250        2.55       2.33         2.35        2.13
# 5% Estimate     0.031        1.76       1.58         1.60        1.40
# 50% Estimate    0.237        2.39       2.18         2.20        1.98
# 97.5% Estimate  0.412        3.06       2.82         2.85        2.59
# Median-True    -0.013       -0.16      -0.15        -0.15       -0.15
#                mean_v.high.false mean_v.high.true mean_v.low.false
# True                        1.82             2.10             1.50
# 5% Estimate                 1.24             1.63             0.24
# 50% Estimate                1.71             1.99             1.11
# 97.5% Estimate              2.24             2.41             2.09
# Median-True                -0.11            -0.11            -0.39
#                mean_v.low.true mean_v.moderate.false mean_v.moderate.true
# True                     4.100                1.7000                3.100
# 5% Estimate              3.514                1.1433                2.562
# 50% Estimate             4.054                1.6932                3.029
# 97.5% Estimate           4.784                2.3218                3.609
# Median-True             -0.046               -0.0068               -0.071
#                sd_v.true    t0
# True              0.2000 0.050
# 5% Estimate       0.1658 0.011
# 50% Estimate      0.2063 0.082
# 97.5% Estimate    0.2518 0.227
# Median-True       0.0063 0.032
