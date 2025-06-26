# q(save = "no")
cat("\n\n--------------------BV Model--------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ggdmc_data6.rda"
load(fn)

model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = c("S", "COLOR"), mean_v = c("NOISE", "M"), sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(left = "z_key", right = "x_key")),
    factors = list(
        S = c("left", "right"),
        COLOR = c("red", "blue"),
        NOISE = c("high", "moderate", "low")
    ),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("z_key", "x_key"),
    type = "lba"
)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.06, thin = 8, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.02, thin = 4, seed = 9032)
fits2 <- RestartSampling_subject(fits1, sub_migration_prob = 0.00, thin = 32, seed = 9032)

fits <- fits2
fit <- RebuildPosterior(fits)

options(digits = 2)
est_phi <- compare(fit, ps = p_vector)

save_fn <- "~/Documents/ggdmc/tests/testthat/Group2/data/6_Bv_model.rda"
save(fits0, fits1, fits2, est_phi, file = save_fn)
#                    A B.left.blue B.left.red B.right.blue B.right.red
# True           0.250        2.55       2.33         2.35        2.13
# 5% Estimate    0.123        2.43       2.25         2.24        2.03
# 50% Estimate   0.342        3.23       3.00         2.99        2.72
# 97.5% Estimate 0.542        4.26       3.96         3.95        3.61
# Median-True    0.092        0.68       0.67         0.64        0.59
#                mean_v.high.false mean_v.high.true mean_v.low.false
# True                        1.82             2.50              1.5
# 5% Estimate                 1.45             2.44              1.4
# 50% Estimate                2.29             3.17              2.9
# 97.5% Estimate              3.34             4.14              4.6
# Median-True                 0.47             0.67              1.4
#                mean_v.low.true mean_v.moderate.false mean_v.moderate.true
# True                       4.5                  1.70                 3.50
# 5% Estimate                4.3                  1.26                 3.35
# 50% Estimate               5.5                  2.37                 4.33
# 97.5% Estimate             7.2                  3.70                 5.64
# Median-True                1.0                  0.67                 0.83
#                sd_v.true      t0
# True               0.200  0.0500
# 5% Estimate        0.165  0.0019
# 50% Estimate       0.222  0.0235
# 97.5% Estimate     0.304  0.0969
# Median-True        0.022 -0.0265



#                     A  B.left.blue B.left.red B.right.blue B.right.red
# True           0.2500      2.5500     2.3300       2.3500      2.1300
# 2.5% Estimate  0.0158      1.7933     1.5628       1.5629      1.3884
# 50% Estimate   0.2627      2.7997     2.5208       2.5258      2.2992
# 97.5% Estimate 0.4796      3.9261     3.5863       3.5949      3.3135
# Median-True    0.0127      0.2497     0.1908       0.1758      0.1692
#                mean_v.high.false mean_v.high.true mean_v.low.false
# True                      1.8200           2.5000           1.5000
# 2.5% Estimate             1.4646           2.0796           0.0650
# 50% Estimate              2.3455           2.8276           1.3159
# 97.5% Estimate            3.2737           3.6609           3.3624
# Median-True               0.5255           0.3276          -0.1841
#                mean_v.low.true mean_v.moderate.false mean_v.moderate.true
# True                    4.5000                1.7000               3.5000
# 2.5% Estimate           4.3391                1.1604               3.1546
# 50% Estimate            5.5209                2.2630               4.1155
# 97.5% Estimate          6.8297                3.4052               5.1906
# Median-True             1.0209                0.5630               0.6155
#                sd_v.true     t0
# True              0.2000 0.0500
# 2.5% Estimate     0.2081 0.0100
# 50% Estimate      0.2768 0.1088
# 97.5% Estimate    0.3556 0.2403
# Median-True       0.0768 0.0588
