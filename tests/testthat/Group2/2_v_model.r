# q(save = "no")
cat("\n\n--------------------Drift Rate Model--------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data2.rda"
load(fn)

model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = c("POLITICAL_VIEW", "M"), sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2"), POLITICAL_VIEW = c("liberal", "conservative")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "lba"
)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.00, thin = 8, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 4, seed = 9032)
fits2 <- RestartSampling_subject(fits1, sub_migration_prob = 0.00, thin = 32, seed = 9032)
fits3 <- RestartSampling_subject(fits2, sub_migration_prob = 0.00, thin = 32, seed = 9032)

fits <- fits3
fit <- RebuildPosterior(fits)

options(digits = 2)
est_phi <- compare(fit, ps = p_vector)
#                   A    B mean_v.conservative.false mean_v.conservative.true
# True           0.75 1.25                      1.50                     2.15
# 5% Estimate    0.57 0.79                      0.87                     1.58
# 50% Estimate   0.86 1.51                      1.87                     2.39
# 97.5% Estimate 1.22 2.32                      2.86                     3.30
# Median-True    0.11 0.26                      0.37                     0.24
#                mean_v.liberal.false mean_v.liberal.true sd_v.true    t0
# True                           1.50                2.50    0.1000  0.15
# 5% Estimate                    0.59                1.88    0.0536  0.01
# 50% Estimate                   1.65                2.77    0.0987  0.10
# 97.5% Estimate                 2.71                3.78    0.1785  0.26
# Median-True                    0.15                0.27   -0.0013 -0.05
plot(fits[[1]])
plot(fits[[2]])
plot(fits[[3]])
plot(fit, den = T, pll = F)
# gdmc:::plot(fits[[1]])
# gdmc:::plot(fits[[2]])
# gdmc:::plot(fits[[3]])
# gdmc:::plot(fits[[1]], pll = F, den = T)
# gdmc:::plot(fits[[2]], pll = F, den = T)
# gdmc:::plot(fits[[3]], pll = F, den = T)
# gdmc:::plot(fit, pll = F, den = T)

# est <- gdmc:::summary(fit, start = fit@nmc * 0.5, recovery = TRUE, ps = p_vector, verbose = TRUE)
# est1 <- gdmc:::summary(fits[[1]], start = 1, recovery = TRUE, ps = p_vector, verbose = TRUE)
# est2 <- gdmc:::summary(fits[[2]], start = 1, recovery = TRUE, ps = p_vector, verbose = TRUE)
# est3 <- gdmc:::summary(fits[[3]], start = 1, recovery = TRUE, ps = p_vector, verbose = TRUE)
#                     A       B  mean_v.conservative.false
# True            0.7500  1.2500                    1.5000
# 2.5% Estimate   0.4397  0.4883                    0.1896
# 50% Estimate    0.6998  1.0331                    1.2874
# 97.5% Estimate  1.1271  2.2026                    2.6059
# Median-True    -0.0502 -0.2169                   -0.2126
#
#                mean_v.conservative.true mean_v.liberal.false
# True                             2.1500               1.5000
# 2.5% Estimate                    1.2464               0.1659
# 50% Estimate                     1.9935               1.2724
# 97.5% Estimate                   3.2094               2.6683
# Median-True                     -0.1565              -0.2276
#
#                mean_v.liberal.true sd_v.true     t0
# True                        2.5000    0.1000 0.1500
# 2.5% Estimate               1.5131    0.0762 0.0188
# 50% Estimate                2.3460    0.1313 0.2091
# 97.5% Estimate              3.6602    0.2193 0.3435
# Median-True                -0.1540    0.0313 0.0591
