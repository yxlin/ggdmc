# q(save = "no")
cat("\n\n-------------------- 6 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data1.rda"
load(fn)

model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = "M", sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "lba"
)
# sub_dmis[[1]]@is_positive_drift
# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.06, thin = 1, nmc = 3, seed = 9032)
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.06, thin = 8, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 8, seed = 9032)

fits <- fits1
fit <- RebuildPosterior(fits)

options(digits = 2)
est_phi <- compare(fit, ps = p_vector)
#                    A     B mean_v.false mean_v.true sd_v.true     t0
# True            0.750 1.250         1.50       2.500     0.100  0.150
# 5% Estimate     0.415 0.572         0.39       1.493     0.070  0.013
# 50% Estimate    0.699 1.325         1.66       2.489     0.117  0.136
# 97.5% Estimate  1.074 2.342         2.91       3.681     0.193  0.310
# Median-True    -0.051 0.075         0.16      -0.011     0.017 -0.014


#                   A    B mean_v.false mean_v.true sd_v.true      t0
# True           0.75 1.25         1.50        2.50    0.1000  0.1500
# 5% Estimate    0.51 0.56         0.42        1.65    0.0226  0.0127
# 50% Estimate   0.97 1.49         1.96        2.95    0.0989  0.1406
# 97.5% Estimate 2.00 3.44         4.78        5.77    0.2902  0.3874
# Median-True    0.22 0.24         0.46        0.45   -0.0011 -0.0094

#  A   B mean_v.false mean_v.true sd_v.true      t0
# True           0.75 1.2          1.5         2.5       0.1  0.1500
# 5% Estimate    0.99 3.0          6.3         6.2       1.0  0.0091
# 50% Estimate   2.43 4.9          8.7         8.6       1.6  0.1053
# 97.5% Estimate 3.71 6.6          9.9         9.9       2.8  0.3098
# Median-True    1.68 3.7          7.2         6.1       1.5 -0.0447




# ggdmc::plot(fits[[1]], start = fits[[1]]@nmc * 0.5)
# ggdmc::plot(fits[[2]], start = fits[[1]]@nmc * 0.5)
# ggdmc::plot(fits[[3]], start = fits[[1]]@nmc * 0.5)
# ggdmc::plot(fit, pll = F, den = T)
# gdmc:::plot(fits[[1]])
# gdmc:::plot(fits[[2]])
# gdmc:::plot(fits[[3]])
# gdmc:::plot(fits[[1]], pll = F, den = T)
# gdmc:::plot(fits[[2]], pll = F, den = T)
# gdmc:::plot(fits[[3]], pll = F, den = T)
