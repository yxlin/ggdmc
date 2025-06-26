# q(save = "no")
cat("\n\n---------------B Model (13 parameters)-----------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ggdmc_data3.rda"
load(fn)

model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = c("S", "政黨傾向"), mean_v = "M", sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(S = c("紅", "黃", "藍", "綠"), 政黨傾向 = c("自由派", "保守派")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "lba"
)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.06, thin = 8, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.02, thin = 4, seed = 9032)
fits2 <- RestartSampling_subject(fits1, sub_migration_prob = 0.00, thin = 32, seed = 9032)

fits <- fits2
fit <- RebuildPosterior(fits)

options(digits = 2)
est_phi <- compare(fit, ps = p_vector)

#                    A  B.紅.保守派  B.紅.自由派   B.綠.保守派  B.綠.自由派
# True           0.250        2.40        1.20        4.20        2.10
# 5% Estimate    0.025        2.05        0.91        3.79        1.80
# 50% Estimate   0.269        2.64        1.39        4.55        2.36
# 97.5% Estimate 0.780        3.16        1.78        5.28        2.85
# Median-True    0.019        0.24        0.19        0.35        0.26
#                B.藍.保守派 B.藍.自由派 B.黃.保守派 B.黃.自由派 mean_v.false
# True                  3.60        1.80        3.00        1.50         1.15
# 5% Estimate           3.30        1.50        2.62        1.23         0.99
# 50% Estimate          4.01        2.04        3.26        1.74         1.36
# 97.5% Estimate        4.69        2.50        3.84        2.16         1.73
# Median-True           0.41        0.24        0.26        0.24         0.21
#                mean_v.true sd_v.true      t0
# True                  2.80     0.800  0.1000
# 5% Estimate           2.66     0.711  0.0097
# 50% Estimate          2.98     0.783  0.0614
# 97.5% Estimate        3.34     0.873  0.1605
# Median-True           0.18    -0.017 -0.0386
