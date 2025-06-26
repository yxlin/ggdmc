#  q(save = "no")
cat("\n\n--------------------Testing A Model--------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ggdmc_data4.rda"
load(fn)

model <- ggdmcModel::BuildModel(
    p_map = list(B = "1", A = c("S", "政黨傾向"), mean_v = "M", sd_v = "M", st0 = "1", t0 = "1"),
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

#                 A.紅.保守派   A.紅.自由派  A.綠.保守派  A.綠.自由派   A.藍.保守派
# True                1.2300      0.2500      1.8900      0.5500      1.6700
# 2.5% Estimate       1.1794      0.0200      1.7097      0.2186      1.5452
# 50% Estimate        1.7683      0.2695      2.4424      0.5827      2.2787
# 97.5% Estimate      2.5932      0.7415      3.4335      1.0873      3.3097
# Median-True         0.5383      0.0195      0.5524      0.0327      0.6087
#                 A.藍.自由派   A.黃.保守派   A.黃.自由派     B  mean_v.false
# True                0.4500      1.4500      0.3500 1.2500       1.1500
# 2.5% Estimate       0.3163      1.2844      0.0886 1.1985       1.2056
# 50% Estimate        0.6773      1.8949      0.4225 1.7197       1.8197
# 97.5% Estimate      1.1783      2.7140      0.9390 2.1732       2.4192
# Median-True         0.2273      0.4449      0.0725 0.4697       0.6697
#                mean_v.true sd_v.true      t0
# True                2.8000    0.8000  0.1000
# 2.5% Estimate       2.8901    0.7583  0.0031
# 50% Estimate        3.5166    0.9074  0.0541
# 97.5% Estimate      4.2135    1.0879  0.1409
# Median-True         0.7166    0.1074 -0.0459
