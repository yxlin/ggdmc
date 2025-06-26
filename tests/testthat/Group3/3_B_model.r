# q(save = "no")
cat("\n\n-------------------- B model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ggdmc_data3.rda"
save_fn <- "~/Documents/ggdmc/tests/testthat/Group3/data/3_B_model.rda"
load(fn)

hyper_model <- ggdmcModel:::BuildModel(
    p_map = list(A = "1", B = c("S", "政黨傾向"), t0 = "1", mean_v = "M", sd_v = "M", st0 = "1"),
    factors = list(S = c("紅", "黃", "藍", "綠"), 政黨傾向 = c("自由派", "保守派")),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "hyper",
    verbose = FALSE
)

hyper_model <- ggdmcModel:::BuildModel(
    p_map = list(A = "1", B = c("S", "政黨傾向"), mean_v = "M", sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(S = c("紅", "黃", "藍", "綠"), 政黨傾向 = c("自由派", "保守派")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "hyper"
)

hyper_model@pnames

ls()
hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)


fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors, sub_migration_prob = 0.01, thin = 8, ncore = 3)
fits1 <- RestartSampling_hyper(fits0, thin = 16)
fits2 <- RestartSampling_hyper(fits1, thin = 4)

fits <- fits2
fit <- RebuildPosterior(fits)

options(digits = 2)
est <- compare(fit, ps = true_vector)
save(fits0, fits1, fits2, est, file = save_fn)

#   loc_A loc_B.紅.保守派 loc_B.紅.自由派 loc_B.綠.保守派
# True            0.2500          2.4000          1.2000          4.2000
# 2.5% Estimate   0.2264          2.3760          1.1617          4.1127
# 50% Estimate    0.2442          2.4130          1.1785          4.1826
# 97.5% Estimate  0.2618          2.4508          1.1951          4.2521
# Median-True    -0.0058          0.0130         -0.0215         -0.0174
#                loc_B.綠.自由派 loc_B.藍.保守派 loc_B.藍.自由派 loc_B.黃.保守派
# True                    2.1000          3.6000          1.8000          3.0000
# 2.5% Estimate           2.0489          3.5411          1.7732          2.9453
# 50% Estimate            2.0844          3.5923          1.7909          2.9948
# 97.5% Estimate          2.1197          3.6455          1.8086          3.0423
# Median-True            -0.0156         -0.0077         -0.0091         -0.0052
#                loc_B.黃.自由派 loc_mean_v.false loc_mean_v.true loc_sd_v.true
# True                    1.5000           1.1500          2.8000        0.8000
# 2.5% Estimate           1.4930           1.1304          2.7606        0.7922
# 50% Estimate            1.5116           1.1477          2.7914        0.8098
# 97.5% Estimate          1.5306           1.1650          2.8221        0.8279
# Median-True             0.0116          -0.0023         -0.0086        0.0098
#                 loc_t0  sca_A sca_B.紅.保守派 sca_B.紅.自由派 sca_B.綠.保守派
# True            0.1000 0.1000          0.2000          0.1000          0.4000
# 2.5% Estimate   0.0982 0.0906          0.1892          0.0858          0.3590
# 50% Estimate    0.0999 0.1018          0.2127          0.0966          0.4034
# 97.5% Estimate  0.1016 0.1153          0.2429          0.1095          0.4593
# Median-True    -0.0001 0.0018          0.0127         -0.0034          0.0034
#                sca_B.綠.自由派 sca_B.藍.保守派 sca_B.藍.自由派 sca_B.黃.保守派
# True                    0.2000          0.3000          0.1000          0.3000
# 2.5% Estimate           0.1821          0.2696          0.0918          0.2488
# 50% Estimate            0.2042          0.3036          0.1034          0.2802
# 97.5% Estimate          0.2315          0.3456          0.1173          0.3191
# Median-True             0.0042          0.0036          0.0034         -0.0198
#                sca_B.黃.自由派 sca_mean_v.false sca_mean_v.true sca_sd_v.true
# True                    0.1000           0.1000          0.2000        0.1000
# 2.5% Estimate           0.0971           0.0870          0.1590        0.0909
# 50% Estimate            0.1095           0.0979          0.1785        0.1023
# 97.5% Estimate          0.1242           0.1113          0.2026        0.1168
# Median-True             0.0095          -0.0021         -0.0215        0.0023
#                 sca_t0
# True            0.0100
# 2.5% Estimate   0.0088
# 50% Estimate    0.0099
# 97.5% Estimate  0.0113
# Median-True    -0.0001
