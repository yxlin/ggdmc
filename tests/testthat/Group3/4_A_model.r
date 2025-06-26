# q(save = "no")
cat("\n\n-------------------- A model --------------------")

rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ggdmc_data4.rda"
save_fn <- "~/Documents/ggdmc/tests/testthat/Group3/data/4_A_model.rda"
load(fn)

hyper_model <- ggdmcModel:::BuildModel(
    p_map = list(A = c("S", "政黨傾向"), B = "1", t0 = "1", mean_v = "M", sd_v = "M", st0 = "1"),
    factors = list(S = c("紅", "黃", "藍", "綠"), 政黨傾向 = c("自由派", "保守派")),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "hyper", verbose = FALSE
)


hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)


fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors, sub_migration_prob = 0.01, thin = 8, ncore = 3)
fits1 <- RestartSampling_hyper(fits0, thin = 16)
fits2 <- RestartSampling_hyper(fits1, thin = 4)

fits <- fits2
fit <- RebuildPosterior(fits)

options(digits = 2)
est <- compare(fit, ps = true_vector)
save(fits0, fits1, fits2, est, file = save_fn)
#        A.紅.保守派_loc A.紅.保守派_sca A.紅.自由派_loc A.紅.自由派_sca
# True                    1.2300          0.2000          0.2500          0.1000
# 2.5% Estimate           1.1759          0.1539          0.2189          0.0814
# 50% Estimate            1.2222          0.1840          0.2425          0.0952
# 97.5% Estimate          1.2683          0.2227          0.2656          0.1148
# Median-True            -0.0078         -0.0160         -0.0075         -0.0048

#                A.綠.保守派_loc A.綠.保守派_sca A.綠.自由派_loc A.綠.自由派_sca
# True                    1.8900          0.2000          0.5500          0.1000
# 2.5% Estimate           1.8085          0.1675          0.5202          0.0852
# 50% Estimate            1.8568          0.1994          0.5452          0.1008
# 97.5% Estimate          1.9060          0.2385          0.5695          0.1217
# Median-True            -0.0332         -0.0006         -0.0048          0.0008

#                A.藍.保守派_loc A.藍.保守派_sca A.藍.自由派_loc A.藍.自由派_sca
# True                    1.6700          0.2000          0.4500          0.1000
# 2.5% Estimate           1.5866          0.1586          0.4209          0.0767
# 50% Estimate            1.6344          0.1879          0.4447          0.0905
# 97.5% Estimate          1.6850          0.2285          0.4679          0.1100
# Median-True            -0.0356         -0.0121         -0.0053         -0.0095

#                A.黃.保守派_loc A.黃.保守派_sca A.黃.自由派_loc A.黃.自由派_sca
# True                    1.4500          0.2000          0.3500          0.1000
# 2.5% Estimate           1.3613          0.1658          0.3341          0.0904
# 50% Estimate            1.4097          0.1963          0.3620          0.1057
# 97.5% Estimate          1.4558          0.2373          0.3895          0.1255
# Median-True            -0.0403         -0.0037          0.0120          0.0057

#                 B_loc  B_sca mean_v.false_loc mean_v.false_sca mean_v.true_loc
# True           1.2500 0.1000           1.1500           0.1000          2.8000
# 2.5% Estimate  1.2254 0.0994           1.1103           0.0876          2.7730
# 50% Estimate   1.2536 0.1174           1.1347           0.1036          2.7944
# 97.5% Estimate 1.2841 0.1419           1.1610           0.1248          2.8153
# Median-True    0.0036 0.0174          -0.0153           0.0036         -0.0056

#                mean_v.true_sca sd_v.true_loc sd_v.true_sca t0_loc  t0_sca
# True                    0.1000        0.8000        0.1000 0.1000  0.1000
# 2.5% Estimate           0.0763        0.7800        0.0825 0.0976  0.0675
# 50% Estimate            0.0902        0.8049        0.0975 0.1182  0.0794
# 97.5% Estimate          0.1089        0.8293        0.1169 0.1382  0.0948
# Median-True            -0.0098        0.0049       -0.0025 0.0182 -0.0206
