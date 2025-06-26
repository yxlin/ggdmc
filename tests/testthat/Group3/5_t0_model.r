# q(save = "no")
cat("\n\n-------------------- t0 model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ggdmc_data5.rda"
save_fn <- "~/Documents/ggdmc/tests/testthat/Group3/data/5_t0_model.rda"
load(fn)

hyper_model <- ggdmcModel:::BuildModel(
    p_map = list(
        A = "1", B = "1",
        t0 = c("HANDEDNESS", "MOTOR_SKILL"),
        mean_v = c("S", "M"), sd_v = "M", st0 = "1"
    ),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(
        S = c("紅", "黃", "藍", "綠"),
        HANDEDNESS = c("left_hander", "right_hander"),
        MOTOR_SKILL = c("excellent", "good", "poor")
    ),
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

#                A_loc   A_sca  B_loc  B_sca    mean_v.紅.false_loc
# True           0.2500  0.2500 1.3300 0.2000              1.2400
# 2.5% Estimate  0.2874  0.1816 1.3034 0.1867              1.2194
# 50% Estimate   0.3383  0.2151 1.3583 0.2195              1.2471
# 97.5% Estimate 0.3887  0.2569 1.4135 0.2674              1.2767
# Median-True    0.0883 -0.0349 0.0283 0.0195              0.0071

#                mean_v.紅.false_sca mean_v.紅.true_loc mean_v.紅.true_sca
# True                        0.1000             2.5000             0.2000
# 2.5% Estimate               0.1018             2.4001             0.1602
# 50% Estimate                0.1197             2.4438             0.1882
# 97.5% Estimate              0.1448             2.4889             0.2280
# Median-True                 0.0197            -0.0562            -0.0118

#                mean_v.綠.false_loc mean_v.綠.false_sca mean_v.綠.true_loc
# True                        1.7200              0.1000             3.1000
# 2.5% Estimate               1.6978              0.0886             2.9832
# 50% Estimate                1.7220              0.1044             3.0576
# 97.5% Estimate              1.7488              0.1250             3.1358
# Median-True                 0.0020              0.0044            -0.0424

#                mean_v.綠.true_sca mean_v.藍.false_loc mean_v.藍.false_sca
# True                       0.3000              1.4400              0.1000
# 2.5% Estimate              0.2584              1.3989              0.0889
# 50% Estimate               0.3034              1.4240              0.1047
# 97.5% Estimate             0.3668              1.4495              0.1259
# Median-True                0.0034             -0.0160              0.0047

#                mean_v.藍.true_loc mean_v.藍.true_sca mean_v.黃.false_loc
# True                       2.9000             0.2000              1.5200
# 2.5% Estimate              2.8417             0.1705              1.4813
# 50% Estimate               2.8912             0.1995              1.5074
# 97.5% Estimate             2.9405             0.2348              1.5317
# Median-True               -0.0088            -0.0005             -0.0126

#                mean_v.黃.false_sca mean_v.黃.true_loc mean_v.黃.true_sca
# True                        0.1000             2.7000             0.2000
# 2.5% Estimate               0.0843             2.6771             0.1565
# 50% Estimate                0.0991             2.7227             0.1845
# 97.5% Estimate              0.1191             2.7702             0.2245
# Median-True                -0.0009             0.0227            -0.0155

#                sd_v.true_loc sd_v.true_sca t0.left_hander.excellent_loc
# True                  0.2000        0.2000                       0.0300
# 2.5% Estimate         0.2071        0.1551                       0.1353
# 50% Estimate          0.2497        0.1814                       0.1682
# 97.5% Estimate        0.2980        0.2176                       0.1999
# Median-True           0.0497       -0.0186                       0.1382

#                t0.left_hander.excellent_sca t0.left_hander.good_loc
# True                                 0.2000                  0.0500
# 2.5% Estimate                        0.1080                  0.1152
# 50% Estimate                         0.1261                  0.1401
# 97.5% Estimate                       0.1494                  0.1644
# Median-True                         -0.0739                  0.0901

#                t0.left_hander.good_sca t0.left_hander.poor_loc
# True                            0.2000                  0.1200
# 2.5% Estimate                   0.0859                  0.1623
# 50% Estimate                    0.1012                  0.1934
# 97.5% Estimate                  0.1222                  0.2246
# Median-True                    -0.0988                  0.0734

#                t0.left_hander.poor_sca t0.right_hander.excellent_loc
# True                            0.2000                        0.0100
# 2.5% Estimate                   0.1075                        0.1446
# 50% Estimate                    0.1262                        0.1810
# 97.5% Estimate                  0.1533                        0.2153
# Median-True                    -0.0738                        0.1710
#                t0.right_hander.excellent_sca t0.right_hander.good_loc
# True                                  0.2000                   0.0600
# 2.5% Estimate                         0.1206                   0.1705
# 50% Estimate                          0.1418                   0.2065
# 97.5% Estimate                        0.1709                   0.2447
# Median-True                          -0.0582                   0.1465
#                t0.right_hander.good_sca t0.right_hander.poor_loc
# True                             0.2000                   0.1000
# 2.5% Estimate                    0.1275                   0.1454
# 50% Estimate                     0.1508                   0.1771
# 97.5% Estimate                   0.1817                   0.2065
# Median-True                     -0.0492                   0.0771
#                t0.right_hander.poor_sca
# True                             0.2000
# 2.5% Estimate                    0.1108
# 50% Estimate                     0.1310
# 97.5% Estimate                   0.1590
# Median-True                     -0.0690
