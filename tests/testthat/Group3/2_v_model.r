# q(save = "no")
cat("\n\n-------------------- v model --------------------")

rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ggdmc_data2.rda"
save_fn <- "~/Documents/ggdmc/tests/testthat/Group3/data/2_v_model.rda"
load(fn)

hyper_model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = c("POLITICAL_VIEW", "M"), sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2"), POLITICAL_VIEW = c("liberal", "conservative")),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("r1", "r2"),
    type = "hyper",
    verbose = FALSE
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

#                loc_A    loc_B loc_mean_v.conservative.false
# True           0.400  0.50000                         0.150
# 5% Estimate    0.381  0.46193                         0.038
# 50% Estimate   0.423  0.49909                         0.171
# 97.5% Estimate 0.471  0.54461                         0.252
# Median-True    0.023 -0.00091                         0.021
#                loc_mean_v.conservative.true loc_mean_v.liberal.false
# True                                  2.150                     0.15
# 5% Estimate                           2.106                     0.14
# 50% Estimate                          2.163                     0.26
# 97.5% Estimate                        2.231                     0.33
# Median-True                           0.013                     0.11
#                loc_mean_v.liberal.true loc_sd_v.true loc_t0 sca_A sca_B
# True                            2.5000         0.100  0.300 0.100 0.100
# 5% Estimate                     2.4245         0.014  0.294 0.111 0.104
# 50% Estimate                    2.4958         0.084  0.319 0.137 0.126
# 97.5% Estimate                  2.5806         0.137  0.349 0.186 0.168
# Median-True                    -0.0042        -0.016  0.019 0.037 0.026
#                sca_mean_v.conservative.false sca_mean_v.conservative.true
# True                                   0.200                       0.2000
# 5% Estimate                            0.130                       0.1580
# 50% Estimate                           0.182                       0.1933
# 97.5% Estimate                         0.288                       0.2554
# Median-True                           -0.018                      -0.0067
#                sca_mean_v.liberal.false sca_mean_v.liberal.true sca_sd_v.true
# True                              0.200                   0.200         0.100
# 5% Estimate                       0.132                   0.199         0.078
# 50% Estimate                      0.174                   0.244         0.110
# 97.5% Estimate                    0.289                   0.320         0.167
# Median-True                      -0.026                   0.044         0.010
#                sca_t0
# True            0.100
# 5% Estimate     0.069
# 50% Estimate    0.085
# 97.5% Estimate  0.113
# Median-True    -0.015
