# q(save = "no")
cat("\n\n-------------------- v model --------------------")

rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data2.rda")
load(fn)

fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
    sub_migration_prob = 0.00, thin = 1, seed = 123, is_pblocked = TRUE
)

fits1 <- RestartSampling_hyper(fits0, sub_migration_prob = 0.00, thin = 1, seed = 123)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est <- ggdmc::compare(fit, ps = true_vector)
# mpsrf =  1.057453
#                 loc_A loc_B loc_mean_v.conservative.false
# True            0.400 0.500                         0.150
# 5% Estimate     0.360 0.478                         0.013
# 50% Estimate    0.389 0.514                         0.105
# 97.5% Estimate  0.438 0.560                         0.242
# Median-True    -0.011 0.014                        -0.045
#                loc_mean_v.conservative.true loc_mean_v.liberal.false
# True                                  2.150                    0.150
# 5% Estimate                           2.104                    0.005
# 50% Estimate                          2.165                    0.058
# 97.5% Estimate                        2.235                    0.183
# Median-True                           0.015                   -0.092
#                loc_mean_v.liberal.true loc_sd_v.true loc_t0  sca_A  sca_B
# True                             2.500         0.100 0.3000  0.100 0.1000
# 5% Estimate                      2.442         0.011 0.2672  0.072 0.0859
# 50% Estimate                     2.514         0.069 0.3054  0.089 0.1063
# 97.5% Estimate                   2.599         0.134 0.3574  1.057 0.7062
# Median-True                      0.014        -0.031 0.0054 -0.011 0.0063
#                sca_mean_v.conservative.false sca_mean_v.conservative.true
# True                                   0.200                       0.2000
# 5% Estimate                            0.182                       0.1549
# 50% Estimate                           0.249                       0.1904
# 97.5% Estimate                         0.387                       0.2698
# Median-True                            0.049                      -0.0096
#                sca_mean_v.liberal.false sca_mean_v.liberal.true sca_sd_v.true
# True                              0.200                   0.200         0.100
# 5% Estimate                       0.179                   0.180         0.076
# 50% Estimate                      0.232                   0.221         0.107
# 97.5% Estimate                    0.356                   0.333         0.174
# Median-True                       0.032                   0.021         0.007
#                sca_t0
# True            0.100
# 5% Estimate     0.091
# 50% Estimate    0.113
# 97.5% Estimate  0.357
# Median-True     0.013
