# q(save = "no")
cat("\n\n-------------------- 6 parameters --------------------")

rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data1.rda"
load(fn)


fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
    sub_migration_prob = 0.01, thin = 2
)

fits1 <- RestartSampling_hyper(fits0, thin = 8)
# fits2 <- RestartSampling_hyper(fits1, thin = 4)

fits <- fits1
fit <- RebuildPosterior(fits)

options(digits = 2)
est <- compare(fit, ps = true_vector)
# save(fits0, fits1, fits2, est, file = save_fn)
#                 loc_A   loc_B loc_mean_v.false loc_mean_v.true loc_sd_v.true
# True            0.400  0.5000           0.1500          2.5000        0.1000
# 5% Estimate     0.351  0.4675           0.0069          2.4354        0.0201
# 50% Estimate    0.382  0.4988           0.0738          2.4932        0.0902
# 97.5% Estimate  0.419  0.5370           0.2011          2.5636        0.1351
# Median-True    -0.018 -0.0012          -0.0762         -0.0068       -0.0098
#                 loc_t0  sca_A  sca_B sca_mean_v.false sca_mean_v.true
# True            0.3000 0.1000 0.1000             0.20          0.2000
# 5% Estimate     0.2779 0.0844 0.0865             0.19          0.1625
# 50% Estimate    0.2921 0.1029 0.1058             0.25          0.1984
# 97.5% Estimate  0.3094 0.1372 0.1402             0.34          0.2609
# Median-True    -0.0079 0.0029 0.0058             0.05         -0.0016
#                sca_sd_v.true  sca_t0
# True                  0.1000  0.0500
# 5% Estimate           0.0709  0.0395
# 50% Estimate          0.0989  0.0481
# 97.5% Estimate        0.1554  0.0634
# Median-True          -0.0011 -0.0019
