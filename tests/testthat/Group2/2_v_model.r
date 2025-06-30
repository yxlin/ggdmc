#  q(save = "no")
cat("\n\n--------------------Drift Rate Model--------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data2.rda")
load(fn)


fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.01, thin = 2, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.01, thin = 1, seed = 9032)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)

cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- compare(fit, ps = p_vector)
#                    A    B mean_v.conservative.false mean_v.conservative.true
# True            0.75 1.25                     1.500                    2.150
# 5% Estimate     0.47 0.92                     0.745                    1.596
# 50% Estimate    0.65 1.42                     1.467                    2.123
# 97.5% Estimate  0.90 2.02                     2.330                    2.845
# Median-True    -0.10 0.17                    -0.033                   -0.027
#                mean_v.liberal.false mean_v.liberal.true sd_v.true      t0
# True                           1.50               2.500     0.100  0.1500
# 5% Estimate                    0.33               1.837     0.081  0.0081
# 50% Estimate                   1.10               2.409     0.120  0.0710
# 97.5% Estimate                 2.01               3.217     0.192  0.2047
# Median-True                   -0.40              -0.091     0.020 -0.0790
