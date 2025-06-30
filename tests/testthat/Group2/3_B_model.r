# q(save = "no")
cat("\n\n---------------B Model (13 parameters)-----------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data3.rda")
load(fn)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.01, thin = 2, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.02, thin = 2, seed = 9032)

fits <- fits1
fit <- RebuildPosterior(fits)
hat <- gelman(fit)

cat("mpsrf = ", hat$mpsrf, "\n")
options(digits = 2)
est_theta <- compare(fit, ps = p_vector)
# $mpsrf
# [1] 1.021133

#                    A B.紅.保守派 B.紅.自由派 B.綠.保守派 B.綠.自由派
# True           0.250       2.400       1.200        4.20        2.10
# 5% Estimate    0.037       1.906       0.743        3.52        1.55
# 50% Estimate   0.281       2.359       1.129        4.09        1.99
# 97.5% Estimate 0.784       2.856       1.508        4.75        2.45
# Median-True    0.031      -0.041      -0.071       -0.11       -0.11
#                B.藍.保守派 B.藍.自由派 B.黃.保守派 B.黃.自由派 mean_v.false
# True                  3.60       1.800        3.00       1.500        1.150
# 5% Estimate           2.86       1.301        2.37       1.001        0.851
# 50% Estimate          3.38       1.725        2.86       1.404        1.127
# 97.5% Estimate        3.98       2.152        3.40       1.809        1.455
# Median-True          -0.22      -0.075       -0.14      -0.096       -0.023
#                mean_v.true sd_v.true    t0
# True                 2.800     0.800 0.100
# 5% Estimate          2.537     0.752 0.067
# 50% Estimate         2.771     0.814 0.120
# 97.5% Estimate       3.076     0.903 0.207
# Median-True         -0.029     0.014 0.020
