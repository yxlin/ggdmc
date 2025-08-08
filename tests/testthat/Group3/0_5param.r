# q(save = "no")
cat("\n\n-------------------- 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data0.rda")
load(fn)

fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
    sub_migration_prob = 0.01, thin = 2, seed = 123
)

fits1 <- RestartSampling_hyper(fits0, sub_migration_prob = 0.00, thin = 1, seed = 123)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est <- ggdmc::compare(fit, ps = true_vector)

# mpsrf =  1.075406
#                 loc_A loc_B loc_mean_v.false loc_mean_v.true loc_t0  sca_A
# True           0.4000 0.500            0.150            2.50  0.300 0.1000
# 5% Estimate    0.3665 0.465            0.069            2.41  0.306 0.0869
# 50% Estimate   0.4093 0.504            0.241            2.49  0.316 0.1067
# 97.5% Estimate 0.4472 0.547            0.349            2.58  0.331 0.2577
# Median-True    0.0093 0.004            0.091           -0.01  0.016 0.0067
#                 sca_B sca_mean_v.false sca_mean_v.true sca_t0
# True           0.1000           0.2000           0.200  0.050
# 5% Estimate    0.0882           0.1468           0.186  0.029
# 50% Estimate   0.1091           0.2056           0.227  0.036
# 97.5% Estimate 0.1751           0.6052           0.344  0.048
# Median-True    0.0091           0.0056           0.027 -0.014
