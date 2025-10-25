# q(save = "no")
cat("\n\n-------------------- 5 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
fn <- "~/Documents/ggdmc/tests/testthat/Group6/data/ddm_data0.rda"
load(fn)
gplot_path <- "~/Documents/ggdmc/tests/testthat/Group5/"
helper_path <- paste0(gplot_path, "gplot.r")
source(helper_path)

fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors, sub_migration_prob = 0.05, thin = 2)
fits1 <- RestartSampling_hyper(fits0)

fits <- fits1
fit <- RebuildPosterior(fits)

est_phi <- ggdmc::compare(fit, ps = true_vector)
#                  loc_a loc_sz loc_t0  loc_v    loc_z   sca_a sca_sz  sca_t0
# True            1.0000  0.250  0.150  2.500  0.38000  0.0500 0.0100  0.0200
# 5% Estimate     0.9837  0.246  0.147  2.344  0.37663  0.0395 0.0085  0.0155
# 50% Estimate    0.9976  0.249  0.152  2.483  0.37914  0.0482 0.0103  0.0187
# 97.5% Estimate  1.0143  0.253  0.159  2.670  0.38211  0.0637 0.0137  0.0247
# Median-True    -0.0024 -0.001  0.002 -0.017 -0.00086 -0.0018 0.0003 -0.0013
#                  sca_v   sca_z
# True            0.5000  0.0100
# 5% Estimate     0.4074  0.0071
# 50% Estimate    0.4967  0.0087
# 97.5% Estimate  0.6728  0.0113
# Median-True    -0.0033 -0.0013
