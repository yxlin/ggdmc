# q(save = "no")
cat("\n\n--------------------t0 model (17 parameters)--------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/ggdmc_data5.rda"
load(fn)


fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.06, thin = 8, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.02, thin = 4, seed = 9032)
fits2 <- RestartSampling_subject(fits1, sub_migration_prob = 0.00, thin = 8, seed = 9032)

fits <- fits2
fit <- RebuildPosterior(fits)

options(digits = 2)
est_phi <- compare(fit, ps = p_vector)
