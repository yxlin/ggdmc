#!/usr/bin/env Rscript
# q(save = "no")
cat("\n\n----------------- Testing model11 --------------------\n")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina11.rda")
save_path <- file.path(home_dir, "Group10_cdm_subjects/data/dina11.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)

cat(paste0("Load data file from ", data_path), "\n")

# StartSampling_subject
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.03, thin = 2, is_pblocked = TRUE,
    seed = 9032
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2,
    seed = 9032
)
fits <- fits1

cat(paste0("Save fit results to ", save_path), "\n")
save(fits0, fits1, file = save_path)

fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 3)
est_theta <- ggdmc::compare(fit, ps = p_vector)
