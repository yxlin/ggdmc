#!/usr/bin/env Rscript
# q(save = "no")
cat("\n\n----------------- Testing model0 --------------------\n")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"

data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N000100.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N000200.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N000500.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N001000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N002000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N005000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N010000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N020000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N050000.rda")

load(data_path)
cat(paste0("Load data file from ", data_path), "\n")
sub_dmis[[1]]@model

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 1, is_pblocked = FALSE,
    seed = 9032
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 4, is_pblocked = FALSE,
    seed = 9032
)

fits <- fits1

fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

for (i in seq_len(3)) {
    hat <- gelman(fits[[i]])
    cat("Chain", i, "mpsrf = ", hat$mpsrf, "\n")
}


options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
p0 <- ggdmc::plot(fits0[[1]], start = 1)
p0 <- ggdmc::plot(fits0[[2]], start = 1)
p0 <- ggdmc::plot(fits0[[3]], start = 1)

p0 <- ggdmc::plot(fits0[[1]], start = 200)
p0 <- ggdmc::plot(fits0[[2]], start = 200)
p0 <- ggdmc::plot(fits0[[3]], start = 30)
