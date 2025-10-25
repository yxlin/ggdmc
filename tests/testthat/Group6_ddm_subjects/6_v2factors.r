#  q(save = "no")
cat("\n\n-------------------- DDM v x (S & NOISE) model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group5_gen_ddm/data/ddm_data6.rda")
save_path <- file.path(home_dir, "Group6_ddm_subjects/data/v2factor.rda")
gplot_path <- file.path(home_dir, "Group6_ddm_subjects/gplot.r")


load(data_path)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.04,
    thin = 4, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.05, thin = 2, seed = 9032)
fits2 <- RestartSampling_subject(fits1, sub_migration_prob = 0.05, thin = 16, seed = 9032)
save(fits0, fits1, file = save_path)

fits <- fits2
fit <- RebuildPosterior(fits)
options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
