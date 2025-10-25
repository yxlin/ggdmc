# q(save = "no")
cat("\n\n---------------- DDM 6 parameter v model -------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group5_gen_ddm/data/ddm2.rda")
save_path <- file.path(home_dir, "Group6_ddm_subjects/data/v6param.rda")
# gplot_path <- file.path(home_dir, "Group6_ddm_subjects/gplot.r")
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
#                  a    sz     t0 v.s1 v.s2     z
# True           1.0 0.250 0.1500  2.5 2.10  0.38
# 5% Estimate    1.0 0.039 0.1353  2.8 1.95  0.30
# 50% Estimate   1.1 0.310 0.1512  3.5 2.55  0.37
# 97.5% Estimate 1.3 0.714 0.1623  4.5 3.31  0.49
# Median-True    0.1 0.060 0.0012  1.0 0.45 -0.01
