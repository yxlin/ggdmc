# q(save = "no")
cat("\n\n-------------------- 6 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))


home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group5_gen_ddm/data/ddm1.rda")
save_path <- file.path(home_dir, "Group6_ddm_subjects/data/6param.rda")
# gplot_path <- file.path(home_dir, "Group6_ddm_subjects/gplot.r")
load(data_path)


fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05,
    thin = 2, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 5,
    seed = 9032
)
save(fits0, fits1, file = save_path)

fits <- fits1
fit <- RebuildPosterior(fits)
options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
#                    a    sv     sz       t0    v       z
# True          1.0000 0.100  0.250  0.15000 2.50  0.3800
#    5 Estimate 0.9422 0.093  0.017  0.14153 2.20  0.3357
#   50 Estimate 1.0087 0.785  0.174  0.14917 2.63  0.3774
# 97.5 Estimate 1.1378 1.996  0.463  0.15628 3.45  0.4348
# Median-True   0.0087 0.685 -0.076 -0.00083 0.13 -0.0026
