# q(save = "no")
cat("\n\n-------------------- 6 parameters --------------------")
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
wkdir <- file.path(home_dir, "Group2_lba_subjects/")
data_path <- file.path(home_dir, "Group0_gen_lba/data/lba_data1.rda")
cat("\nWorking directory: ", getwd(), "\n")
load(data_path)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.03, thin = 2, seed = 9032, is_pblocked = T)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 2, seed = 9032)
fits <- fits1
fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- compare(fit, ps = p_vector)
#                   A    B mean_v.false mean_v.true sd_v.true      t0
# True          0.750 1.25        1.500       2.500    0.1000  0.1500
#    5 Estimate 0.536 0.89        0.449       1.774    0.0587  0.0047
#   50 Estimate 0.785 1.51        1.465       2.592    0.0986  0.0613
# 97.5 Estimate 1.209 2.45        2.787       3.885    0.1674  0.2088
# Median-True   0.035 0.26       -0.035       0.092   -0.0014 -0.0887

# p0 <- plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- plot(fits[[1]], start = fits[[1]]@nmc * 0.5)
# p1 <- plot(fits[[2]], start = fits[[1]]@nmc * 0.5)
# p1 <- plot(fits[[3]], start = fits[[1]]@nmc * 0.5)
# p1 <- plot(fits[[1]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
# p1 <- plot(fits[[2]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
# p1 <- plot(fits[[3]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
