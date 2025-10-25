# q(save = "no")
cat("\n\n-------------------- 5 DDM parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group5_gen_ddm/data/ddm0.rda")
load(data_path)
# gplot_path <- file.path(home_dir, "Group6_ddm_subjects/gplot.r")
# source(gplot_path)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.02,
    thin = 2, seed = 9032
)
fits1 <- RestartSampling_subject(fits0,
    sub_migration_prob = 0.00,
    thin = 4, seed = 9032
)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
#                    a    sz    t0     v     z
# True           1.000 0.250 0.150  2.50 0.380
# 5% Estimate    0.958 0.044 0.147  2.06 0.363
# 50% Estimate   1.025 0.394 0.161  2.37 0.415
# 97.5% Estimate 1.125 0.701 0.170  2.78 0.487
# Median-True    0.025 0.144 0.011 -0.13 0.035

# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F)
# p1 <- ggdmc::plot(fit, facet_chains = T, start = 20)

# p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = T, start = fit@nmc * 0.5)
# p0 <- ggdmc::plot(fit, pll = F, den = T)
# p0 <- ggdmc::plot(fit, pll = F, den = F)


# p1 <- ggdmc::plot(fits[[1]], facet_chains = F)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F)
