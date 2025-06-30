q(save = "no")
cat("\n\n-------------------- 6 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data1.rda")
load(fn)


fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.03, thin = 2, seed = 9032, is_pblocked = T)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 2, seed = 9032)
fits <- fits1
fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- compare(fit, ps = p_vector)
#                    A     B mean_v.false mean_v.true sd_v.true    t0
# True           0.750 1.250         1.50         2.5     0.100 0.150
# 5% Estimate    0.498 0.542         0.54         1.7     0.078 0.027
# 50% Estimate   0.799 1.274         1.75         2.7     0.117 0.172
# 97.5% Estimate 1.246 2.488         3.14         4.1     0.187 0.341
# Median-True    0.049 0.024         0.25         0.2     0.017 0.022

# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# ggdmc::plot(fits[[1]], start = fits[[1]]@nmc * 0.5)
# ggdmc::plot(fits[[2]], start = fits[[1]]@nmc * 0.5)
# ggdmc::plot(fits[[3]], start = fits[[1]]@nmc * 0.5)
# ggdmc::plot(fits[[1]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
# ggdmc::plot(fits[[2]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
# ggdmc::plot(fits[[3]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
