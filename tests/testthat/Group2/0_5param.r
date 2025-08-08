#   q(save = "no")
cat("\n\n-------------------- Testing model0 - 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data0.rda")
load(fn)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.02, seed = 9032
)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, seed = 9032)
fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
#                    A    B mean_v.false mean_v.true     t0
# True           0.750 1.25        1.500       2.500  0.150
# 5% Estimate    0.103 0.77        0.972       2.112  0.015
# 50% Estimate   0.797 1.43        1.534       2.581  0.118
# 97.5% Estimate 1.576 2.17        2.078       3.103  0.267
# Median-True    0.047 0.18        0.034       0.081 -0.032

# hat <- gelman(fit)
# hat$mpsrf

# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, pll = F, start = fit@nmc * 0.5)

# p1 <- ggdmcPhi::plot(fits[[1]], facet_chains = F)
# p1 <- ggdmcPhi::plot(fits[[2]], facet_chains = F)
# p1 <- ggdmcPhi::plot(fits[[3]], facet_chains = F)

# p1 <- ggdmcPhi::plot(fits[[1]], facet_chains = F, pll = F)
# p1 <- ggdmcPhi::plot(fits[[2]], facet_chains = F, pll = F)
# p1 <- ggdmcPhi::plot(fits[[3]], facet_chains = F, pll = F)
