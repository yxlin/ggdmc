#   q(save = "no")
cat("\n\n-------------------- Testing model0 - 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data0.rda"
load(fn)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.02,
    thin = 2, seed = 9032
)

fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 4, seed = 9032)
fits2 <- RestartSampling_subject(fits1, sub_migration_prob = 0.00, thin = 4, seed = 9032)

fits <- fits2
fit <- RebuildPosterior(fits)

class(fit)

hat <- gelman(fit, verbose = TRUE)
hat <- gelman(fits[[1]], verbose = TRUE)
hat <- gelman(fits[[2]], verbose = TRUE)
hat <- gelman(fits[[3]], verbose = TRUE)

hat
options(digits = 2)
est_phi <- ggdmc::compare(fit, ps = p_vector)
#                    A    B mean_v.false mean_v.true     t0
# True           0.750 1.25          1.5        2.50  0.150
# 5% Estimate    0.111 0.73          1.6        1.58  0.014
# 50% Estimate   0.817 1.37          2.1        2.08  0.112
# 97.5% Estimate 1.618 2.08          2.6        2.61  0.261
# Median-True    0.067 0.12          0.6       -0.42 -0.038

p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
p1 <- ggdmc::plot(fit, facet_chains = F, pll = F, start = fit@nmc * 0.5)

p1 <- ggdmcPhi::plot(fits[[1]], facet_chains = F)
p1 <- ggdmcPhi::plot(fits[[2]], facet_chains = F)
p1 <- ggdmcPhi::plot(fits[[3]], facet_chains = F)

p1 <- ggdmcPhi::plot(fits[[1]], facet_chains = F, pll = F)
p1 <- ggdmcPhi::plot(fits[[2]], facet_chains = F, pll = F)
p1 <- ggdmcPhi::plot(fits[[3]], facet_chains = F, pll = F)

# p1 <- ggdmcPhi::plot(fits[[1]], pll = F, den = T)
# p1 <- ggdmcPhi::plot(fits[[2]], pll = F, den = T)
# p1 <- ggdmcPhi::plot(fits[[3]], pll = F, den = T)

# p2 <- ggdmcPhi::plot(fits[[1]], start = fits[[1]]@nmc * 0.5, facet_chains = FALSE)
# p2 <- ggdmcPhi::plot(fits[[2]], start = fits[[2]]@nmc * 0.5, facet_chains = FALSE)
# p2 <- ggdmcPhi::plot(fits[[3]], start = fits[[3]]@nmc * 0.5, facet_chains = FALSE)
