q(save = "no")
cat("\n\n-------------------- 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "lbaModel", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data0.rda"
load(fn)

fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
    sub_migration_prob = 0.05, thin = 4
)

fits1 <- RestartSampling_hyper(fits0, sub_migration_prob = 0.00, thin = 2)
# fits2 <- RestartSampling_hyper(fits1)

fits <- fits1
fit <- RebuildPosterior(fits)
options(digits = 2)
est <- ggdmcPhi::compare(fit, ps = true_vector)

hat <- gelman(fit, verbose = TRUE)
hat <- gelman(fits[[1]], verbose = TRUE)
hat <- gelman(fits[[2]], verbose = TRUE)
hat <- gelman(fits[[3]], verbose = TRUE)

# save(fits0, fits1, fits2, est, file = save_fn)

#  loc_A   loc_B loc_mean_v.false loc_mean_v.true loc_t0  sca_A
# True           0.4000  0.5000           0.1500           2.500  0.300  0.100
# 5% Estimate    0.3813  0.4640           0.0049           2.451  0.290  0.070
# 50% Estimate   0.4071  0.4957           0.0530           2.519  0.306  0.086
# 97.5% Estimate 0.4392  0.5337           0.1670           2.598  0.326  0.114
# Median-True    0.0071 -0.0043          -0.0970           0.019  0.006 -0.014
#                 sca_B sca_mean_v.false sca_mean_v.true sca_t0
# True           0.1000            0.200           0.200 0.0500
# 5% Estimate    0.0877            0.192           0.187 0.0457
# 50% Estimate   0.1069            0.247           0.229 0.0559
# 97.5% Estimate 0.1417            0.341           0.306 0.0740
# Median-True    0.0069            0.047           0.029 0.0059

# p0 <- ggdmc::plot(fit)
# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE)
# p0 <- ggdmc::plot(fit, pll = FALSE, den = FALSE)

# p0 <- ggdmcPhi::plot(fit, facet_chain = TRUE, hide_legend = FALSE)
# p0 <- ggdmcPhi::plot(fit, facet_chain = FALSE, hide_legend = FALSE)

# p0 <- ggdmcPhi::plot(fits[[1]], facet_chain = FALSE)
# p0 <- ggdmcPhi::plot(fits[[1]], facet_chain = FALSE, subchain = TRUE, chains = 1:3)
# p0 <- ggdmcPhi::plot(fits[[1]], facet_chain = TRUE, subchain = TRUE, chains = 1:3)

# p0 <- ggdmcPhi::plot(fits[[1]], pll = F, den = T)
# p0 <- ggdmcPhi::plot(fits[[1]], pll = F, den = T, subchain = TRUE, chains = 1:3)

# p0 <- ggdmcPhi::plot(fits[[1]], pll = F, den = F)
# p0 <- ggdmcPhi::plot(fits[[1]], pll = F, den = F, subchain = TRUE, chains = 1:3)

# p0 <- ggdmcPhi::plot(fits[[1]], pll = F, den = F, hide_legend = F)
# p0 <- ggdmcPhi::plot(fits[[1]], pll = F, den = F, subchain = TRUE, chains = 1:3)
