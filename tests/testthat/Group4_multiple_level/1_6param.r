# q(save = "no")
cat("\n\n-------------------- 6 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

wkdir <- "~/Documents/ggdmc/tests/testthat/Group4/"
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/data/"

data_path <- paste0(data_dir, "lba_data1.rda")
save_path <- paste0(wkdir, "fit_data/6param.rda")
options(digits = 2)

load(data_path)
load(save_path)

fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.02,
    thin = 10L, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    pop_migration_prob = 0.01,
    sub_migration_prob = 0.00,
    is_hblocked = T,
    thin = 4L, seed = 9032
)
save(fits0, fits1, file = save_path)

# fits6 <- RestartSampling(fits5,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.00,
#     is_pblocked = F,
#     is_hblocked = T,
#     thin = 4L, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, fits4, fits5, fits6, file = save_path)


fits <- fits1
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

est_phi <- compare(phi, ps = true_vector)
#              loc_A loc_B loc_mean_v.false loc_mean_v.true loc_sd_v.true
# True            0.40 0.500            0.150            2.50        0.1000
# 5% Estimate     0.33 0.449            0.013            2.34        0.0058
# 50% Estimate    0.36 0.522            0.127            2.46        0.0518
# 97.5% Estimate  0.40 0.635            0.510            2.72        0.1130
# Median-True    -0.04 0.022           -0.023           -0.04       -0.0482
#                loc_t0  sca_A sca_B sca_mean_v.false sca_mean_v.true
# True            0.300  0.100 0.100            0.200           0.200
# 5% Estimate     0.267  0.063 0.090            0.022           0.029
# 50% Estimate    0.292  0.079 0.142            0.123           0.108
# 97.5% Estimate  0.322  0.108 0.218            0.331           0.213
# Median-True    -0.008 -0.021 0.042           -0.077          -0.092
#                sca_sd_v.true  sca_t0
# True                   0.100  0.0500
# 5% Estimate            0.078  0.0200
# 50% Estimate           0.113  0.0427
# 97.5% Estimate         0.164  0.0717
# Median-True            0.013 -0.0073

# est_theta <- compare_many(thetas, ps = ps)
# hat <- gelman(phi)
# cat("mpsrf = ", hat$mpsrf, "\n")
# hats <- lapply(thetas, gelman)
# sort(unlist(lapply(hats, function(x) x$mpsrf)), decreasing = TRUE)


# DT1 <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
# DT2 <- ggdmc::prepare_thetas_data(fits[[2]]$subject_theta, start = fits[[2]]$phi@nmc * 0.5)
# DT3 <- ggdmc::prepare_thetas_data(fits[[3]]$subject_theta, start = fits[[3]]$phi@nmc * 0.5)
# p1 <- plot_thetas(DT1)
# p1 <- plot_thetas(DT2)
# p1 <- plot_thetas(DT3)

# p2 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
# p2 <- plot(fits[[2]]$phi, facet_chains = F, start = fits[[2]]$phi@nmc * 0.5)
# p2 <- plot(fits[[3]]$phi, facet_chains = F, start = fits[[3]]$phi@nmc * 0.5)

# p2 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
# p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
