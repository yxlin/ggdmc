q(save = "no")
cat("\n\n-------------------- DDM a model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
fn <- "~/Documents/ggdmc/tests/testthat/Group6/data/ddm_data3.rda"
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/"
save_path <- paste0(wkdir, "fit_data/3_a_model.rda")
load(fn)
# load(save_path)

# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.08,
#     rp = 0.1, gamma_precursor = 5.38,
#     thin = 8L, pop_debug = F, seed = 9032
# )

# save(fits0, file = save_path)

# fits1 <- RestartSampling(fits0,
#     sub_migration_prob = 0.03,
#     pop_migration_prob = 0.05,
#     thin = 8L, seed = 9032
# )

# save(fits0, fits1, file = save_path)

# fits2 <- RestartSampling(fits1,
#     pop_migration_prob = 0.05,
#     sub_migration_prob = 0.05,
#     thin = 2L, seed = 9032
# )
# save(fits0, fits1, fits2, file = save_path)

# fits3 <- RestartSampling(fits2,
#     pop_migration_prob = 0.05,
#     sub_migration_prob = 0.05,
#     rp = 0.1, gamma_precursor = 5.38,
#     thin = 8L, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, file = save_path)


fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)
est_phi <- compare(phi, ps = true_vector)
est_theta <- compare_many(thetas, ps = ps)
#        a.s1.保守派 a.s1.自由派 a.s2.保守派 a.s2.自由派      sz     t0       v
# Mean      1.0306     2.50210       1.516      3.5213 0.01215 0.1452  2.1846
# True      1.0156     2.48513       1.501      3.5048 0.24711 0.1496  2.1464
# Diff     -0.0150    -0.01698      -0.015     -0.0164 0.23496 0.0044 -0.0382
# Sd        0.0602     0.00041       0.040      0.0038 0.00087 0.0142  0.3542
# True      0.0591     0.04655       0.049      0.0793 0.01686 0.0152  0.3509
# Diff     -0.0011     0.04614       0.009      0.0755 0.01599 0.0010 -0.0034
#            z
# Mean 0.38356
# True 0.39073
# Diff 0.00716
# Sd   0.00052
# True 0.03067
# Diff 0.03015
rhat <- gelman(phi)
rhat
# est_phi <- compare(fits[[1]]$phi, ps = true_vector)
# est_phi <- compare(fits[[2]]$phi, ps = true_vector)

# DT <- ggdmc::prepare_thetas_data(thetas)

# DT <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
# DT <- ggdmc::prepare_thetas_data(fits[[2]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
# DT <- ggdmc::prepare_thetas_data(fits[[3]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)

# p1 <- plot_thetas(DT)


# p3 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
# p3 <- plot(fits[[2]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
# p3 <- plot(fits[[3]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)

p3 <- plot(phi, pll = F, den = T)
