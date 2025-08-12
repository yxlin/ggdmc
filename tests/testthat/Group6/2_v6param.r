# q(save = "no")
cat("\n\n---------------- DDM 6 parameter v model -------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")
options(digits = 2)
fn <- "~/Documents/ggdmc/tests/testthat/Group6/data/ddm_data2.rda"
wkdir <- "~/Documents/ggdmc/tests/testthat/Group6/"
save_path <- paste0(wkdir, "fit_data/2_v6param.rda")

load(fn)
load(save_path)
# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.05,
#     thin = 4L, pop_debug = F, seed = 9032
# )

# save(fits0, file = save_path)

# fits1 <- RestartSampling(fits0,
#     pop_migration_prob = 0.02,
#     sub_migration_prob = 0.00,
#     thin = 4L, seed = 9032
# )

# save(fits0, fits1, file = save_path)

# fits2 <- RestartSampling(fits1,
#     pop_migration_prob = 0.01,
#     sub_migration_prob = 0.00,
#     thin = 2L, seed = 9032
# )
# save(fits0, fits1, fits2, file = save_path)

# fits3 <- RestartSampling(fits2,
#     pop_migration_prob = 0.00,
#     sub_migration_prob = 0.00,
#     thin = 2L, seed = 9032
# )

fits <- fits2
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

est_phi <- compare(phi, ps = true_vector)
#                loc_a loc_sz  loc_t0 loc_v.s1 loc_v.s2  loc_z   sca_a sca_sz
# True           1.000  0.250  0.1500    2.500     2.10 0.3800  0.0500  0.010
# 5% Estimate    0.994  0.012  0.1405    2.282     2.01 0.3688  0.0252  0.064
# 50% Estimate   1.019  0.124  0.1465    2.472     2.17 0.3842  0.0473  0.151
# 97.5% Estimate 1.051  0.341  0.1536    2.727     2.35 0.4077  0.0837  0.328
# Median-True    0.019 -0.126 -0.0035   -0.028     0.07 0.0042 -0.0027  0.141
#                 sca_t0 sca_v.s1 sca_v.s2 sca_z
# True            0.0200    0.500   0.2500 0.010
# 5% Estimate     0.0130    0.235   0.1024 0.013
# 50% Estimate    0.0169    0.435   0.2567 0.028
# 97.5% Estimate  0.0232    0.670   0.5810 0.051
# Median-True    -0.0031   -0.065   0.0067 0.018


# DT <- ggdmcPhi::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
# p1 <- plot_thetas(DT)


# p3 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
# p3 <- plot(fits[[2]]$phi, facet_chains = F)
# p3 <- plot(fits[[3]]$phi, facet_chains = F)


# p3 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
rhat <- gelman(phi)
rhat

# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     sub_migration_prob = 0.05,
#     thin = 2, seed = 9032
# )
# fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 2, seed = 9032)

# fits <- fits1
# fit <- RebuildPosterior(fits)

# # Hyper only
# fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
#     sub_migration_prob = 0.05, thin = 2
# )
# fits1 <- RestartSampling_hyper(fits0)

# fits <- fits1
# fit <- RebuildPosterior(fits)

# est_phi <- ggdmc::compare(fit, ps = true_vector)
#  loc_a  loc_sz  loc_t0 loc_v.s1 loc_v.s2  loc_z sca_a sca_sz
# True            1.0000 0.25000  0.1500    2.500   2.1000  0.380 0.050 0.0100
# 5% Estimate     0.9688 0.24545  0.1379    2.062   1.9122  0.331 0.043 0.0088
# 50% Estimate    0.9942 0.25049  0.1474    2.423   2.1095  0.379 0.056 0.0112
# 97.5% Estimate  1.0199 0.25610  0.1561    2.832   2.2207  1.342 0.082 0.0217
# Median-True    -0.0058 0.00049 -0.0026   -0.077   0.0095 -0.001 0.006 0.0012
#                  sca_t0 sca_v.s1 sca_v.s2  sca_z
# True            0.02000    0.500    0.250 0.0100
# 5% Estimate     0.01538    0.456    0.215 0.0089
# 50% Estimate    0.01942    0.582    0.271 0.0115
# 97.5% Estimate  0.03381    1.920    0.548 2.4745
# Median-True    -0.00058    0.082    0.021 0.0015
