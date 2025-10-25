#!/usr/bin/env Rscript
# q(save = "no")
cat("\n\n----------------- Testing model0 --------------------\n")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_mvn_N50.rda")
save_path <- file.path(home_dir, "Group10_cdm_subjects/data/subject_dina0_mvn_N50.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)
cat(paste0("Load data file from ", data_path), "\n")
# p_vector
# guess1 guess2 guess3 guess4 guess5  mean1  mean2  sigma  slip1  slip2  slip3
#    0.1    0.2    0.3    0.4    0.5    0.5    0.2    0.2    0.2    0.4    0.6
# slip4  slip5
#        0.8    0.9

ncore <- 1
# load(save_path)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    ncore = ncore,
    sub_migration_prob = 0.0, thin = 1, is_pblocked = TRUE,
    seed = 9032, sub_debug = F,
)

fits <- fits0
p0 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = T, start = 11)

# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = T, start = 31)

# save(fits0, file = save_path)

# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     ncore = ncore,
#     sub_migration_prob = 0.05, thin = 2, is_pblocked = FALSE,
#     seed = 9032
# )
# save(fits0, file = save_path)

# fits1 <- ggdmc:::RestartSampling_subject(fits0,
#     sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE,
#     seed = 9032
# )
# save(fits0, fits1, file = save_path)

# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     nmc = 100, ncore = 1,
#     sub_migration_prob = 0.05, thin = 1, is_pblocked = TRUE,
#     seed = 9032, sub_debug = FALSE
# )

# plot(fits0[[1]])
# plot(fits0[[1]], start = 15)
# slotNames(fits0[[1]])
# str(fits0[[1]]@log_likelihoods)

# head(fits0[[1]]@log_likelihoods)
# tail(fits0[[1]]@log_likelihoods)
# fits0[[1]]@log_likelihoods
# ylim <- c(min(fits0[[1]]@log_likelihoods), max(fits0[[1]]@log_likelihoods))
# ylim <- c(-300, 0)

# plot(fits0[[1]]@log_likelihoods[14, ], type = "l", ylim = ylim)
# lines(fits0[[1]]@log_likelihoods[38, ], type = "l")

# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     ncore = 4,
#     sub_migration_prob = 0.05, thin = 8, is_pblocked = TRUE,
#     seed = 9032
# )

# fits1 <- ggdmc:::RestartSampling_subject(fits0,
#     sub_migration_prob = 0.00, thin = 2, is_pblocked = FALSE,
#     seed = 9032
# )
# save(fits0, fits1, file = save_path)

# fits <- fits1
# cat(paste0("Save fit results to ", save_path), "\n")

# fit <- RebuildPosterior(fits)
# hat <- gelman(fit)
# cat("mpsrf = ", hat$mpsrf, "\n")

# hat <- gelman(fits[[1]])
# cat("mpsrf = ", hat$mpsrf, "\n")

# hat <- gelman(fits[[2]])
# cat("mpsrf = ", hat$mpsrf, "\n")

# hat <- gelman(fits[[3]])
# cat("mpsrf = ", hat$mpsrf, "\n")



# options(digits = 2)
# options(scipen = 999)
# est_theta <- ggdmc::compare(fits0[[1]], ps = p_vector)

# est_theta <- ggdmc::compare(fit, ps = p_vector)
# n = 50
#    guess1 guess2 guess3 guess4 guess5  mean1   mean2 sigma slip1
# True            0.10  0.200  0.300  0.400  0.500  0.500 0.20000 0.200 0.200
#    5 Estimate   0.07  0.054  0.087  0.009  0.083  0.050 0.00064 0.065 0.041
#   50 Estimate   0.54  0.499  0.523  0.452  0.566  0.492 0.50669 0.597 0.574
# 97.5 Estimate   1.07  1.073  1.087  1.077  1.078  1.093 1.07312 1.043 1.066
# Median-True     0.44  0.299  0.223  0.052  0.066 -0.008 0.30669 0.397 0.374
#               slip2  slip3  slip4  slip5
# True          0.400  0.600  0.800  0.900
#    5 Estimate 0.089  0.117  0.066  0.078
#   50 Estimate 0.608  0.572  0.482  0.539
# 97.5 Estimate 1.072  1.079  1.068  1.083
# Median-True   0.208 -0.028 -0.318 -0.361

# n = 500
# True           0.100  0.200  0.300  0.400  0.5000 0.5000  0.20 0.200  0.20
#    5 Estimate  0.068  0.055  0.087  0.018  0.0826 0.0408  0.07 0.062  0.05
#   50 Estimate  0.457  0.499  0.555  0.480  0.4982 0.5067  0.63 0.566  0.46
# 97.5 Estimate  1.075  1.059  1.092  1.074  1.0784 1.0866  1.07 1.048  1.07
# Median-True    0.357  0.299  0.255  0.080 -0.0018 0.0067  0.43 0.366  0.26
#               slip2  slip3  slip4  slip5
# True          0.400  0.600  0.800  0.900
#    5 Estimate 0.089  0.122  0.066  0.087
#   50 Estimate 0.539  0.551  0.524  0.580
# 97.5 Estimate 1.068  1.083  1.068  1.083
# Median-True   0.139 -0.049 -0.276 -0.320

# n = 5000
# True           0.100  0.200  0.300  0.400  0.5000 0.500  0.20 0.200 0.200 0.400
#    5 Estimate  0.068  0.055  0.087  0.039  0.0826 0.045  0.07 0.065 0.064 0.089
#   50 Estimate  0.457  0.499  0.541  0.488  0.4982 0.515  0.63 0.566 0.477 0.547
# 97.5 Estimate  1.075  1.059  1.068  1.074  1.0784 1.087  1.07 1.045 1.066 1.072
# Median-True    0.357  0.299  0.241  0.088 -0.0018 0.015  0.43 0.366 0.277 0.147
#                slip3  slip4  slip5
# True           0.600  0.800  0.900
#    5 Estimate  0.122  0.066  0.085
#   50 Estimate  0.539  0.524  0.578
# 97.5 Estimate  1.083  1.068  1.083
# Median-True   -0.061 -0.276 -0.322



#               guess1 guess2 guess3  guess4    guess5  mean1 mean2    sigma
# True           0.100  0.200 0.3000  0.4000  0.500000 0.5000 0.200  0.20000
#    5 Estimate  0.015  0.018 0.0097  0.0042  0.000055 0.0015 0.014  0.00011
#   50 Estimate  0.344  0.540 0.7010  0.3599  0.095015 0.6076 0.856  0.16373
# 97.5 Estimate  1.093  1.083 1.1000  1.0897  1.097774 1.0983 1.092  0.99221
# Median-True    0.244  0.340 0.4010 -0.0401 -0.404985 0.1076 0.656 -0.03627
#               slip1 slip2  slip3 slip4  slip5
# True          0.200 0.400  0.600  0.80  0.900
#    5 Estimate 0.089 0.023  0.051  0.11  0.026
#   50 Estimate 0.586 0.436  0.470  0.75  0.652
# 97.5 Estimate 1.081 1.048  0.947  1.10  1.094
# Median-True   0.386 0.036 -0.130 -0.05 -0.248
# est_theta <- ggdmc::compare(fits[[1]], ps = p_vector)
# est_theta <- ggdmc::compare(fits[[2]], ps = p_vector)
# est_theta <- ggdmc::compare(fits[[3]], ps = p_vector)
# est_theta <- ggdmc::compare(fits[[4]], ps = p_vector)

# fits2 <- ggdmc:::RestartSampling_subject(fits1,
#     sub_migration_prob = 0.05, thin = 8, is_pblocked = TRUE,
#     seed = 9032
# )
# fits3 <- ggdmc:::RestartSampling_subject(fits2,
#     sub_migration_prob = 0.00, thin = 2, is_pblocked = TRUE,
#     seed = 9032
# )

# fits4 <- ggdmc:::RestartSampling_subject(fits3,
#     sub_migration_prob = 0.00, thin = 8, is_pblocked = TRUE,
#     seed = 9032
# )




# save(fits0, fits1, fits2, fits3, fits4, file = save_path)




# cat("Number of chains:", fits[[1]]@nchain, "\n")
# cat("Any NA/Inf:", any(!is.finite(fits[[1]]@theta)), "\n")
# cat("Theta variance:", var(as.vector(fits[[1]]@theta)), "\n")
# summary(fits[[1]])



# options(digits = 2)
# est_theta <- ggdmc::compare(fit, ps = p_vector)

# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, start = 1)


# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = T, start = 1)

# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, start = fits[[1]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, start = fits[[2]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, start = fits[[3]]@nmc * 0.5)

# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, den = TRUE, start = fits[[1]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, den = TRUE, start = fits[[2]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, den = TRUE, start = fits[[3]]@nmc * 0.5)
