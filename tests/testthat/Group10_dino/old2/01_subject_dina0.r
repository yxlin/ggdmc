#!/usr/bin/env Rscript
q(save = "no")
cat("\n\n----------------- Testing model0 --------------------\n")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N100.rda")
save_path <- file.path(home_dir, "Group10_cdm_subjects/data/subject_dina0_N100.rda")

# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N1000.rda")
# save_path <- file.path(home_dir, "Group10_cdm_subjects/data/subject_dina0_N1000.rda")

# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N5000.rda")
# save_path <- file.path(home_dir, "Group10_cdm_subjects/data/subject_dina0_N5000.rda")

# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N10000.rda")
# save_path <- file.path(home_dir, "Group10_cdm_subjects/data/subject_dina0_N10000.rda")

# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N50000.rda")
# save_path <- file.path(home_dir, "Group10_cdm_subjects/data/subject_dina0_N50000.rda")


source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)
cat(paste0("Load data file from ", data_path), "\n")

ncore <- 1
# load(save_path)

# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     ncore = ncore,
#     sub_migration_prob = 0.05, thin = 2, is_pblocked = FALSE,
#     seed = 9032
# )
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    ncore = ncore, nmc = 2, nchain = 3,
    sub_migration_prob = 0.00, thin = 1, is_pblocked = FALSE,
    seed = 9032
)

fits <- fits0
p0 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = T, start = 1)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 4, is_pblocked = FALSE,
    seed = 9032
)

# fits2 <- ggdmc:::RestartSampling_subject(fits1,
#     sub_migration_prob = 0.00, thin = 2, is_pblocked = TRUE,
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


fits <- fits1

# cat(paste0("Save fit results to ", save_path), "\n")
# save(fits0, fits1, fits2, fits3, fits4, file = save_path)
# save(fits0, fits1, file = save_path)

fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

for (i in seq_len(ncore)) {
    hat <- gelman(fits[[i]])
    cat("Chain", i, "mpsrf = ", hat$mpsrf, "\n")
}


options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)

est_theta <- ggdmc::compare(fits[[1]], ps = p_vector)
est_theta <- ggdmc::compare(fits[[4]], ps = p_vector)

est_theta <- ggdmc::compare(fits[[2]], ps = p_vector)
# n = 100
#               guess1 guess2 guess3 slip1 slip2 slip3
# True           0.100  0.300   0.50  0.20 0.400 0.600
#    5 Estimate  0.069  0.071   0.44  0.30 0.089 0.046
#   50 Estimate  0.422  0.686   0.67  0.66 0.713 0.626
# 97.5 Estimate  0.823  1.069   0.86  1.06 1.077 1.071
# Median-True    0.322  0.386   0.17  0.46 0.313 0.026
# n = 1000
#    5 Estimate  0.047   0.24  0.420  0.13 0.205  0.110
#   50 Estimate  0.615   0.51  0.577  0.70 0.478  0.574
# 97.5 Estimate  0.946   0.88  0.719  1.08 0.846  1.022
# Median-True    0.515   0.21  0.077  0.50 0.078 -0.026
# n = 5000
# True           0.100 0.3000  0.500  0.20  0.400  0.60
#    5 Estimate  0.042 0.0663  0.462  0.11  0.047  0.18
#   50 Estimate  0.275 0.3087  0.519  0.34  0.289  0.35
# 97.5 Estimate  0.642 0.7107  0.740  0.71  0.691  1.07
# Median-True    0.175 0.0087  0.019  0.14 -0.111 -0.25
# n = 10000
#    5 Estimate   0.13   0.23 0.4203  0.21  0.213  0.071
#   50 Estimate   0.33   0.40 0.5014  0.42  0.384  0.311
# 97.5 Estimate   0.64   0.66 0.7352  0.72  0.644  1.073
# Median-True     0.23   0.10 0.0014  0.22 -0.016 -0.289
# n = 50000
# True           0.100   0.30   0.50  0.20  0.40  0.60
#    5 Estimate  0.077   0.30   0.47  0.16  0.28  0.21
#   50 Estimate  0.559   0.56   0.73  0.64  0.54  1.00
# 97.5 Estimate  0.579   0.57   0.74  0.66  0.55  1.09
# Median-True    0.459   0.26   0.23  0.44  0.14  0.40
# chain 2 and 4
#    5 Estimate   0.04  0.285  0.460 0.125  0.264  0.18
#   50 Estimate   0.14  0.335  0.488 0.225  0.313  0.26
# 97.5 Estimate   0.25  0.390  0.506 0.334  0.368  0.31
# Median-True     0.04  0.035 -0.012 0.025 -0.087 -0.34
# pdf("Rplot.pdf")
p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
p0 <- ggdmc::plot(fit,
    pll = FALSE, den = TRUE, start = fit@nmc * 0.5,
    hide_legend = FALSE
)


# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = 1)
# p0 <- ggdmc::plot(fit, pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[1]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[2]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[3]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[4]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[5]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[6]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[7]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[8]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[9]], pll = TRUE, den = FALSE, start = 1)
p0 <- ggdmc::plot(fits[[10]], pll = TRUE, den = FALSE, start = 1)

p0 <- ggdmc::plot(fits[[1]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[2]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[3]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[4]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[5]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[6]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[7]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[8]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[9]], pll = TRUE, den = FALSE, start = 51)
p0 <- ggdmc::plot(fits[[10]], pll = TRUE, den = FALSE, start = 91)

p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, start = fits[[1]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, start = fits[[2]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, start = fits[[3]]@nmc * 0.5)

p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, den = TRUE, start = fits[[1]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, den = TRUE, start = fits[[2]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, den = TRUE, start = fits[[3]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[4]], facet_chains = F, pll = F, den = TRUE, start = fits[[2]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[5]], facet_chains = F, pll = F, den = TRUE, start = fits[[3]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[6]], facet_chains = F, pll = F, den = TRUE, start = fits[[1]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[7]], facet_chains = F, pll = F, den = TRUE, start = fits[[2]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[8]], facet_chains = F, pll = F, den = TRUE, start = fits[[3]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[9]], facet_chains = F, pll = F, den = TRUE, start = fits[[2]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[10]], facet_chains = F, pll = F, den = TRUE, start = fits[[3]]@nmc * 0.5)

# p0 <- ggdmc::plot(fits[[1]], pll = F, den = T, start = 1)
# p0 <- ggdmc::plot(fits[[2]], pll = F, den = T, start = 1)
# p0 <- ggdmc::plot(fits[[3]], pll = F, den = T, start = 51)
# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = 1)


# dev.off()

# n = 100
# mpsrf =  1.008118
# mpsrf =  1.043129
# mpsrf =  1.070978
# mpsrf =  1.054143
#               guess1 guess2 guess3    pi1    pi2    pi3 slip1 slip2 slip3
# True           0.100  0.300  0.500  0.250  0.250  0.250  0.20  0.40  0.60
#    5 Estimate  0.084  0.098  0.391  0.019  0.020  0.022  0.18  0.15  0.08
#   50 Estimate  0.395  0.440  0.565  0.212  0.226  0.224  0.63  0.58  0.47
# 97.5 Estimate  0.984  1.001  0.821  0.705  0.694  0.687  1.06  1.05  1.05
# Median-True    0.295  0.140  0.065 -0.038 -0.024 -0.026  0.43  0.18 -0.13

# N = 1000
# mpsrf =  1.035788
# mpsrf =  1.467893
# mpsrf =  1.569553
# mpsrf =  1.389407
#               guess1 guess2 guess3    pi1    pi2    pi3 slip1 slip2 slip3
# True            0.10   0.30   0.50  0.250  0.250  0.250  0.20  0.40  0.60
#    5 Estimate   0.52   0.51   0.60  0.026  0.023  0.025  0.60  0.52  0.57
#   50 Estimate   0.69   0.61   0.71  0.183  0.189  0.187  0.72  0.60  0.80
# 97.5 Estimate   1.07   1.05   1.06  0.569  0.499  0.574  1.05  0.82  1.08
# Median-True     0.59   0.31   0.21 -0.067 -0.061 -0.063  0.52  0.20  0.20

# N = 5000
# mpsrf =  4.854002
# mpsrf =  1.772545
# mpsrf =  1.961914
# mpsrf =  1.773462
#               guess1 guess2 guess3   pi1    pi2    pi3 slip1  slip2  slip3
# True           0.100   0.30  0.500 0.250  0.250  0.250 0.200  0.400  0.600
#    5 Estimate  0.062   0.12  0.460 0.043  0.023  0.022 0.073  0.052  0.037
#   50 Estimate  0.360   0.43  0.524 0.341  0.217  0.156 0.355  0.316  0.273
# 97.5 Estimate  0.918   0.94  0.945 0.639  0.513  0.530 0.647  0.622  0.838
# Median-True    0.260   0.13  0.024 0.091 -0.033 -0.094 0.155 -0.084 -0.327


# n = 50
#                guess1 guess2 guess3   pi1    pi2    pi3 slip1 slip2  slip3
# True            0.10   0.30   0.50  0.250  0.250  0.250  0.20  0.40  0.600
#    5 Estimate   0.12   0.11   0.31  0.022  0.027  0.021  0.11  0.12  0.060
#   50 Estimate   0.51   0.49   0.51  0.225  0.231  0.215  0.55  0.53  0.551
# 97.5 Estimate   1.02   1.02   0.81  0.642  0.666  0.664  1.03  1.04  1.064
# Median-True     0.41   0.19   0.01 -0.025 -0.019 -0.035  0.35  0.13 -0.049
# n = 500
# True            0.10   0.30   0.50  0.250  0.250 0.250  0.20  0.40  0.60
#    5 Estimate   0.51   0.47   0.54  0.012  0.035 0.050  0.53  0.42  0.51
#   50 Estimate   0.68   0.65   0.64  0.154  0.206 0.283  0.70  0.55  0.77
# 97.5 Estimate   1.06   1.07   0.98  0.543  0.577 0.665  1.07  0.88  1.08
# Median-True     0.58   0.35   0.14 -0.096 -0.044 0.033  0.50  0.15  0.17
# n = 1000
# True            0.10   0.30   0.50  0.250  0.250  0.250  0.20  0.40  0.600
#    5 Estimate   0.11   0.24   0.40  0.026  0.020  0.026  0.17  0.21  0.072
#   50 Estimate   0.54   0.59   0.60  0.193  0.199  0.218  0.55  0.48  0.535
# 97.5 Estimate   1.00   1.01   0.91  0.544  0.572  0.593  0.95  0.78  1.013
# Median-True     0.44   0.29   0.10 -0.057 -0.051 -0.032  0.35  0.08 -0.065
# n = 2000
#                guess1 guess2 guess3   pi1    pi2  pi3  slip1  slip2  slip3
# True           0.100  0.300  0.500  0.250  0.250 0.25  0.200  0.400  0.600
#    5 Estimate  0.013  0.027  0.383  0.015  0.105 0.21  0.015  0.139  0.013
#   50 Estimate  0.169  0.250  0.477  0.159  0.216 0.36  0.151  0.323  0.134
# 97.5 Estimate  0.345  0.426  0.521  0.343  0.391 0.58  0.376  0.442  0.306
# Median-True    0.069 -0.050 -0.023 -0.091 -0.034 0.11 -0.049 -0.077 -0.466
# n = 5000
#   guess1 guess2  guess3    pi1     pi2   pi3 slip1 slip2  slip3
# True            0.10   0.30  0.5000  0.250  0.2500 0.250 0.200 0.400  0.600
#    5 Estimate   0.15   0.10  0.3765  0.028  0.0081 0.031 0.044 0.277  0.043
#   50 Estimate   0.36   0.42  0.4917  0.220  0.0844 0.292 0.314 0.422  0.277
# 97.5 Estimate   1.02   0.93  1.0600  0.516  0.3090 0.593 0.704 0.541  0.819
# Median-True     0.26   0.12 -0.0083 -0.030 -0.1656 0.042 0.114 0.022 -0.323





# options(digits = 2)
#  est_theta <- ggdmc::compare(fits[[1]], ps = p_vector)
#  est_theta <- ggdmc::compare(fits[[2]], ps = p_vector)
#  est_theta <- ggdmc::compare(fits[[3]], ps = p_vector)
