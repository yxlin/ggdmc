#!/usr/bin/env Rscript
# q(save = "no")
cat("\n\n----------------- Testing model0 --------------------\n")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"

# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N000100.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N000200.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N000500.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N001000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N002000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N005000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N010000.rda")
# data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N020000.rda")
data_path <- file.path(home_dir, "Group9_gen_cdm/data/subject_dina0_N050000.rda")

load(data_path)
#               guess1 guess2 guess3 slip1 slip2 slip3
# True           0.100  0.200  0.300  0.20  0.40  0.60
#    5 Estimate  0.078  0.044  0.039  0.24  0.15  0.44
#   50 Estimate  0.496  0.443  0.227  0.46  0.47  0.73
# 97.5 Estimate  0.960  0.965  0.472  0.71  0.89  0.98
# Median-True    0.396  0.243 -0.073  0.26  0.07  0.13
#               guess1 guess2 guess3 slip1 slip2 slip3
#
# True           0.100   0.20  0.300  0.20 0.400  0.60
#    5 Estimate  0.061   0.35  0.086  0.17 0.421  0.19
#   50 Estimate  0.889   0.44  0.596  0.53 0.492  0.87
# 97.5 Estimate  0.995   0.56  0.705  0.61 0.577  0.99
# Median-True    0.789   0.24  0.296  0.33 0.092  0.27
# 5000
#               guess1 guess2 guess3  slip1 slip2 slip3
# True          0.1000  0.200  0.300  0.200  0.40  0.60
#    5 Estimate 0.0159  0.066  0.271  0.146  0.21  0.37
#   50 Estimate 0.1073  0.161  0.312  0.189  0.28  0.43
# 97.5 Estimate 0.2668  0.273  0.345  0.261  0.36  0.47
# Median-True   0.0073 -0.039  0.012 -0.011 -0.12 -0.17
#
# N = 20,000
# mpsrf =  29.54512
# Chain 1 mpsrf =  1.009577
# Chain 2 mpsrf =  1.01442
# Chain 3 mpsrf =  1.021254
#               guess1 guess2 guess3   slip1  slip2 slip3
# True           0.100  0.200 0.3000  0.2000  0.400  0.60
#    5 Estimate  0.049  0.158 0.2768  0.1528  0.283  0.41
#   50 Estimate  0.145  0.215 0.3025  0.1952  0.324  0.44
# 97.5 Estimate  0.892  0.650 0.7167  0.5278  0.640  0.99
# Median-True    0.045  0.015 0.0025 -0.0048 -0.076 -0.16

# N = 50,000
# mpsrf =  45.98172
# Chain 1 mpsrf =  1.008902
# Chain 2 mpsrf =  1.009017
# Chain 3 mpsrf =  1.023759
#               guess1  guess2 guess3  slip1  slip2 slip3
# True           0.100  0.2000  0.300  0.200  0.400  0.60
#    5 Estimate  0.022  0.1646  0.299  0.140  0.291  0.42
#   50 Estimate  0.080  0.1989  0.313  0.166  0.316  0.44
# 97.5 Estimate  0.822  0.5952  0.735  0.497  0.604  1.00
# Median-True   -0.020 -0.0011  0.013 -0.034 -0.084 -0.16


cat(paste0("Load data file from ", data_path), "\n")

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 8, is_pblocked = FALSE,
    seed = 9032
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 4, is_pblocked = FALSE,
    seed = 9032
)

fits <- fits1

fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

for (i in seq_len(3)) {
    hat <- gelman(fits[[i]])
    cat("Chain", i, "mpsrf = ", hat$mpsrf, "\n")
}

# plot(fit, pll = F, den = TRUE)
options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
# N = 100       guess1 guess2 guess3   pi1   pi2   pi3 slip1 slip2 slip3
# True            0.10  0.200  0.300 0.158 0.263 0.151  0.20 0.400 0.600
#    5 Estimate   0.11  0.086  0.066 0.020 0.035 0.017  0.14 0.082 0.131
#   50 Estimate   0.51  0.428  0.214 0.205 0.300 0.181  0.45 0.425 0.641
# 97.5 Estimate   0.92  0.917  0.378 0.568 0.716 0.668  0.85 0.900 0.973
# Median-True     0.41  0.228 -0.086 0.047 0.037 0.031  0.25 0.025 0.041
# N = 200
# True           0.100   0.20  0.300 0.158 0.263 0.151  0.20 0.400 0.6000
#    5 Estimate  0.091   0.21  0.220 0.022 0.035 0.031  0.15 0.114 0.1450
#   50 Estimate  0.473   0.52  0.347 0.211 0.279 0.229  0.49 0.475 0.6099
# 97.5 Estimate  0.884   0.90  0.480 0.678 0.652 0.625  0.93 0.895 0.9675
# Median-True    0.373   0.32  0.047 0.053 0.016 0.078  0.29 0.075 0.0099

#    5 Estimate  0.043   0.13  0.049 0.041  0.013 0.033 0.029 0.308  0.031
#   50 Estimate  0.254   0.40  0.220 0.246  0.130 0.232 0.216 0.435  0.269
# 97.5 Estimate  0.425   0.63  0.319 0.615  0.357 0.596 0.394 0.554  0.483
# Median-True    0.154   0.20 -0.080 0.088 -0.133 0.081 0.016 0.035 -0.331

#    5 Estimate  0.055   0.14  0.094 0.018  0.013 0.049 0.045 0.351  0.049
#   50 Estimate  0.363   0.44  0.278 0.178  0.129 0.338 0.294 0.491  0.371
# 97.5 Estimate  0.968   0.64  0.659 0.531  0.362 0.661 0.974 0.581  0.985
# Median-True    0.263   0.24 -0.022 0.020 -0.134 0.187 0.094 0.091 -0.229

#    5 Estimate  0.031  0.075  0.252 0.20  0.050 0.056  0.017  0.091  0.019
#   50 Estimate  0.254  0.301  0.323 0.37  0.206 0.207  0.134  0.222  0.187
# 97.5 Estimate  0.427  0.434  0.363 0.61  0.405 0.419  0.318  0.369  0.407
# Median-True    0.154  0.101  0.023 0.22 -0.057 0.056 -0.066 -0.178 -0.413

#               guess1 guess2 guess3  pi1    pi2   pi3   slip1 slip2  slip3
# True           0.100  0.200 0.3000 0.16  0.263 0.151 0.20000  0.40  0.600
#    5 Estimate  0.077  0.057 0.2536 0.21  0.022 0.061 0.08343  0.19  0.061
#   50 Estimate  0.309  0.301 0.3068 0.37  0.135 0.214 0.20037  0.28  0.260
# 97.5 Estimate  0.430  0.434 0.3444 0.67  0.338 0.428 0.33012  0.37  0.433
# Median-True    0.209  0.101 0.0068 0.22 -0.128 0.064 0.00037 -0.12 -0.340

#    guess1 guess2 guess3  pi1    pi2   pi3  slip1  slip2  slip3
# True           0.100  0.200  0.300 0.16  0.263 0.151  0.200  0.400  0.600
#    5 Estimate  0.077  0.093  0.278 0.21  0.032 0.077  0.026  0.257  0.087
#   50 Estimate  0.298  0.324  0.313 0.37  0.135 0.240  0.136  0.315  0.299
# 97.5 Estimate  0.429  0.441  0.338 0.69  0.312 0.442  0.258  0.385  0.450
# Median-True    0.198  0.124  0.013 0.22 -0.128 0.089 -0.064 -0.085 -0.301

#  guess1 guess2  guess3  pi1    pi2   pi3 slip1  slip2 slip3
# True            0.10   0.20  0.3000 0.16  0.263 0.151  0.20  0.400  0.60
#    5 Estimate   0.11   0.11  0.2637 0.20  0.075 0.093  0.12  0.279  0.11
#   50 Estimate   0.28   0.32  0.2901 0.33  0.186 0.225  0.19  0.319  0.27
# 97.5 Estimate   0.41   0.43  0.3153 0.50  0.343 0.422  0.26  0.360  0.42
# Median-True     0.18   0.12 -0.0099 0.17 -0.077 0.074 -0.01 -0.081 -0.33

# est_theta <- ggdmc::compare(fits[[1]], ps = p_vector)
# est_theta <- ggdmc::compare(fits[[2]], ps = p_vector)
# est_theta <- ggdmc::compare(fits[[3]], ps = p_vector)
