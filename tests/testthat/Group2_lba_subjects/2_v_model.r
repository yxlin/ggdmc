# q(save = "no")
cat("\n\n--------------------Drift Rate Model--------------------")
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat"
wkdir <- file.path(home_dir, "Group2_lba_subjects/")
data_path <- file.path(home_dir, "Group0_gen_lba/data/lba_data2.rda")
cat("\nWorking directory: ", getwd(), "\n")
load(data_path)


fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.01, thin = 2, seed = 9032)
fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.01, thin = 1, seed = 9032)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)

cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- compare(fit, ps = p_vector)
#    A    B mean_v.conservative.false mean_v.conservative.true
# True          0.75 1.25                      1.50                     2.15
#    5 Estimate 0.60 1.06                      1.29                     1.95
#   50 Estimate 0.87 1.83                      2.16                     2.74
# 97.5 Estimate 1.26 2.71                      3.23                     3.75
# Median-True   0.12 0.58                      0.66                     0.59
#               mean_v.liberal.false mean_v.liberal.true sd_v.true     t0
# True                          1.50                2.50     0.100  0.150
#    5 Estimate                 1.37                2.33     0.116  0.011
#   50 Estimate                 2.25                3.18     0.171  0.094
# 97.5 Estimate                 3.43                4.29     0.259  0.241
# Median-True                   0.75                0.68     0.071 -0.056

#                    A    B mean_v.conservative.false mean_v.conservative.true
# True            0.75 1.25                     1.500                    2.150
# 5% Estimate     0.47 0.92                     0.745                    1.596
# 50% Estimate    0.65 1.42                     1.467                    2.123
# 97.5% Estimate  0.90 2.02                     2.330                    2.845
# Median-True    -0.10 0.17                    -0.033                   -0.027
#                mean_v.liberal.false mean_v.liberal.true sd_v.true      t0
# True                           1.50               2.500     0.100  0.1500
# 5% Estimate                    0.33               1.837     0.081  0.0081
# 50% Estimate                   1.10               2.409     0.120  0.0710
# 97.5% Estimate                 2.01               3.217     0.192  0.2047
# Median-True                   -0.40              -0.091     0.020 -0.0790
