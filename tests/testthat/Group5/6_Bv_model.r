q(save = "no")
cat("\n\n-------------------- Bv model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmc/tests/testthat/Group1/data/lba_data6.rda"
save_path <- "~/Documents/ggdmc/tests/testthat/Group5/fit_data/6_Bv_model.rda"
options(digits = 2)

load(fn)
load(save_path)

# fits0 <- StartSampling(pop_dmis, pop_priors, sub_migration_prob = 0.06, thin = 8L, seed = 9032)
# save(fits0, file = save_path)

# fits1 <- RestartSampling(fits0,
#     sub_migration_prob = 0.00, pop_migration_prob = 0.05,
#     nmc = 1000L, report_length = 200L,
#     thin = 8L, seed = 9032
# )
# save(fits0, fits1, file = save_path)

# Check its size
obj_size_mb <- as.numeric(object.size(fits1)) / (1024^2)
print(paste("Size:", round(obj_size_mb, 2), "MB"))

fits <- fits1
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)
est_phi <- compare(phi, ps = true_vector)
est_theta <- compare_many(thetas, ps = ps)
#          A B.left.blue B.left.red B.right.blue B.right.red mean_v.high.false
# Mean  0.36        2.35       2.07         2.05        2.04             1.681
# True  0.26        2.58       2.28         2.27        2.25             1.756
# Diff -0.11        0.23       0.22         0.22        0.21             0.075
# Sd    0.48        0.52       0.42         0.43        0.44             0.069
# True  0.17        0.29       0.18         0.26        0.24             0.146
# Diff -0.31       -0.23      -0.24        -0.17       -0.20             0.078
#      mean_v.high.true mean_v.low.false mean_v.low.true mean_v.moderate.false
# Mean            1.960            1.228           4.177                 1.570
# True            2.016            1.474           4.100                 1.706
# Diff            0.056            0.246          -0.077                 0.137
# Sd              0.299            0.254           0.919                 0.139
# True            0.248            0.175           0.459                 0.154
# Diff           -0.051           -0.079          -0.460                 0.015
#      mean_v.moderate.true sd_v.true     t0
# Mean                 3.22    0.2635  0.097
# True                 3.08    0.2549  0.059
# Diff                -0.15   -0.0085 -0.037
# Sd                   1.47    0.1789  0.110
# True                 0.44    0.1669  0.042
# Diff                -1.03   -0.0120 -0.068

rhat <- gelman(phi)
plot(phi, pll = F, den = T)
rhat
