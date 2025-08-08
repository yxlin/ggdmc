# q(save = "no")
cat("\n\n-------------------- A model --------------------")

rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data4.rda")
load(fn)

fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
    sub_migration_prob = 0.00, thin = 1, seed = 123, is_pblocked = TRUE
)

fits1 <- RestartSampling_hyper(fits0, sub_migration_prob = 0.00, thin = 1, seed = 123)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est <- ggdmc::compare(fit, ps = true_vector)
mpsrf <- 1.089713
#                loc_A.紅.保守派 loc_A.紅.自由派 loc_A.綠.保守派 loc_A.綠.自由派
# True                      1.23          0.2500           1.890           0.550
# 5% Estimate               1.18          0.2216           1.871           0.507
# 50% Estimate              1.26          0.2555           1.953           0.535
# 97.5% Estimate            1.36          0.2981           2.043           0.573
# Median-True               0.03          0.0055           0.063          -0.015
#                loc_A.藍.保守派 loc_A.藍.自由派 loc_A.黃.保守派 loc_A.黃.自由派
# True                     1.670           0.450           1.450          0.3500
# 5% Estimate              1.637           0.410           1.471          0.3196
# 50% Estimate             1.696           0.451           1.536          0.3561
# 97.5% Estimate           1.769           0.498           1.622          0.4053
# Median-True              0.026           0.001           0.086          0.0061
#                 loc_B loc_mean_v.false loc_mean_v.true loc_sd_v.true  loc_t0
# True            1.250           1.1500          2.8000        0.8000  0.1000
# 5% Estimate     1.206           1.1255          2.7761        0.7603  0.0076
# 50% Estimate    1.231           1.1548          2.8081        0.7908  0.0601
# 97.5% Estimate  1.261           1.1917          2.8443        0.8299  0.1400
# Median-True    -0.019           0.0048          0.0081       -0.0092 -0.0399
#                sca_A.紅.保守派 sca_A.紅.自由派 sca_A.綠.保守派 sca_A.綠.自由派
# True                     0.200          0.1000           0.200           0.100
# 5% Estimate              0.205          0.0744           0.201           0.071
# 50% Estimate             0.249          0.0926           0.248           0.088
# 97.5% Estimate           0.361          0.3603           0.345           0.390
# Median-True              0.049         -0.0074           0.048          -0.012
#                sca_A.藍.保守派 sca_A.藍.自由派 sca_A.黃.保守派 sca_A.黃.自由派
# True                    0.2000           0.100           0.200           0.100
# 5% Estimate             0.1573           0.094           0.182           0.093
# 50% Estimate            0.1922           0.116           0.222           0.115
# 97.5% Estimate          0.2670           0.469           0.336           0.984
# Median-True            -0.0078           0.016           0.022           0.015
#                 sca_B sca_mean_v.false sca_mean_v.true sca_sd_v.true sca_t0
# True            0.100           0.1000          0.1000        0.1000  0.100
# 5% Estimate     0.063           0.0747          0.0757        0.0795  0.094
# 50% Estimate    0.076           0.0919          0.0924        0.0978  0.130
# 97.5% Estimate  0.121           0.1430          0.3627        0.5790  0.353
# Median-True    -0.024          -0.0081         -0.0076       -0.0022  0.030
