# q(save = "no")
cat("\n\n-------------------- B model --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(data_dir, "data/lba_data3.rda")
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
mpsrf <- 1.051379
#                loc_A loc_B.紅.保守派 loc_B.紅.自由派 loc_B.綠.保守派
# True            0.25           2.400          1.2000           4.200
# 5% Estimate     0.23           2.348          1.1635           4.115
# 50% Estimate    0.27           2.422          1.1944           4.255
# 97.5% Estimate  0.32           2.507          1.2312           4.423
# Median-True     0.02           0.022         -0.0056           0.055
#                loc_B.綠.自由派 loc_B.藍.保守派 loc_B.藍.自由派 loc_B.黃.保守派
# True                    2.1000           3.600           1.800           3.000
# 5% Estimate             2.0538           3.470           1.778           2.851
# 50% Estimate            2.1091           3.551           1.812           2.962
# 97.5% Estimate          2.1830           3.647           1.849           3.093
# Median-True             0.0091          -0.049           0.012          -0.038
#                loc_B.黃.自由派 loc_mean_v.false loc_mean_v.true loc_sd_v.true
# True                     1.500           1.1500           2.800         0.800
# 5% Estimate              1.473           1.1276           2.790         0.754
# 50% Estimate             1.501           1.1557           2.845         0.789
# 97.5% Estimate           1.539           1.1898           2.901         0.836
# Median-True              0.001           0.0057           0.045        -0.011
#                loc_t0 sca_A sca_B.紅.保守派 sca_B.紅.自由派 sca_B.綠.保守派
# True           0.1000 0.100           0.200          0.1000            0.40
# 5% Estimate    0.0981 0.090           0.186          0.0772            0.37
# 50% Estimate   0.1021 0.112           0.228          0.0947            0.46
# 97.5% Estimate 0.1103 0.184           0.320          0.5435            0.64
# Median-True    0.0021 0.012           0.028         -0.0053            0.06
#                sca_B.綠.自由派 sca_B.藍.保守派 sca_B.藍.自由派 sca_B.黃.保守派
# True                     0.200           0.300          0.1000           0.300
# 5% Estimate              0.149           0.201          0.0821           0.280
# 50% Estimate             0.183           0.248          0.1009           0.344
# 97.5% Estimate           0.263           0.352          0.1664           0.464
# Median-True             -0.017          -0.052          0.0009           0.044
#                sca_B.黃.自由派 sca_mean_v.false sca_mean_v.true sca_sd_v.true
# True                    0.1000            0.100           0.200         0.100
# 5% Estimate             0.0742            0.070           0.128         0.091
# 50% Estimate            0.0906            0.086           0.157         0.112
# 97.5% Estimate          0.3542            0.127           0.231         0.229
# Median-True            -0.0094           -0.014          -0.043         0.012
#                sca_t0
# True           0.0100
# 5% Estimate    0.0093
# 50% Estimate   0.0114
# 97.5% Estimate 0.5571
# Median-True    0.0014
