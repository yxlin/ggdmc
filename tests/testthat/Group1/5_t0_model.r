#  q(save = "no")
cat("\n\n--------------------Testing t0 Model--------------------")
rm(list = ls())
pkg <- c("ggdmcDE", "ggdmcModel", "ggdmcPrior", "ggdmcPhi", "ggdmcLikelihood", "ggplot2")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "tests/testthat/Group1/data/gdmc_data5.rda"
# fn <- "data/gdmc_data5.rda"
load(fn)

model <- ggdmcModel:::BuildModel(
    p_map = list(
        A = "1", B = "1",
        t0 = c("HANDEDNESS", "MOTOR_SKILL"),
        mean_v = c("S", "M"), sd_v = "M", st0 = "1"
    ),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(
        S = c("紅", "黃", "藍", "綠"),
        HANDEDNESS = c("left_hander", "right_hander"),
        MOTOR_SKILL = c("excellent", "good", "poor")
    ),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "lba"
)


dmis <- ggdmcModel::BuildDMI(dat, model)
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
priors <- set_priors(p_prior = p_prior)

nmc <- 1000
nchain <- model@npar * 3
thin <- 4
theta_input <- setThetaInput(nmc = nmc, nchain = nchain, pnames = model@pnames, thin = thin)
samples0 <- initialise_theta(theta_input, priors, dmis[[1]], seed = 929726)

# Burn-in ----------------------
de_input <- ggdmcDE::setDEInput(
    migration_prob_Hu = 0.05,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)
configs <- ggdmcDE::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]

fit0 <- run_subject(cfg, dmis[[1]], samples0, seed = 123)

# gdmc:::plot(fit0, start = fit0@nmc * 0.5)
# gdmc:::plot(fit1)



# Sampling ----------------------
de_input <- ggdmcDE::setDEInput(
    migration_prob_Hu = 0.00,
    nparameter = as.integer(theta_input@nparameter), nchain = as.integer(theta_input@nchain)
)

configs <- ggdmcDE::set_configs(prior = priors, theta_input = theta_input, de_input = de_input)
cfg <- configs[[1]]
fit1 <- run_subject(cfg, dmis[[1]], fit0, seed = 123)


# Cairo::CairoPDF("Rplot.pdf", family = "Noto Sans CJK TC")
# p3 <- gdmc:::plot(fit4)
# p4 <- gdmc:::plot(fit4, pll = FALSE, den = TRUE)
# dev.off()

est0 <- gdmc:::summary(fit4, start = 1, recovery = TRUE, ps = p_vector[model@pnames], verbose = TRUE)
#   A       B mean_v.紅.false mean_v.紅.true mean_v.綠.false
# True            0.2500  1.3300          1.2400         2.5000          1.7200
# 2.5% Estimate   0.1240  1.1036          0.8190         2.2162          1.4088
# 50% Estimate    0.2483  1.3291          1.1983         2.5166          1.7911
# 97.5% Estimate  0.3215  1.4942          1.5234         2.7813          2.1327
# Median-True    -0.0017 -0.0009         -0.0417         0.0166          0.0711
#                mean_v.綠.true mean_v.藍.false mean_v.藍.true mean_v.黃.false
# True                   3.1000          1.4400         2.9000          1.5200
# 2.5% Estimate          2.7978          1.0336         2.5973          1.1302
# 50% Estimate           3.1361          1.4144         2.9261          1.4957
# 97.5% Estimate         3.4561          1.7535         3.2221          1.8182
# Median-True            0.0361         -0.0256         0.0261         -0.0243
#                mean_v.黃.true sd_v.true t0.left_hander.excellent
# True                   2.7000    0.2000                   0.1100
# 2.5% Estimate          2.4058    0.1722                   0.1036
# 50% Estimate           2.7184    0.2038                   0.1189
# 97.5% Estimate         3.0042    0.2414                   0.1536
# Median-True            0.0184    0.0038                   0.0089
#                t0.left_hander.good t0.left_hander.poor
# True                        0.1200              0.1300
# 2.5% Estimate               0.1063              0.1179
# 50% Estimate                0.1222              0.1337
# 97.5% Estimate              0.1560              0.1678
# Median-True                 0.0022              0.0037
#                t0.right_hander.excellent t0.right_hander.good
# True                              0.0100               0.0200
# 2.5% Estimate                     0.0005               0.0077
# 50% Estimate                      0.0131               0.0236
# 97.5% Estimate                    0.0462               0.0567
# Median-True                       0.0031               0.0036
#                t0.right_hander.poor
# True                         0.0300
# 2.5% Estimate                0.0174
# 50% Estimate                 0.0330
# 97.5% Estimate               0.0664
# Median-True                  0.0030


# expected_values <- "tests/testthat/Group1/data/expected_model5.rda"
expected_values <- "data/expected_model5.rda"
samples <- samples0
est <- est0
save(samples, est, fit0, fit1, fit2, fit3, fit4, file = expected_values)

# load(expected_values)
# testthat::expect_true(all(samples0@theta[, , 1] == samples@theta[, , 1]))
# testthat::expect_true(all(est0 == est))
