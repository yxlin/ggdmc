# q(save = "no")
## TODO: old test code to be updated
cat("\n\n--------------------Testing B Model--------------------")
rm(list = ls())
cat("\nWorking directory: ", getwd(), "\n")
pkg <- c("lbaModel", "ggdmcPrior", "ggdmc", "ggdmcModel")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))


model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = c("S", "政黨傾向"), mean_v = "M", sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(S = c("紅", "黃", "藍", "綠"), 政黨傾向 = c("自由派", "保守派")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "lba"
)

pop_mean <- c(
    A = 0.25, B.紅.自由派 = 1.2, B.黃.自由派 = 1.5, B.藍.自由派 = 1.8, B.綠.自由派 = 2.1,
    B.紅.保守派 = 2.4, B.黃.保守派 = 3.0, B.藍.保守派 = 3.6, B.綠.保守派 = 4.2,
    mean_v.true = 2.80, mean_v.false = 1.15, sd_v.true = 0.8, t0 = 0.1
)

pop_scale <- c(
    A = 0.1, B.紅.自由派 = .1, B.黃.自由派 = .1, B.藍.自由派 = .1, B.綠.自由派 = .2,
    B.紅.保守派 = .2, B.黃.保守派 = .3, B.藍.保守派 = .3, B.綠.保守派 = .4,
    mean_v.true = .2, mean_v.false = .1, sd_v.true = 0.1, t0 = 0.01
)

pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(FALSE, model@npar)
)


# ---------------------------------------
pop_model <- lbaModel::setLBA(model, population_distribution = pop_dist)

hdat <- simulate(pop_model, nsim = 1024, n_subject = 32)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)

true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)


# Generate pop samples -------------------------
p0 <- runif(model@npar)
names(p0) <- model@pnames

model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("tnorm", model@npar),
    log_p = rep(TRUE, model@npar)
)

# Prior log likelihoods
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
l_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
s_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

pop_priors <- ggdmcPrior::set_priors(p_prior = model_likelihood, l_prior = l_prior, s_prior = s_prior)

# The following steps will take a while.
fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.02,
    thin = 2L, seed = 9032
)

fits1 <- RestartSampling(fits0,
    pop_migration_prob = 0.01,
    sub_migration_prob = 0.01,
    thin = 2L, seed = 9032
)
# save(fits0, fits1, file = save_path)


#                      A    B.紅.保守派 B.紅.自由派   B.綠.保守派 B.綠.自由派
# True            0.2500      2.4000      1.2000      4.2000      2.1000
# 2.5% Estimate   0.0047      1.6774      0.6351      3.1318      1.4086
# 50% Estimate    0.1233      2.4144      1.2296      4.1593      2.0982
# 97.5% Estimate  0.6184      3.2029      1.7997      5.1885      2.8129
# Median-True    -0.1267      0.0144      0.0296     -0.0407     -0.0018
#                  B.藍.保守派   B.藍.自由派  B.黃.保守派  B.黃.自由派 mean_v.false
# True                3.6000      1.8000      3.0000      1.5000       1.1500
# 2.5% Estimate       2.6822      1.1314      2.2000      0.9273       0.6324
# 50% Estimate        3.5775      1.8189      3.0532      1.5497       1.1933
# 97.5% Estimate      4.5049      2.4884      3.9173      2.1907       1.7097
# Median-True        -0.0225      0.0189      0.0532      0.0497       0.0433
#                mean_v.true sd_v.true      t0
# True                2.8000    0.8000  0.1000
# 2.5% Estimate       2.3742    0.6575  0.0072
# 50% Estimate        2.8530    0.7684  0.0968
# 97.5% Estimate      3.3391    0.8635  0.2022
# Median-True         0.0530   -0.0316 -0.0032
