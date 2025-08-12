# q(save = "no")
cat("\n\n--------------------Testing Drift Rate Model--------------------")
rm(list = ls())

pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood", "lbaModel")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))


model <- BuildModel(
    p_map = list(A = "1", B = "1", mean_v = c("POLITICAL_VIEW", "M"), sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2"), POLITICAL_VIEW = c("liberal", "conservative")),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("r1", "r2"),
    type = "lba"
)

pop_mean <- c(
    A = .4, B = .5, mean_v.liberal.false = .15, mean_v.liberal.true = 2.5,
    mean_v.conservative.false = .15, mean_v.conservative.true = 2.15,
    sd_v.true = 0.1, t0 = 0.3
)

pop_scale <- c(
    A = .1, B = .1,
    mean_v.conservative.false = .2, mean_v.conservative.true = .2,
    mean_v.liberal.false = .2,
    mean_v.liberal.true = .2,
    sd_v.true = 0.1, t0 = 0.1
)

pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(FALSE, model@npar)
)


true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)


pop_model <- lbaModel::setLBA(model, population_distribution = pop_dist)
hdat <- simulate(pop_model, nsim = 256, n_subject = 32)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)


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


fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.02,
    thin = 16L, seed = 9032
)
# save(fits0, file = save_path)

# Has not converged yet
fits <- fits0
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

options(digits = 2)
est_phi <- compare(phi, ps = true_vector)
#    loc_A loc_B loc_mean_v.conservative.false
# True            0.40  0.50                          0.15
# 5% Estimate     0.08  0.11                          0.13
# 50% Estimate    0.51  0.68                          0.96
# 97.5% Estimate  7.56  5.98                          7.69
# Median-True     0.11  0.18                          0.81
#                loc_mean_v.conservative.true loc_mean_v.liberal.false
# True                                   2.15                     0.15
# 5% Estimate                            0.66                     0.15
# 50% Estimate                           2.96                     1.07
# 97.5% Estimate                         8.93                     7.72
# Median-True                            0.81                     0.92
#                loc_mean_v.liberal.true loc_sd_v.true loc_t0 sca_A sca_B
# True                              2.50         0.100  0.300  0.10  0.10
# 5% Estimate                       0.76         0.014  0.045  0.16  0.21
# 50% Estimate                      3.41         0.142  0.275  0.69  0.73
# 97.5% Estimate                    9.06         7.689  7.043  8.55  8.76
# Median-True                       0.91         0.042 -0.025  0.59  0.63
#                sca_mean_v.conservative.false sca_mean_v.conservative.true
# True                                    0.20                         0.20
# 5% Estimate                             0.33                         0.22
# 50% Estimate                            1.05                         2.44
# 97.5% Estimate                          8.78                         9.65
# Median-True                             0.85                         2.24
#                sca_mean_v.liberal.false sca_mean_v.liberal.true sca_sd_v.true
# True                               0.20                    0.20          0.10
# 5% Estimate                        0.35                    0.31          0.11
# 50% Estimate                       1.12                    3.05          0.32
# 97.5% Estimate                     8.82                    9.71          8.97
# Median-True                        0.92                    2.85          0.22
#                sca_t0
# True             0.10
# 5% Estimate      0.10
# 50% Estimate     0.28
# 97.5% Estimate   8.10
# Median-True      0.18

hat <- gelman(phi)
cat("mpsrf = ", hat$mpsrf, "\n")

plot(phi, pll = F, den = T)

DT <- prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.5)
p1 <- plot_thetas(DT)
p2 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
p2 <- plot(fits[[2]]$phi, facet_chains = F, start = fits[[2]]$phi@nmc * 0.5)
p2 <- plot(fits[[3]]$phi, facet_chains = F, start = fits[[3]]$phi@nmc * 0.5)
p2 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
p3 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)
