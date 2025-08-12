# q(save = "no")
cat("\n\n-------------------- 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood", "lbaModel")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

hyper_model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = "M", sd_v = "1", st0 = "1", t0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(sd_v = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "hyper",
    verbose = FALSE
)

model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    accumulators = c("r1", "r2"),
    type = "lba"
)

pop_mean <- c(A = .4, B = .5, mean_v.false = 0.15, mean_v.true = 2.5, t0 = 0.3)
pop_scale <- c(A = .1, B = .1, mean_v.false = .2, mean_v.true = .2, t0 = 0.05)
pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = c(0, 0, 0, 0, 0),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

pop_model <- setLBA(model, population_distribution = pop_dist)
hdat <- simulate(pop_model, nsim = 128, n_subject = 32)

pop_dmis <- ggdmcModel::BuildDMI(hdat, model)
hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)


p0 <- runif(model@npar)
names(p0) <- model@pnames
model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = c(0, 0, 0, 0, 0),
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

true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)



fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
    sub_migration_prob = 0.01, thin = 2, seed = 123
)

fits1 <- RestartSampling_hyper(fits0, sub_migration_prob = 0.00, thin = 1, seed = 123)

fits <- fits1
fit <- RebuildPosterior(fits)

hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est <- ggdmc::compare(fit, ps = true_vector)

