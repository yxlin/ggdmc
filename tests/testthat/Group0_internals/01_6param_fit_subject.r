q(save = "no")
cat("\n\n-------------------- Run 6 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood", "lbaModel")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "M", st0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "lba"
)

p_vector <- c(
    A = .75, B = 1.25, mean_v.false = 1.5, mean_v.true = 2.5,
    sd_v.true = 0.1, t0 = 0.15
)

sub_model <- setLBA(model)
dat <- simulate(sub_model, nsim = 256, parameter_vector = p_vector, n_subject = 1)


sub_dmis <- ggdmcModel::BuildDMI(dat, model)

p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

sub_priors <- set_priors(p_prior = p_prior)


fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors, sub_migration_prob = 0.03, thin = 2, seed = 9032, is_pblocked = T)

fits1 <- RestartSampling_subject(fits0, sub_migration_prob = 0.00, thin = 2, seed = 9032)
fits <- fits1
fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- compare(fit, ps = p_vector)
#     A    B mean_v.false mean_v.true sd_v.true     t0
# True            0.75 1.25         1.50        2.50    0.1000  0.150
# 5% Estimate     0.46 0.84         0.29        1.64    0.0646  0.006
# 50% Estimate    0.65 1.34         1.18        2.27    0.1022  0.062
# 97.5% Estimate  0.97 2.09         2.28        3.28    0.1669  0.187
# Median-True    -0.10 0.09        -0.32       -0.23    0.0022 -0.088

p0 <- plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
p1 <- plot(fits[[1]], start = fits[[1]]@nmc * 0.5)
p1 <- plot(fits[[2]], start = fits[[1]]@nmc * 0.5)
p1 <- plot(fits[[3]], start = fits[[1]]@nmc * 0.5)
p1 <- plot(fits[[1]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
p1 <- plot(fits[[2]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
p1 <- plot(fits[[3]], pll = F, den = T, start = fits[[1]]@nmc * 0.5)
