# q(save = "no")
cat("\n\n-------------------- Generate B model --------------------")
rm(list = ls())
cat("\nWorking directory: ", getwd(), "\n")
pkg <- c("lbaModel", "ggdmcPrior", "ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ggdmc/tests/testthat/Group1/data/"
helper_path <- paste0(wkdir, "helpers.r")
save_path <- paste0(wkdir, "lba_data3.rda")
source(helper_path)


hyper_model <- ggdmcModel:::BuildModel(
    p_map = list(A = "1", B = c("S", "政黨傾向"), mean_v = "M", sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(S = c("紅", "黃", "藍", "綠"), 政黨傾向 = c("自由派", "保守派")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "hyper",
    verbose = FALSE
)

model <- ggdmcModel:::BuildModel(
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
sub_model <- lbaModel::setLBA(model)
pop_model <- lbaModel::setLBA(model, population_distribution = pop_dist)

p_vector <- c(
    A = 0.25, B.紅.自由派 = 1.2, B.黃.自由派 = 1.5, B.藍.自由派 = 1.8, B.綠.自由派 = 2.1,
    B.紅.保守派 = 2.4, B.黃.保守派 = 3.0, B.藍.保守派 = 3.6, B.綠.保守派 = 4.2,
    mean_v.false = 1.15, mean_v.true = 2.80, sd_v.true = 0.8, t0 = 0.1
)
dat <- simulate(sub_model, nsim = 1024, parameter_vector = p_vector, n_subject = 1)
hdat <- simulate(pop_model, nsim = 1024, n_subject = 32)

sub_dmis <- ggdmcModel::BuildDMI(dat, model)
pop_dmis <- ggdmcModel::BuildDMI(hdat, model)
hyper_dmi <- ggdmcModel::BuildDMI(hdat, hyper_model)

res <- complex_get_accuracy(dat)
res <- complex_get_accuracy(hdat)
ps <- attr(hdat, "parameters")

true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)

# Generate subj samples -------------------------
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

nmc <- 500
sub_theta_input <- setThetaInput(nmc = nmc, pnames = model@pnames)
sub_samples <- initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]], seed = 483128, verbose = FALSE)

save(hyper_model, model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi, sub_priors, sub_samples, sub_theta_input,
    file = save_path
)


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


pop_theta_input <- setThetaInput(nmc = nmc, pnames = pop_priors@pnames)

print(sub_theta_input)
print(pop_theta_input)

pop_samples <- initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671, verbose = FALSE)


save(hyper_model, model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi, sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input,
    file = save_path
)
