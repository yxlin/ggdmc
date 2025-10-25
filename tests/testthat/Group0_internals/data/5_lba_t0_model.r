# q(save = "no")
cat("\n\n-------------------- Generate t0 data --------------------")
rm(list = ls())
pkg <- c("lbaModel", "ggdmcPrior", "ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))


cat("\nWorking directory: ", getwd(), "\n")
wkdir <- "~/Documents/ggdmc/tests/testthat/Group1/data/"
helper_path <- paste0(wkdir, "helpers.r")
save_path <- paste0(wkdir, "lba_data5.rda")
source(helper_path)

hyper_model <- ggdmcModel:::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = c("S", "M"), sd_v = "M", st0 = "1", t0 = c("HANDEDNESS", "MOTOR_SKILL")),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(S = c("紅", "黃", "藍", "綠"), HANDEDNESS = c("left_hander", "right_hander"), MOTOR_SKILL = c("excellent", "good", "poor")),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "hyper"
)


model <- ggdmcModel:::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = c("S", "M"), sd_v = "M", st0 = "1", t0 = c("HANDEDNESS", "MOTOR_SKILL")),
    match_map = list(M = list(紅 = "反應東", 黃 = "反應南", 藍 = "反應西", 綠 = "反應北")),
    factors = list(S = c("紅", "黃", "藍", "綠"), HANDEDNESS = c("left_hander", "right_hander"), MOTOR_SKILL = c("excellent", "good", "poor")),
    constants = c(st0 = 0, sd_v.false = 1),
    accumulators = c("反應東", "反應南", "反應西", "反應北"),
    type = "lba"
)
pop_mean <- c(
    A = 0.25, B = 1.33,
    t0.left_hander.excellent = 0.11, t0.right_hander.excellent = 0.01,
    t0.left_hander.good = 0.12, t0.right_hander.good = 0.02,
    t0.left_hander.poor = 0.13, t0.right_hander.poor = 0.03,
    mean_v.紅.true = 2.5, mean_v.黃.true = 2.7,
    mean_v.藍.true = 2.9, mean_v.綠.true = 3.1,
    mean_v.紅.false = 1.24, mean_v.黃.false = 1.52,
    mean_v.藍.false = 1.44, mean_v.綠.false = 1.72,
    sd_v.true = 0.2
)

pop_scale <- c(
    A = 0.25, B = .2,
    t0.left_hander.excellent = 0.2, t0.right_hander.excellent = 0.2,
    t0.left_hander.good = 0.2, t0.right_hander.good = 0.2,
    t0.left_hander.poor = 0.2, t0.right_hander.poor = 0.2,
    mean_v.紅.true = .2, mean_v.黃.true = .2,
    mean_v.藍.true = .2, mean_v.綠.true = .3,
    mean_v.紅.false = .1, mean_v.黃.false = .1,
    mean_v.藍.false = .1, mean_v.綠.false = .1,
    sd_v.true = 0.2
)

pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(FALSE, model@npar)
)


p_vector <- c(
    A = 0.25, B = 1.33,
    t0.left_hander.excellent = 0.11, t0.right_hander.excellent = 0.01,
    t0.left_hander.good = 0.12, t0.right_hander.good = 0.02,
    t0.left_hander.poor = 0.13, t0.right_hander.poor = 0.03,
    mean_v.紅.true = 2.5, mean_v.黃.true = 2.7,
    mean_v.藍.true = 2.9, mean_v.綠.true = 3.1,
    mean_v.紅.false = 1.24, mean_v.黃.false = 1.52,
    mean_v.藍.false = 1.44, mean_v.綠.false = 1.72,
    sd_v.true = 0.2
)

# ---------------------------------------
sub_model <- lbaModel::setLBA(model)
pop_model <- lbaModel::setLBA(model, population_distribution = pop_dist)


dat <- simulate(sub_model, nsim = 1536, parameter_vector = p_vector, n_subject = 1)
hdat <- simulate(pop_model, nsim = 1536, n_subject = 32)

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

nmc <- 1000
thin <- 2
sub_theta_input <- setThetaInput(nmc = nmc, pnames = model@pnames, thin = thin, report_length = 200)

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

pop_theta_input <- setThetaInput(nmc = nmc, pnames = pop_priors@pnames, thin = thin, report_length = 200)

print(sub_theta_input)
print(pop_theta_input)

pop_samples <- initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671, verbose = FALSE)



save(hyper_model, model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis, hyper_dmi, sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_samples, pop_theta_input,
    file = save_path
)
