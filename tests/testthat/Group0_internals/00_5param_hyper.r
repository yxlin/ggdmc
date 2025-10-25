# q(save = "no")
cat("\n\n-------------------- 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood", "lbaModel")
ok <- sapply(pkg, require, character.only = TRUE)
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

# head(hdat)
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


# hyper_dmi
dmis <- pop_dmis
priors <- pop_priors
samples_list <- NULL
nmc <- 500L
nchain <- NULL
thin <- 1L
report_length <- 100L
max_init_attempts <- 1000L
is_print <- TRUE
sub_migration_prob <- 0.00
gamma_precursor <- 2.38
rp <- 0.001
is_pblocked <- FALSE
sub_debug <- FALSE
ncore <- 1L
seed <- 123

nparameter <- priors@nparameter # 10
pnames <- priors@pnames
nchain <- ggdmc:::.get_nchain(priors, nchain)

if (is.null(priors@h_prior)) {
    stop("h_prior has not found")
}

theta_input <- setThetaInput(
    nmc = nmc, nchain = nchain, thin = thin,
    pnames = pnames, report_length = report_length, max_init_attempts = max_init_attempts,
    is_print = is_print
)

de_input <- setDEInput(
    sub_migration_prob = sub_migration_prob,
    gamma_precursor = gamma_precursor, rp = rp,
    is_pblocked = is_pblocked,
    nparameter = as.integer(nparameter), nchain = as.integer(nchain),
    sub_debug = sub_debug
)
config_list <- set_configs(
    prior = priors, theta_input = theta_input, de_input = de_input,
    ncore = ncore, seed = seed
)

if (is.null(samples_list)) {
    message("Initialise ", ncore, " independent sets of new samples")
    samples_list <- lapply(seq_len(ncore), function(i) {
        tmp_list <- initialise_phi(theta_input, priors, dmis, seed = config_list[[i]]@seed)
        tmp_list$phi
    })
}

message("\nLaunching ", ncore, " independent chains")


res <- ggdmc:::parallel_lapply(
    ncore = ncore, config_list = config_list,
    fun = run_hyper, hyper_dmi = hyper_dmi, samples_list = samples_list
)


seq_list <- seq_len(ncore)
res <- lapply(seq_list, function(i) {
    set.seed(config_list[[i]]@seed)
    run_hyper(config_list[[i]], hyper_dmi, samples_list[[i]])
})

# hyper_dmi@data
# names(samples_list[[1]])
# slotNames(samples_list[[1]])

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
