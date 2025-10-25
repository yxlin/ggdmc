#!/usr/bin/env Rscript
# Debug likelihood comparison
cat("\n\n---------------- Debug Likelihood Comparison ----------------\n")
rm(list = ls())

pkg <- c("lbaModel", "ggdmcPrior", "ggdmc")
sapply(pkg, require, character.only = TRUE)

# Build simple LBA model
model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    accumulators = c("r1", "r2"),
    type = "lba"
)

# True parameters
p_vector <- c(A = .75, B = 1.25, mean_v.false = 1.5, mean_v.true = 2.5, t0 = .15)
sub_model <- setLBA(model)
dat <- simulate(sub_model, nsim = 100, parameter_vector = p_vector, n_subject = 1)
sub_dmis <- ggdmcModel::BuildDMI(dat, model)

# Test both methods with same parameters
test_params <- p_vector  # Use true parameters

cat("\nTesting with true parameters:\n")
print(test_params)

# Method 2: Fast (via run_subject_typeconv internals)
# We need to call the C++ directly
# For now, let's just compare by running a single iteration

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
sub_priors <- set_priors(p_prior = p_prior)

# Initialize with known parameters
nmc <- 2
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames, thin = 1)
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]],
                                      seed = 123, verbose = FALSE)

# Set the first chain to our test parameters
sub_samples@theta[1, , 1] <- test_params

de_input <- ggdmc::setDEInput(
    sub_migration_prob = 0.0,  # No migration, just likelihood evaluation
    gamma_precursor = 2.38, rp = 0.001,
    is_pblocked = FALSE,
    nparameter = as.integer(model@npar),
    nchain = as.integer(model@npar * 3),
    sub_debug = FALSE
)

config_list <- ggdmc::set_configs(
    prior = sub_priors, theta_input = sub_theta_input,
    de_input = de_input,
    ncore = 1, seed = 123
)

cat("\nInitial log-likelihood from samples object:\n")
cat("Chain 1: ", sub_samples@log_likelihoods[1, 1], "\n")

cat("\n=== Running typeconv version ===\n")
fit_typeconv <- ggdmc::run_subject_typeconv(
    config_r = config_list[[1]], dmi = sub_dmis[[1]],
    samples = sub_samples
)

cat("\nFinal log-likelihood from typeconv:\n")
cat("Chain 1 sample 1: ", fit_typeconv@log_likelihoods[1, 1], "\n")
cat("Chain 1 sample 2: ", fit_typeconv@log_likelihoods[1, 2], "\n")
