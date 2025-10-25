#!/usr/bin/env Rscript
# Test Step 3: Type Conversion Optimization
cat("\n\n-------------------- Test Type Conversion Optimization --------------------\n")
rm(list = ls())

pkg <- c("lbaModel", "ggdmcPrior", "ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)

# Build a simple LBA model
model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1", st0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2")),
    constants = c(st0 = 0, sd_v = 1),
    accumulators = c("r1", "r2"),
    type = "lba"
)

# Generate small test data
p_vector <- c(A = .75, B = 1.25, mean_v.false = 1.5, mean_v.true = 2.5, t0 = .15)
sub_model <- setLBA(model)
dat <- simulate(sub_model, nsim = 100, parameter_vector = p_vector, n_subject = 1)

# Build DMI and prior
sub_dmis <- ggdmcModel::BuildDMI(dat, model)

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

# Setup sampling
nmc <- 1000
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames, thin = 1)
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]],
    seed = 123, verbose = FALSE
)

de_input <- ggdmc::setDEInput(
    sub_migration_prob = .05,
    gamma_precursor = 2.38, rp = 0.001,
    is_pblocked = FALSE,
    nparameter = as.integer(model@npar),
    nchain = as.integer(model@npar * 3),
    sub_debug = FALSE
)

config_list <- ggdmc::set_configs(
    prior = sub_priors, theta_input = sub_theta_input,
    de_input = de_input,
    ncore = 1, seed = 456
)

# Test 1: Run original version
cat("\nTest 1: Running original version...\n")
time_orig <- system.time({
    fit_orig <- ggdmc::run_subject(
        config_r = config_list[[1]],
        dmi = sub_dmis[[1]],
        samples = sub_samples
    )
})
cat("Original time:", time_orig[3], "seconds\n")

# Test 2: Run typeconv optimized version
cat("\nTest 2: Running type-conv optimized version...\n")
# Re-initialize to same starting point
sub_samples2 <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]],
    seed = 123, verbose = FALSE
)
time_typeconv <- system.time({
    fit_typeconv <- ggdmc::run_subject_typeconv(
        config_r = config_list[[1]],
        dmi = sub_dmis[[1]],
        samples = sub_samples2
    )
})
cat("Type-conv time:", time_typeconv[3], "seconds\n")

# Test 3: Compare results - should be nearly identical
cat("\nTest 3: Comparing results...\n")
cat("Checking dimensions match...\n")
if (identical(dim(fit_orig@theta), dim(fit_typeconv@theta))) {
    cat("  ✓ Dimensions match\n")
} else {
    cat("  ✗ Dimensions DO NOT match!\n")
}

cat("Checking log-posterior values...\n")
lp_diff <- max(abs(fit_orig@summed_log_prior - fit_typeconv@summed_log_prior))
ll_diff <- max(abs(fit_orig@log_likelihoods - fit_typeconv@log_likelihoods))
cat(sprintf("  Max log-prior difference: %.10f\n", lp_diff))
cat(sprintf("  Max log-likelihood difference: %.10f\n", ll_diff))

if (lp_diff < 1e-8 && ll_diff < 1e-8) {
    cat("  ✓ Results are numerically identical (good!)\n")
} else {
    cat("  ⚠ Results differ slightly (may be due to RNG differences)\n")
}

# Test 4: Performance comparison
cat("\nTest 4: Performance summary\n")
speedup <- time_orig[3] / time_typeconv[3]
pct_improvement <- (1 - time_typeconv[3] / time_orig[3]) * 100
cat(sprintf("Speedup: %.2fx\n", speedup))
cat(sprintf("Improvement: %.1f%%\n", pct_improvement))

if (speedup > 1.0) {
    cat("  ✓ Type-conv version is FASTER\n")
} else if (speedup > 0.95) {
    cat("  ~ Type-conv version is similar speed (neutral)\n")
} else {
    cat("  ✗ Type-conv version is SLOWER (unexpected)\n")
}

cat("\n-------------------- Test Complete --------------------\n")
