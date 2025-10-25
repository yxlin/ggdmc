library(ggdmcPrior)

# Simulate the user's scenario with 3 skills (8 profiles, 7 to estimate)
# Model with 5 items + 7 profile probs
model_npar <- 5 + 7 + 5  # guess + pi + slip = 17

# Create parameters
p0 <- rep(0, model_npar)
names(p0) <- c(
    paste0("guess", 1:5),
    "pi_000", "pi_001", "pi_010", "pi_011", "pi_100", "pi_101", "pi_110",
    paste0("slip", 1:5)
)

p1 <- rep(1.0, model_npar)
names(p1) <- names(p0)

# Mark profile probs as NA (as they should be for Dirichlet)
pi_params <- grepl("^pi_", names(p0))
p1[pi_params] <- NA

cat("=== Testing with WRONG parameter name ===\n")
cat("Using 'dist' instead of 'dists'\n\n")

# THIS IS THE BUG - using 'dist' instead of 'dists'
tryCatch({
    p_prior_wrong <- BuildPrior(
        p0 = p0,
        p1 = p1,
        lower = rep(NA, model_npar),
        upper = rep(NA, model_npar),
        dist = ifelse(pi_params, "dirichlet", "unif"),  # ❌ WRONG: 'dist' not 'dists'
        log_p = rep(TRUE, model_npar)
    )
}, error = function(e) {
    cat("ERROR:\n")
    cat(conditionMessage(e), "\n\n")
})

cat("=== Testing with CORRECT parameter name ===\n")
cat("Using 'dists' (plural)\n\n")

# THIS IS CORRECT - using 'dists'
tryCatch({
    p_prior_correct <- BuildPrior(
        p0 = p0,
        p1 = p1,
        lower = rep(NA, model_npar),
        upper = rep(NA, model_npar),
        dists = ifelse(pi_params, "dirichlet", "unif"),  # ✓ CORRECT: 'dists'
        log_p = rep(TRUE, model_npar)
    )
    cat("✓ BuildPrior succeeded!\n")
    cat("Prior created with", length(p_prior_correct), "parameters\n")
}, error = function(e) {
    cat("ERROR:\n")
    cat(conditionMessage(e), "\n")
})

