#!/usr/bin/env Rscript
# Debug parameter count issue

cat("\n\n========== Parameter Count Debug ==========\n\n")
rm(list = ls())

library(ggdmc)
library(ggdmcPrior)
library(ggdmcModel)
library(cdModel)

Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1
), ncol = 2, byrow = TRUE)

colnames(Q) <- c("Algebra", "Geometry")
rownames(Q) <- c("Item1", "Item2", "Item3")

cat("Q matrix (J=3, K=2):\n")
print(Q)
cat("\n")

# Build model with discrete profile probabilities
model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        pi1 = "1", pi2 = "1", pi3 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

cat("\nR Model Info:\n")
cat("  model@npar =", model@npar, "\n")
cat("  Parameter names:", paste(model@pnames, collapse=", "), "\n")
cat("  Expected: 3 guess + 3 pi + 3 slip = 9 parameters\n\n")

# Try setting CDM with use_mvn = FALSE
cat("Calling setCDM with use_mvn = FALSE...\n")
cat("(Watch for debug output from C++)\n\n")

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINA",
    use_mvn = FALSE
)

cat("\nAfter setCDM:\n")
cat("  sub_model@type =", sub_model@type, "\n")

# Try to create a test parameter vector
p_vector <- c(
    guess1 = .1, guess2 = .2, guess3 = .3,
    pi1 = .1577, pi2 = .2629, pi3 = .1507,
    slip1 = .2, slip2 = .4, slip3 = .6
)

cat("\nTest parameter vector:\n")
print(p_vector)
cat("  Length:", length(p_vector), "\n\n")

# Simulate minimal data to create DMI
cat("Simulating data...\n")
dat <- simulate(sub_model,
    nsim = 10,
    parameter_vector = p_vector,
    nschool = 1,
    debug = FALSE,
    seed = 123
)

cat("Creating DMI...\n")
sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINA",
    use_mvn = FALSE
)

cat("\nDMI created successfully!\n")
cat("  sub_dmis[[1]] exists\n\n")

# Try to compute likelihood
cat("Testing likelihood calculation...\n")
cat("(This is where the [expected, input] error should appear)\n\n")

tryCatch({
    lik <- ggdmcLikelihood::compute_subject_likelihood(
        sub_dmis[[1]],
        p_vector,
        debug = FALSE
    )
    cat("\n✓ SUCCESS! Likelihood calculated without error.\n")
    cat("  Sum log-likelihood:", sum(log(lik[[1]])), "\n")
}, error = function(e) {
    cat("\n✗ ERROR:\n")
    cat("  ", e$message, "\n\n")
    cat("This confirms the parameter count mismatch.\n")
})

cat("\n========== Debug Complete ==========\n\n")
cat("Expected behavior:\n")
cat("  - C++ should print 'm_expected_nparameter 9' during setCDM or BuildDMI\n")
cat("  - Likelihood calculation should succeed with 9-parameter vector\n\n")

cat("If you see [expected, input] = 6, 9:\n")
cat("  - The C++ code still thinks expected=6\n")
cat("  - Possible causes:\n")
cat("    1. init_rule_parameters() called before set_use_mvn(FALSE)\n")
cat("    2. set_use_mvn() doesn't trigger recalculation\n")
cat("    3. Wrong code path (m_use_mvn is TRUE when it should be FALSE)\n\n")
