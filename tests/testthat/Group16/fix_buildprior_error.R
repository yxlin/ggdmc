library(ggdmcPrior)

# Simulate 3-skill CDM model: 5 items + 7 profile probs + 5 slips = 17 params
model_npar <- 17

# Identify which parameters are profile probabilities
pnames <- c(
    paste0("guess", 1:5),
    "pi_000", "pi_001", "pi_010", "pi_011", "pi_100", "pi_101", "pi_110",
    paste0("slip", 1:5)
)
pi_params <- grepl("^pi_", pnames)

cat("=== WRONG: Setting p0=0 for Dirichlet ===\n")
p0_wrong <- rep(0, model_npar)
names(p0_wrong) <- pnames
p1 <- ifelse(pi_params, NA, 1.0)
names(p1) <- pnames

cat("p0 for Dirichlet params:", p0_wrong[pi_params], "\n")
cat("p1 for Dirichlet params:", p1[pi_params], "\n\n")

tryCatch({
    BuildPrior(
        p0 = p0_wrong,
        p1 = p1,
        lower = rep(NA, model_npar),
        upper = rep(NA, model_npar),
        dists = ifelse(pi_params, "dirichlet", "unif"),
        log_p = rep(TRUE, model_npar)
    )
}, error = function(e) {
    cat("❌ ERROR (expected):\n")
    cat(conditionMessage(e), "\n\n")
})


cat("=== CORRECT: Setting p0 > 0 for Dirichlet ===\n")
# For Dirichlet, p0 = alpha (concentration parameter) MUST BE > 0
p0_correct <- rep(0, model_npar)
p0_correct[pi_params] <- 2.0  # Alpha = 2 for symmetric Dirichlet
names(p0_correct) <- pnames

cat("p0 for Dirichlet params:", p0_correct[pi_params], "\n")
cat("p1 for Dirichlet params:", p1[pi_params], "\n\n")

tryCatch({
    p_prior <- BuildPrior(
        p0 = p0_correct,
        p1 = p1,
        lower = rep(NA, model_npar),
        upper = rep(NA, model_npar),
        dists = ifelse(pi_params, "dirichlet", "unif"),
        log_p = rep(TRUE, model_npar)
    )
    cat("✓ SUCCESS! BuildPrior created with", length(p_prior), "parameters\n\n")
    
    # Show the Dirichlet parameters
    cat("Dirichlet parameters:\n")
    for (name in names(p_prior)[pi_params]) {
        p <- p_prior[[name]]
        cat(sprintf("  %s: alpha=%.1f, p1=%s\n", name, p$p0, 
                   ifelse(is.na(p$p1), "NA", as.character(p$p1))))
    }
}, error = function(e) {
    cat("ERROR:\n")
    cat(conditionMessage(e), "\n")
})


cat("\n=== KEY INSIGHT ===\n")
cat("For Dirichlet distribution:\n")
cat("  - p0 = alpha (concentration parameter) MUST be > 0\n")
cat("  - p1 = NA (or alpha_K for last parameter)\n")
cat("  - Common choices:\n")
cat("    * alpha = 1: Uniform over simplex\n")
cat("    * alpha = 2: Symmetric, slightly concentrated\n")
cat("    * alpha = 0.5: Favors sparse solutions\n")

