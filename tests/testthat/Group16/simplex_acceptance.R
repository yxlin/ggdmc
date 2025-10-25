# Calculate probability that K uniform [0,1] random variables sum to <= 1
# This is the volume of a K-dimensional simplex

calculate_simplex_volume <- function(K) {
    # Volume of simplex = 1/K!
    volume <- 1 / factorial(K)
    return(volume)
}

cat("=== Simplex Acceptance Rates ===\n\n")

# 2 skills: 2^2 = 4 profiles, estimate 3
K_2skill <- 3
vol_2 <- calculate_simplex_volume(K_2skill)
cat("2 skills (", K_2skill, "profile probs to estimate):\n", sep="")
cat("  Probability sum <= 1: ", vol_2, "\n", sep="")
cat("  Acceptance rate: ", sprintf("%.2f%%", vol_2*100), "\n", sep="")
cat("  Rejection rate: ", sprintf("%.2f%%", (1-vol_2)*100), "\n\n", sep="")

# 3 skills: 2^3 = 8 profiles, estimate 7
K_3skill <- 7
vol_3 <- calculate_simplex_volume(K_3skill)
cat("3 skills (", K_3skill, "profile probs to estimate):\n", sep="")
cat("  Probability sum <= 1: ", vol_3, "\n", sep="")
cat("  Acceptance rate: ", sprintf("%.4f%%", vol_3*100), "\n", sep="")
cat("  Rejection rate: ", sprintf("%.4f%%", (1-vol_3)*100), "\n\n", sep="")

# Ratio
cat("Ratio: 3-skill is ", vol_2/vol_3, "x HARDER than 2-skill!\n\n", sep="")

# Simulation to verify
set.seed(123)
n_sim <- 100000

cat("=== Monte Carlo Verification ===\n\n")

# 2-skill case
sims_2 <- matrix(runif(n_sim * K_2skill), ncol=K_2skill)
valid_2 <- rowSums(sims_2) <= 1
cat("2 skills - Simulated acceptance rate: ", 
    sprintf("%.2f%%", mean(valid_2)*100), 
    " (expected: ", sprintf("%.2f%%", vol_2*100), ")\n", sep="")

# 3-skill case
sims_3 <- matrix(runif(n_sim * K_3skill), ncol=K_3skill)
valid_3 <- rowSums(sims_3) <= 1
cat("3 skills - Simulated acceptance rate: ", 
    sprintf("%.4f%%", mean(valid_3)*100), 
    " (expected: ", sprintf("%.4f%%", vol_3*100), ")\n\n", sep="")

# Expected number of proposals needed
cat("=== Expected Proposals Needed ===\n\n")
cat("2 skills: ~", round(1/vol_2), " proposals per acceptance\n", sep="")
cat("3 skills: ~", format(round(1/vol_3), big.mark=","), 
    " proposals per acceptance\n\n", sep="")

# MCMC perspective
cat("=== MCMC Perspective (1000 iterations) ===\n\n")
n_iter <- 1000
cat("2 skills: ~", round(n_iter * vol_2), 
    " acceptances (", sprintf("%.1f%%", vol_2*100), ")\n", sep="")
cat("3 skills: ~", round(n_iter * vol_3, 2), 
    " acceptances (", sprintf("%.3f%%", vol_3*100), ")\n", sep="")
cat("  → Chain appears FLAT (almost no movement)\n")

