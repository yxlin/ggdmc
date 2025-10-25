# Visualize why 2-skill works but 3-skill fails with uniform priors

library(ggplot2)

# Generate random samples
set.seed(123)
n <- 10000

# 2-skill case: 3 parameters
cat("=== Sampling 2-skill case (3 profile probs) ===\n")
pi1 <- runif(n)
pi2 <- runif(n)
pi3 <- runif(n)
sum_2skill <- pi1 + pi2 + pi3
valid_2skill <- sum_2skill <= 1

cat("Valid samples:", sum(valid_2skill), "/", n, 
    "(", sprintf("%.1f%%", mean(valid_2skill)*100), ")\n")
cat("Invalid samples:", sum(!valid_2skill), "/", n,
    "(", sprintf("%.1f%%", mean(!valid_2skill)*100), ")\n\n")

# 3-skill case: 7 parameters
cat("=== Sampling 3-skill case (7 profile probs) ===\n")
sims_3 <- matrix(runif(n * 7), ncol=7)
sum_3skill <- rowSums(sims_3)
valid_3skill <- sum_3skill <= 1

cat("Valid samples:", sum(valid_3skill), "/", n,
    "(", sprintf("%.2f%%", mean(valid_3skill)*100), ")\n")
cat("Invalid samples:", sum(!valid_3skill), "/", n,
    "(", sprintf("%.2f%%", mean(!valid_3skill)*100), ")\n\n")

# Create histogram comparison
df <- data.frame(
    sum = c(sum_2skill, sum_3skill),
    case = rep(c("2 skills (3 params)", "3 skills (7 params)"), each=n),
    valid = c(valid_2skill, valid_3skill)
)

cat("=== Summary Statistics ===\n")
cat("2 skills - Sum range: [", round(min(sum_2skill), 2), ", ", 
    round(max(sum_2skill), 2), "]\n", sep="")
cat("2 skills - Mean sum: ", round(mean(sum_2skill), 2), "\n", sep="")
cat("3 skills - Sum range: [", round(min(sum_3skill), 2), ", ",
    round(max(sum_3skill), 2), "]\n", sep="")
cat("3 skills - Mean sum: ", round(mean(sum_3skill), 2), "\n\n", sep="")

# Save plot
png("/tmp/simplex_comparison.png", width=800, height=500)
p <- ggplot(df, aes(x=sum, fill=valid)) +
    geom_histogram(bins=50, alpha=0.7, position="identity") +
    geom_vline(xintercept=1, color="red", linewidth=1.5, linetype="dashed") +
    facet_wrap(~case, scales="free_y", ncol=1) +
    scale_fill_manual(values=c("TRUE"="green", "FALSE"="red"),
                      labels=c("TRUE"="Valid (sum ≤ 1)", "FALSE"="Invalid (sum > 1)"),
                      name="") +
    labs(x="Sum of Profile Probabilities", y="Count",
         title="Why 2-Skill Works but 3-Skill Doesn't with Uniform Priors",
         subtitle="Red line: constraint boundary (sum must be ≤ 1)") +
    theme_minimal(base_size=14) +
    theme(legend.position="top")
print(p)
dev.off()

cat("Plot saved to /tmp/simplex_comparison.png\n")

# Show MCMC trace simulation
cat("\n=== Simulated MCMC Traces (500 iterations) ===\n\n")

simulate_mcmc <- function(n_iter, n_params) {
    accepted <- numeric(n_iter)
    current_sum <- 0.5  # Start at a reasonable value
    n_accept <- 0
    
    for(i in 1:n_iter) {
        # Propose new values
        proposal <- runif(n_params)
        prop_sum <- sum(proposal)
        
        if(prop_sum <= 1) {
            # Valid proposal - accept it
            current_sum <- prop_sum
            n_accept <- n_accept + 1
        }
        accepted[i] <- current_sum
    }
    
    list(trace=accepted, acceptance_rate=n_accept/n_iter)
}

result_2 <- simulate_mcmc(500, 3)
result_3 <- simulate_mcmc(500, 7)

cat("2-skill acceptance rate:", sprintf("%.1f%%", result_2$acceptance_rate*100), "\n")
cat("3-skill acceptance rate:", sprintf("%.2f%%", result_3$acceptance_rate*100), "\n")

