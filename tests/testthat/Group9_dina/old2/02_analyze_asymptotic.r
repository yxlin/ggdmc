#!/usr/bin/env Rscript
# Analyze Asymptotic Behavior Across Different Sample Sizes

cat("\n\n-------------------- Asymptotic Behavior Analysis --------------------\n")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggplot2", "dplyr", "tidyr")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
data_dir <- paste0(home_dir, "Group9_gen_cdm/data/")
results_dir <- paste0(home_dir, "Group9_gen_cdm/results/")

# Create results directory
if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
    cat("Created directory:", results_dir, "\n")
}

# =============================================================================
# Load Metadata
# =============================================================================

metadata_file <- paste0(data_dir, "generation_metadata.rda")
if (!file.exists(metadata_file)) {
    stop("Metadata file not found. Run 01_gen_subject_dina0.r first.")
}

load(metadata_file)

cat("\nLoaded metadata for", nrow(all_metadata), "datasets\n")
cat("Sample sizes:", paste(N_values, collapse = ", "), "\n")
cat("\nTrue parameter values:\n")
print(p_vector)
cat("\n")

# =============================================================================
# Fit Models and Extract Estimates
# =============================================================================

# Storage for results
all_results <- data.frame(
    N = integer(),
    parameter = character(),
    true_value = numeric(),
    estimate = numeric(),
    se = numeric(),
    lower_95 = numeric(),
    upper_95 = numeric(),
    bias = numeric(),
    coverage = logical(),
    stringsAsFactors = FALSE
)

# MCMC settings (adjust as needed)
nmc <- 500      # iterations per chain
nchain <- 3     # number of chains
burnin <- 0.5   # proportion to discard

cat(rep("=", 80), "\n", sep = "")
cat("MCMC Settings:\n")
cat("  Iterations per chain:", nmc, "\n")
cat("  Number of chains:", nchain, "\n")
cat("  Burn-in proportion:", burnin, "\n")
cat(rep("=", 80), "\n\n", sep = "")

# Process each dataset
for (i in seq_along(N_values)) {
    N <- N_values[i]

    cat(sprintf("[%d/%d] Processing N = %d", i, length(N_values), N))
    cat("\n", rep("-", 80), "\n", sep = "")

    # Load dataset
    data_file <- paste0(data_dir, all_metadata$filename[i])
    load(data_file)

    cat("  Loaded data with N =", N, "\n")

    # Run MCMC
    cat("  Running MCMC (", nmc, " iterations × ", nchain, " chains)...\n", sep = "")

    fits <- StartSampling_subject(sub_dmis[[1]], sub_priors,
        ncore = nchain,
        nmc = nmc,
        nchain = nchain,
        sub_migration_prob = 0.05,
        thin = 1,
        is_pblocked = FALSE,
        seed = 9000 + i
    )

    cat("  MCMC complete\n")

    # Extract posterior samples (after burn-in)
    start_idx <- ceiling(nmc * burnin) + 1

    # Combine chains
    theta_samples <- NULL
    for (ch in 1:nchain) {
        theta_ch <- fits[[ch]]@theta[start_idx:nmc, , ]
        theta_samples <- rbind(theta_samples, theta_ch)
    }

    # Calculate statistics for each parameter
    param_names <- dimnames(theta_samples)[[2]]

    for (pname in param_names) {
        samples <- theta_samples[, pname]

        estimate <- mean(samples)
        se <- sd(samples)
        q <- quantile(samples, probs = c(0.025, 0.975))
        true_val <- p_vector[pname]

        bias <- estimate - true_val
        coverage <- (true_val >= q[1]) && (true_val <= q[2])

        # Store results
        all_results <- rbind(all_results, data.frame(
            N = N,
            parameter = pname,
            true_value = true_val,
            estimate = estimate,
            se = se,
            lower_95 = q[1],
            upper_95 = q[2],
            bias = bias,
            coverage = coverage,
            stringsAsFactors = FALSE
        ))

        cat(sprintf("    %6s: estimate = %6.4f (true = %6.4f), bias = %+7.4f, SE = %6.4f\n",
                   pname, estimate, true_val, bias, se))
    }

    cat("\n")
}

# =============================================================================
# Save Results
# =============================================================================

results_file <- paste0(results_dir, "asymptotic_results.rda")
save(all_results, N_values, p_vector, file = results_file)
cat("Results saved to:", results_file, "\n\n")

# =============================================================================
# Summary Statistics
# =============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("SUMMARY STATISTICS\n")
cat(rep("=", 80), "\n\n", sep = "")

# Calculate summary by parameter
summary_by_param <- all_results %>%
    group_by(parameter) %>%
    summarise(
        mean_bias = mean(bias),
        mean_abs_bias = mean(abs(bias)),
        mean_se = mean(se),
        coverage_rate = mean(coverage),
        .groups = "drop"
    )

cat("Overall statistics by parameter:\n")
print(as.data.frame(summary_by_param))
cat("\n")

# Calculate summary by N
summary_by_N <- all_results %>%
    group_by(N) %>%
    summarise(
        mean_abs_bias = mean(abs(bias)),
        mean_se = mean(se),
        coverage_rate = mean(coverage),
        .groups = "drop"
    )

cat("Overall statistics by sample size:\n")
print(as.data.frame(summary_by_N))
cat("\n")

# =============================================================================
# Create Plots
# =============================================================================

cat("Creating diagnostic plots...\n\n")

# Plot 1: Bias vs N
p1 <- ggplot(all_results, aes(x = N, y = bias, color = parameter)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_log10(breaks = N_values) +
    labs(
        title = "Bias vs Sample Size",
        subtitle = "Bias should approach 0 as N increases",
        x = "Sample Size (N)",
        y = "Bias (Estimate - True)",
        color = "Parameter"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(results_dir, "bias_vs_N.pdf"), p1, width = 10, height = 6)

# Plot 2: Standard Error vs N
p2 <- ggplot(all_results, aes(x = N, y = se, color = parameter)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_x_log10(breaks = N_values) +
    scale_y_log10() +
    labs(
        title = "Standard Error vs Sample Size",
        subtitle = "SE should decrease as O(1/√N)",
        x = "Sample Size (N)",
        y = "Standard Error",
        color = "Parameter"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(results_dir, "se_vs_N.pdf"), p2, width = 10, height = 6)

# Plot 3: Coverage Rate vs N
coverage_by_N <- all_results %>%
    group_by(N) %>%
    summarise(coverage_rate = mean(coverage), .groups = "drop")

p3 <- ggplot(coverage_by_N, aes(x = N, y = coverage_rate)) +
    geom_line(size = 1, color = "blue") +
    geom_point(size = 3, color = "blue") +
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
    scale_x_log10(breaks = N_values) +
    ylim(0, 1) +
    labs(
        title = "95% Credible Interval Coverage Rate vs Sample Size",
        subtitle = "Should be approximately 0.95 (dashed line)",
        x = "Sample Size (N)",
        y = "Coverage Rate"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(results_dir, "coverage_vs_N.pdf"), p3, width = 10, height = 6)

# Plot 4: Estimates vs N for each parameter
p4 <- ggplot(all_results, aes(x = N, y = estimate, color = parameter)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_hline(data = all_results %>%
                   select(parameter, true_value) %>%
                   distinct(),
               aes(yintercept = true_value, color = parameter),
               linetype = "dashed") +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95, fill = parameter),
                alpha = 0.2, color = NA) +
    scale_x_log10(breaks = N_values) +
    labs(
        title = "Parameter Estimates vs Sample Size",
        subtitle = "Estimates should converge to true values (dashed lines)",
        x = "Sample Size (N)",
        y = "Estimate",
        color = "Parameter",
        fill = "Parameter"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(results_dir, "estimates_vs_N.pdf"), p4, width = 10, height = 6)

cat("Plots saved to:", results_dir, "\n")
cat("  - bias_vs_N.pdf\n")
cat("  - se_vs_N.pdf\n")
cat("  - coverage_vs_N.pdf\n")
cat("  - estimates_vs_N.pdf\n\n")

# =============================================================================
# Asymptotic Normality Check
# =============================================================================

cat(rep("=", 80), "\n", sep = "")
cat("ASYMPTOTIC NORMALITY CHECK\n")
cat(rep("=", 80), "\n\n", sep = "")

# For large N, √N × (estimate - true) should be approximately N(0, σ²)
# We'll check this for the largest sample size

N_large <- max(N_values)
results_large_N <- all_results %>% filter(N == N_large)

cat("Checking asymptotic normality for N =", N_large, "\n\n")

for (pname in unique(results_large_N$parameter)) {
    result <- results_large_N %>% filter(parameter == pname)

    standardized <- sqrt(N_large) * result$bias

    cat(sprintf("  %6s: √N × bias = %+7.4f (SE = %6.4f)\n",
               pname, standardized, sqrt(N_large) * result$se))
}

cat("\n")
cat("Interpretation:\n")
cat("  - For large N, these standardized values should be O(1)\n")
cat("  - If they continue to grow with N, there may be bias issues\n")
cat("  - Standard errors should also be O(1) after √N scaling\n\n")

cat(rep("=", 80), "\n", sep = "")
cat("Analysis Complete!\n")
cat(rep("=", 80), "\n\n", sep = "")
