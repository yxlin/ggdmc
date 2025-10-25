#!/usr/bin/env Rscript
# Analysis and visualization of asymptotic behavior test results

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║      Asymptotic Behavior Analysis & Visualization             ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

# Setup ====================================================================
rm(list = ls())
pkg <- c("ggdmc", "ggplot2")
pkg_ok <- sapply(pkg, require, character.only = TRUE)

if (!all(pkg_ok)) {
  stop("Failed to load required packages: ", names(pkg_ok)[!pkg_ok])
}

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
results_dir <- file.path(home_dir, "Group10_cdm_subjects/results")
plots_dir <- file.path(results_dir, "plots")

# Create plots directory if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Load Results =============================================================
summary_path <- file.path(results_dir, "systematic_test_summary.rda")

if (!file.exists(summary_path)) {
  stop("Summary file not found. Run 01_fit_subject_dina0_systematic.r first.")
}

cat("Loading results from:", summary_path, "\n\n")
load(summary_path)

cat("Available sample sizes:", paste(names(all_results), collapse = ", "), "\n\n")

# 1. Convergence Analysis ==================================================
cat("1. CONVERGENCE ANALYSIS\n")
cat("══════════════════════════════════════════════════════════════\n\n")

print(summary_table[, c("N", "Overall_MPSRF", "Converged")])

# Plot MPSRF vs N
p_convergence <- ggplot(summary_table, aes(x = N, y = Overall_MPSRF)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(aes(color = Converged), size = 3) +
  geom_hline(yintercept = 1.1, linetype = "dashed", color = "red") +
  scale_color_manual(values = c("TRUE" = "green", "FALSE" = "red"),
                     labels = c("TRUE" = "Converged", "FALSE" = "Not Converged")) +
  scale_x_log10(labels = scales::comma) +
  labs(title = "Convergence Diagnostics (Gelman-Rubin MPSRF)",
       x = "Sample Size (N, log scale)",
       y = "MPSRF",
       color = "Status") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_convergence)
ggsave(file.path(plots_dir, "convergence_mpsrf.pdf"), p_convergence, width = 8, height = 6)
cat("✓ Saved: convergence_mpsrf.pdf\n\n")

# 2. Estimation Accuracy ===================================================
cat("2. ESTIMATION ACCURACY (ASYMPTOTIC BEHAVIOR)\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# Extract bias and MAE for each parameter across sample sizes
param_names <- colnames(all_results[[1]]$estimates$bias)
n_params <- length(param_names)

# Create data frame for plotting
bias_df <- data.frame()
mae_df <- data.frame()

for (n_str in names(all_results)) {
  res <- all_results[[n_str]]
  n <- res$n

  for (param in param_names) {
    bias_df <- rbind(bias_df, data.frame(
      N = n,
      Parameter = param,
      Bias = res$estimates$bias[param]
    ))

    mae_df <- rbind(mae_df, data.frame(
      N = n,
      Parameter = param,
      MAE = res$estimates$mae[param]
    ))
  }
}

# Plot bias vs N for each parameter
p_bias <- ggplot(bias_df, aes(x = N, y = Bias, color = Parameter)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "Parameter Bias vs Sample Size",
       subtitle = "Bias = Median Estimate - True Value",
       x = "Sample Size (N, log scale)",
       y = "Bias") +
  theme_minimal() +
  theme(legend.position = "right")

print(p_bias)
ggsave(file.path(plots_dir, "bias_by_parameter.pdf"), p_bias, width = 10, height = 6)
cat("✓ Saved: bias_by_parameter.pdf\n\n")

# Plot MAE vs N for each parameter
p_mae <- ggplot(mae_df, aes(x = N, y = MAE, color = Parameter)) +
  geom_line(size = 0.8) +
  geom_point(size = 2) +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(title = "Median Absolute Error vs Sample Size",
       subtitle = "MAE = |Median Estimate - True Value|",
       x = "Sample Size (N, log scale)",
       y = "MAE (log scale)") +
  theme_minimal() +
  theme(legend.position = "right")

print(p_mae)
ggsave(file.path(plots_dir, "mae_by_parameter.pdf"), p_mae, width = 10, height = 6)
cat("✓ Saved: mae_by_parameter.pdf\n\n")

# Plot mean MAE
p_mean_mae <- ggplot(summary_table, aes(x = N, y = Mean_MAE)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(size = 3, color = "steelblue") +
  scale_x_log10(labels = scales::comma) +
  scale_y_log10() +
  labs(title = "Mean Absolute Error vs Sample Size",
       subtitle = "Average across all parameters",
       x = "Sample Size (N, log scale)",
       y = "Mean MAE (log scale)") +
  theme_minimal()

print(p_mean_mae)
ggsave(file.path(plots_dir, "mean_mae.pdf"), p_mean_mae, width = 8, height = 6)
cat("✓ Saved: mean_mae.pdf\n\n")

# 3. Computational Time Analysis ===========================================
cat("3. COMPUTATIONAL TIME ANALYSIS\n")
cat("══════════════════════════════════════════════════════════════\n\n")

print(summary_table[, c("N", "Stage1_Time", "Stage2_Time", "Total_Time")])

# Reshape for plotting
time_df <- data.frame(
  N = rep(summary_table$N, 2),
  Stage = rep(c("Stage 1", "Stage 2"), each = nrow(summary_table)),
  Time = c(summary_table$Stage1_Time, summary_table$Stage2_Time)
)

p_time <- ggplot(time_df, aes(x = N, y = Time, fill = Stage)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_x_log10(labels = scales::comma) +
  scale_fill_manual(values = c("Stage 1" = "steelblue", "Stage 2" = "coral")) +
  labs(title = "Computational Time vs Sample Size",
       x = "Sample Size (N, log scale)",
       y = "Time (seconds)",
       fill = "Sampling Stage") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_time)
ggsave(file.path(plots_dir, "computational_time.pdf"), p_time, width = 8, height = 6)
cat("✓ Saved: computational_time.pdf\n\n")

# 4. Detailed Parameter Comparison =========================================
cat("4. DETAILED PARAMETER COMPARISON\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# Print comparison tables for smallest and largest N
if (length(all_results) >= 2) {
  n_values <- sort(as.numeric(names(all_results)))
  smallest_n <- as.character(n_values[1])
  largest_n <- as.character(n_values[length(n_values)])

  cat(sprintf("Smallest sample size (N = %s):\n", smallest_n))
  print(all_results[[smallest_n]]$estimates$comparison_table)

  cat(sprintf("\n\nLargest sample size (N = %s):\n", largest_n))
  print(all_results[[largest_n]]$estimates$comparison_table)
}

# 5. Asymptotic Rate Analysis ==============================================
cat("\n\n5. ASYMPTOTIC RATE ANALYSIS\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# Fit power law: MAE ~ N^(-alpha)
# log(MAE) = log(c) - alpha * log(N)

if (nrow(summary_table) >= 3) {
  log_n <- log10(summary_table$N)
  log_mae <- log10(summary_table$Mean_MAE)

  fit_lm <- lm(log_mae ~ log_n)
  alpha <- -coef(fit_lm)[2]
  intercept <- coef(fit_lm)[1]

  cat("Power law fit: MAE ~ N^(-alpha)\n")
  cat(sprintf("  Estimated alpha: %.4f\n", alpha))
  cat(sprintf("  R-squared: %.4f\n", summary(fit_lm)$r.squared))
  cat("\nInterpretation:\n")
  cat(sprintf("  MAE decreases at rate ~ 1/N^%.2f\n", alpha))

  if (abs(alpha - 0.5) < 0.1) {
    cat("  → Consistent with parametric rate (1/√N) ✓\n")
  } else if (alpha > 0.3) {
    cat("  → Reasonably good convergence rate\n")
  } else {
    cat("  → Convergence slower than expected\n")
  }

  # Add fitted line to mean MAE plot
  summary_table$Fitted_MAE <- 10^(intercept) * summary_table$N^(-alpha)

  p_mean_mae_fitted <- ggplot(summary_table, aes(x = N)) +
    geom_line(aes(y = Mean_MAE, color = "Observed"), size = 1) +
    geom_point(aes(y = Mean_MAE, color = "Observed"), size = 3) +
    geom_line(aes(y = Fitted_MAE, color = "Fitted"), size = 1, linetype = "dashed") +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    scale_color_manual(values = c("Observed" = "steelblue", "Fitted" = "red")) +
    labs(title = sprintf("Mean MAE vs Sample Size (α = %.3f)", alpha),
         subtitle = "Log-log plot with power law fit",
         x = "Sample Size (N, log scale)",
         y = "Mean MAE (log scale)",
         color = "") +
    theme_minimal() +
    theme(legend.position = "bottom")

  print(p_mean_mae_fitted)
  ggsave(file.path(plots_dir, "mean_mae_fitted.pdf"), p_mean_mae_fitted,
         width = 8, height = 6)
  cat("\n✓ Saved: mean_mae_fitted.pdf\n")
}

# Summary ==================================================================
cat("\n\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║                  ANALYSIS SUMMARY                              ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n\n")

cat("All plots saved to:", plots_dir, "\n")
cat("\nGenerated plots:\n")
cat("  1. convergence_mpsrf.pdf - Gelman-Rubin diagnostics\n")
cat("  2. bias_by_parameter.pdf - Parameter bias vs N\n")
cat("  3. mae_by_parameter.pdf - MAE by parameter vs N\n")
cat("  4. mean_mae.pdf - Average MAE vs N\n")
cat("  5. mean_mae_fitted.pdf - MAE with power law fit\n")
cat("  6. computational_time.pdf - Timing analysis\n\n")

cat("✓ Analysis completed successfully\n\n")
