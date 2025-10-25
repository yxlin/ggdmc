#!/usr/bin/env Rscript
# Test to identify exact MVN generation bug

cat("\n\n=========== MVN Generation Bug Identification ===========\n\n")

# The C++ code does:
# Z = N x K matrix of N(0,1)
# L = chol(Sigma, "lower")  # K x K lower triangular
# X = Z * L.t()             # N x K
# X.each_row() += mean.t()  # Add mean to each row

# Let's verify what this produces

library(MASS)

N <- 10000
K <- 2
true_sigma <- 0.2
true_means <- c(0.5, 0.2)

# Build correlation matrix
Sigma <- matrix(true_sigma, nrow=K, ncol=K)
diag(Sigma) <- 1

cat("Target correlation matrix Sigma:\n")
print(Sigma)
cat("\n")

# Method 1: CORRECT R implementation
cat("Method 1: CORRECT (using mvrnorm):\n")
set.seed(123)
X1 <- mvrnorm(n=N, mu=true_means, Sigma=Sigma)
cat("  Empirical mean:", colMeans(X1), "\n")
cat("  Empirical cor:\n")
print(cor(X1))
cat("  Off-diagonal:", cor(X1)[1,2], "\n\n")

# Method 2: Manual with t(chol())
cat("Method 2: Manual with L = t(chol(Sigma)):\n")
set.seed(123)
Z2 <- matrix(rnorm(N*K), nrow=N, ncol=K)
L2 <- t(chol(Sigma))  # Lower triangular
X2 <- Z2 %*% t(L2)  # Equivalent to Z * L^T in C++
X2 <- sweep(X2, 2, true_means, "+")
cat("  L (lower triangular):\n")
print(L2)
cat("  Empirical mean:", colMeans(X2), "\n")
cat("  Empirical cor:\n")
print(cor(X2))
cat("  Off-diagonal:", cor(X2)[1,2], "\n\n")

# Method 3: What if C++ chol("lower") is DIFFERENT?
# In Armadillo, chol(X, "lower") might return something different
# Let's test what happens if we use the UPPER chol but call it lower
cat("Method 3: WRONG - Using UPPER as if it were lower:\n")
set.seed(123)
Z3 <- matrix(rnorm(N*K), nrow=N, ncol=K)
U <- chol(Sigma)  # Upper triangular (R default)
# Pretend this is lower and use it as L
X3 <- Z3 %*% t(U)  # Using transpose of upper as if it were transpose of lower
X3 <- sweep(X3, 2, true_means, "+")
cat("  U (UPPER triangular from R's chol):\n")
print(U)
cat("  Empirical mean:", colMeans(X3), "\n")
cat("  Empirical cor:\n")
print(cor(X3))
cat("  Off-diagonal:", cor(X3)[1,2], "\n\n")

# Method 4: What if the matrix multiplication is backwards?
cat("Method 4: WRONG - Backwards multiplication:\n")
set.seed(123)
Z4 <- matrix(rnorm(N*K), nrow=N, ncol=K)
L4 <- t(chol(Sigma))
X4 <- t(L4) %*% t(Z4)  # Backwards: L^T * Z^T instead of Z * L^T
X4 <- t(X4)
X4 <- sweep(X4, 2, true_means, "+")
cat("  Empirical mean:", colMeans(X4), "\n")
cat("  Empirical cor:\n")
print(cor(X4))
cat("  Off-diagonal:", cor(X4)[1,2], "\n\n")

# Method 5: What if using L instead of L^T?
cat("Method 5: WRONG - Using L instead of L.t():\n")
set.seed(123)
Z5 <- matrix(rnorm(N*K), nrow=N, ncol=K)
L5 <- t(chol(Sigma))
X5 <- Z5 %*% L5  # Using L instead of L^T
X5 <- sweep(X5, 2, true_means, "+")
cat("  Empirical mean:", colMeans(X5), "\n")
cat("  Empirical cor:\n")
print(cor(X5))
cat("  Off-diagonal:", cor(X5)[1,2], "\n\n")

cat(rep("=", 70), "\n", sep = "")
cat("COMPARING TO PACKAGE OUTPUT\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("From test 17, package gave binary alpha correlation: 0.123 (N=10000)\n")
cat("From test 18, package data is best fit by sigma = 0.15\n\n")

cat("Which method matches?\n")
methods <- c("Correct", "Manual-correct", "Upper-as-lower", "Backwards", "No-transpose")
offdiag <- c(cor(X1)[1,2], cor(X2)[1,2], cor(X3)[1,2], cor(X4)[1,2], cor(X5)[1,2])

for (i in 1:5) {
  # Threshold to binary
  alpha_i <- (get(paste0("X", i)) > 0) * 1
  binary_cor <- cor(alpha_i)[1,2]
  cat(sprintf("%d. %s: X_cor=%.3f, alpha_cor=%.3f\n",
             i, methods[i], offdiag[i], binary_cor))
}

cat("\n")
