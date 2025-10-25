#!/usr/bin/env Rscript
cat("\n\n-------------------- Generate CDM Data (DINO + MVN, Multiple N's) --------------------\n")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
data_dir <- paste0(home_dir, "Group9_gen_cdm/data/")

# Ensure data directory exists
if (!dir.exists(data_dir)) {
    dir.create(data_dir, recursive = TRUE)
    cat("Created directory:", data_dir, "\n")
}

use_mvn <- TRUE
Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1,
    1, 0,
    0, 1
), ncol = 2, byrow = TRUE)

colnames(Q) <- c("Algebra", "Geometry")
rownames(Q) <- c("Item1", "Item2", "Item3", "Item4", "Item5")

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1", guess4 = "1", guess5 = "1",
        mean1 = "1", mean2 = "1", sigma = "1",
        slip1 = "1", slip2 = "1", slip3 = "1", slip4 = "1", slip5 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability,
    rule = "DINO",
    use_mvn = use_mvn
)

# True parameter values
means <- c(0.5, 0.2)
sigma <- 0.2

p_vector <- c(
    guess1 = .1, guess2 = .2, guess3 = .3, guess4 = .4, guess5 = .5,
    mean1 = means[1], mean2 = means[2], sigma = sigma,
    slip1 = .2, slip2 = .4, slip3 = .6, slip4 = .8, slip5 = .9
)

# Prior setup (same for all cases)
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(1.0, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
sub_priors <- set_priors(p_prior = p_prior)


# Generate and save data for multiple sample sizes
N_values <- c(50, 500, 5000, 50000, 500000)

cat("\n=======================================================================\n")
cat("Generating CDM data for", length(N_values), "different sample sizes\n")
cat("=======================================================================\n\n")

for (i in seq_along(N_values)) {
    N <- N_values[i]

    cat("\n-----------------------------------------------------------------------\n")
    cat("Case ", i, "/", length(N_values), ": N = ", N, "\n")
    cat("-----------------------------------------------------------------------\n")

    # Simulate data
    set.seed(123) # Same seed for reproducibility
    dat <- simulate(sub_model,
        nsim = N,
        parameter_vector = p_vector,
        nschool = 1,
        debug = FALSE
    )

    cat("Simulated", N, "observations\n")
    cat("Response dimensions:", dim(dat$responses), "\n")

    # Build DMI
    sub_dmis <- BuildDMI(dat$responses, model,
        q_matrix = model@cdm_info$q_matrix,
        profile_probability = model@cdm_info$profile_probability,
        rule = "DINO",
        use_mvn = use_mvn
    )

    cat("Built DMI with", length(sub_dmis), "subject(s)\n")
    cat("use_mvn:", sub_dmis[[1]]@use_mvn, "\n")

    # Save to file
    save_path <- paste0(data_dir, "subject_dina0_mvn_N", N, ".rda")
    save(model, dat, p_vector, sub_dmis, Q, sub_priors, N, file = save_path)
    cat("Saved to:", save_path, "\n")

    # Verify file size
    file_info <- file.info(save_path)
    cat("File size:", round(file_info$size / 1024, 2), "KB\n")
}

cat("\n=======================================================================\n")
cat("Data generation complete!\n")
cat("=======================================================================\n\n")

# Print summary
cat("Generated files:\n")
for (N in N_values) {
    save_path <- paste0(data_dir, "subject_dina0_mvn_N", N, ".rda")
    if (file.exists(save_path)) {
        file_info <- file.info(save_path)
        cat(sprintf(
            "  - N=%6d: %s (%.2f KB)\n", N, basename(save_path),
            file_info$size / 1024
        ))
    }
}
cat("\n")
