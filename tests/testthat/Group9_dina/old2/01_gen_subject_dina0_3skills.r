#!/usr/bin/env Rscript
# Generate CDM Data with Multiple Sample Sizes for Asymptotic Behaviour Study
# q(save = "no")
cat("\n\n--------- Generate CDM Data (DINA + No MVN, Multiple N's) ----------\n")
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

# =============================================================================
# Model Setup
# =============================================================================
Q <- as.matrix(CDM::data.ecpe$q.matrix[8:10, ])
rownames(Q) <- c("Item1", "Item2", "Item3")

K <- ncol(Q)
means <- rep(0, ncol(Q))
means
Sigma <- cdModel::build_correlation_matrix(K, sigma)

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        mean1 = "1", mean2 = "1", mean3 = "1",
        sigma = "1",
        slip1 = "1", slip2 = "1", slip3 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

cat("Model built with", model@npar, "parameters:\n")
cat(" ", paste(model@pnames, collapse = ", "), "\n\n")

# Using discrete profile probabilities (no MVN)
use_mvn <- TRUE

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = use_mvn
)


# =============================================================================
# True Parameter Values
# =============================================================================
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)

sim_p_vector <- c(
    guess1 = .10, guess2 = .20, guess3 = .30,
    mean1 = means[1], mean2 = means[2], mean3 = means[3], sigma = sigma,
    slip1 = .12, slip2 = .15, slip3 = .20
)


true_vector <- c(
    guess1 = .10, guess2 = .20, guess3 = .30,
    pi_000 = skill_probs$probability[1],
    pi_100 = skill_probs$probability[2],
    pi_010 = skill_probs$probability[3],
    pi_110 = skill_probs$probability[4],
    pi_001 = skill_probs$probability[5],
    pi_101 = skill_probs$probability[6],
    pi_011 = skill_probs$probability[7],
    slip1 = .12, slip2 = .15, slip3 = .20
)


# =============================================================================
# Generate Data for Multiple Sample Sizes
# =============================================================================

# Sample sizes for asymptotic behavior study
# From small (100) to very large (100,000)
N_values <- c(100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)

cat("Generating datasets for N =", paste(N_values, collapse = ", "), "\n")
cat(rep("=", 80), "\n\n", sep = "")
i <- 1

# Generate data for each sample size
for (i in seq_along(N_values)) {
    N <- N_values[i]

    cat(sprintf("[%d/%d] N = %d", i, length(N_values), N))
    cat("\n", rep("-", 80), "\n", sep = "")

    # Use different seed for each N to ensure independence
    seed_i <- 12300 + i

    # Simulate data
    cat("  Simulating data (seed =", seed_i, ")...\n")

    dat <- simulate(sub_model,
        nsim = N,
        parameter_vector = sim_p_vector,
        nschool = 1,
        debug = FALSE,
        seed = seed_i
    )

    cat("  Generated", N, "observations\n")
    cat("  Response dimensions:", dim(dat$responses), "\n")
    params_dina <- test_parameter_names(Q, "DINA", FALSE)

    model <- BuildModel(
        p_map = list(
            guess1 = "1", guess2 = "1", guess3 = "1",
            pi_00 = "1", pi_10 = "1", pi_01 = "1",
            slip1 = "1", slip2 = "1", slip3 = "1"
        ),
        factors = NULL,
        constants = NULL,
        match_map = NULL,
        accumulators = Q,
        type = "cdm",
        verbose = TRUE
    )


    model <- BuildModel(
        p_map = list(
            guess1 = "1", guess2 = "1", guess3 = "1",
            pi_000 = "1", pi_100 = "1", pi_010 = "1",
            pi_110 = "1", pi_001 = "1", pi_101 = "1",
            pi_011 = "1",
            slip1 = "1", slip2 = "1", slip3 = "1"
        ),
        factors = NULL,
        constants = NULL,
        match_map = NULL,
        accumulators = Q,
        type = "cdm",
        verbose = TRUE
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


    # dat2 <- CDM::data.ecpe$data[, c("id", paste0("E", 8:10))]
    # source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
    # dat_long <- wide2long_base(dat2)
    # dat_long$s <- 1
    # dat_responses <- dat_long[, 2:5]
    # head(dat$responses)
    # head(dat_responses)

    sub_dmis <- BuildDMI(dat$responses, model,
        q_matrix = model@cdm_info$q_matrix,
        rule = "DINA",
        use_mvn = FALSE
    )
    # sub_dmis <- BuildDMI(dat_responses, model,
    #     q_matrix = model@cdm_info$q_matrix,
    #     rule = "DINA",
    #     use_mvn = FALSE
    # )


    cat("Built DMI with", length(sub_dmis), "subject(s)\n")
    cat("use_mvn:", sub_dmis[[1]]@use_mvn, "\n")

    # Save dataset
    save_filename <- sprintf("subject_dina0_N%06d.rda", N)
    save_path <- paste0(data_dir, save_filename)

    true_vector <- c(
        guess1 = .10, guess2 = .20, guess3 = .30,
        slip1 = .12, slip2 = .15, slip3 = .20
    )


    save(model, sub_dmis, sub_priors, dat, Q, N, true_vector, file = save_path)
    cat("  Saved to:", save_filename, "\n")
}

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 1, is_pblocked = FALSE,
    seed = 9032
)

p0 <- ggdmc::plot(fits0[[1]], start = 1)
p0 <- ggdmc::plot(fits0[[2]], start = 1)
p0 <- ggdmc::plot(fits0[[3]], start = 1)

p0 <- ggdmc::plot(fits0[[1]], start = 200)
p0 <- ggdmc::plot(fits0[[2]], start = 200)
p0 <- ggdmc::plot(fits0[[3]], start = 30)
