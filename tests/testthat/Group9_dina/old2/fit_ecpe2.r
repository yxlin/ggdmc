q(save = "no")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
Q <- as.matrix(CDM::data.ecpe$q.matrix[8:10, ])

rownames(Q) <- c("Item1", "Item2", "Item3")

K <- ncol(Q)
means <- rep(0, ncol(Q))
sigma <- 0
Sigma <- cdModel::build_correlation_matrix(K, sigma)
Sigma

skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        pi_000 = "1", pi_100 = "1", pi_010 = "1",
        pi_110 = "1", pi_001 = "1", pi_101 = "1",
        pi_011 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1"
    ),
    factors = NULL,
    constants = c(
        pi_000 = skill_probs$probability[1],
        pi_100 = skill_probs$probability[2],
        pi_010 = skill_probs$probability[3],
        pi_110 = skill_probs$probability[4],
        pi_001 = skill_probs$probability[5],
        pi_101 = skill_probs$probability[6],
        pi_011 = skill_probs$probability[7]
    ),
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
    constants = c(
        pi_000 = .089,
        pi_100 = .031,
        pi_010 = .133,
        pi_110 = .084,
        pi_001 = .121,
        pi_101 = .052,
        pi_011 = .096
    ),
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
model@pnames

# p0 <- rep(0, model@npar)
p0 <- c(
    0, 0, 0, 0.0001, 0.0001, 0.0001,
    0.0001, 0.0001, 0.0001,
    0.0001, 0, 0, 0
)

names(p0) <- model@pnames

p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = c(1, 1, 1, NA, NA, NA, NA, NA, NA, NA, 1, 1, 1), # ✅ NA for Dirichlet params
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = c(
        "unif", "unif", "unif", # guess parameters
        "dirichlet", "dirichlet", "dirichlet",
        "dirichlet", "dirichlet", "dirichlet",
        "dirichlet",
        # ✅ profile probabilities
        "unif", "unif", "unif" # slip parameters
    ),
    log_p = rep(TRUE, model@npar)
)

# p_prior <- ggdmcPrior::BuildPrior(
#     p0 = p0,
#     p1 = rep(1.0, model@npar),
#     lower = rep(NA, model@npar),
#     upper = rep(NA, model@npar),
#     dist = rep("unif", model@npar),
#     log_p = rep(TRUE, model@npar)
# )
sub_priors <- set_priors(p_prior = p_prior)

dat2 <- CDM::data.ecpe$data[, c("id", paste0("E", 8:10))]
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
dat_long <- wide2long_base(dat2)
dat_long$s <- 1
dat_responses <- dat_long[, 2:5]

sub_dmis <- BuildDMI(dat_responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = FALSE
)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 1, is_pblocked = FALSE,
    seed = 9032
)
fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 1, is_pblocked = FALSE
)

# p0 <- ggdmc::plot(fits0[[1]], start = 1)
# p0 <- ggdmc::plot(fits0[[2]], start = 1)
# p0 <- ggdmc::plot(fits0[[3]], start = 1)

# p0 <- ggdmc::plot(fits0[[1]], start = 200)
# p0 <- ggdmc::plot(fits0[[2]], start = 200)
# p0 <- ggdmc::plot(fits0[[3]], start = 200)

# p0 <- ggdmc::plot(fits1[[1]])
# p0 <- ggdmc::plot(fits1[[2]])
# p0 <- ggdmc::plot(fits1[[3]])


fits <- fits1
fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")
ncore <- 3
for (i in seq_len(ncore)) {
    hat <- gelman(fits[[i]])
    cat("Chain", i, "mpsrf = ", hat$mpsrf, "\n")
}


options(digits = 2)

est_theta <- ggdmc::summary(fit)
est_theta$quantiles
est_theta$statistics
#           5%  50% 97.5%
# guess1 0.807 0.90  1.00
# guess2 0.438 0.71  0.99
# guess3 0.350 0.64  0.98
# slip1  0.011 0.10  0.20
# slip2  0.033 0.31  0.58
# slip3  0.036 0.32  0.66
# > est_theta$statistics
#        Mean    SD Naive SE Time-series SE
# guess1 0.90 0.059   0.0011         0.0057
# guess2 0.71 0.175   0.0032         0.0364
# guess3 0.65 0.195   0.0036         0.0268
# slip1  0.10 0.059   0.0011         0.0058
# slip2  0.31 0.176   0.0032         0.0379
# slip3  0.33 0.195   0.0036         0.0268

# Matching profile probabilities
# ggdmc
#            5%   50% 97.5%
# guess1 0.6372 0.682 0.768
# guess2 0.4465 0.509 0.568
# guess3 0.0033 0.035 0.162
# slip1  0.0010 0.012 0.048
# slip2  0.1016 0.149 0.196
# slip3  0.0024 0.023 0.091
# CDM::din
#            5%   50% 97.5%
# guess1  0.499 0.641 0.782
# guess2  0.140 0.337 0.535
# guess3  0.065 0.269 0.473
# slip1  -0.034 0.013 0.061
# slip2  -0.038 0.079 0.197
# slip3  -0.056 0.081 0.217

# dat <- CDM::data.ecpe$data[, -1]
# Q <- CDM::data.ecpe$q.matrix
ecpe <- CDM::din(data = dat2[, -1], q.matrix = Q)
param <- CDM::IRT.se(ecpe, extended = TRUE)
p <- split(param, param$partype)
cdm_result <- tibble::tibble(
    guess = p$guess$est,
    slip = p$slip$est,
    guess005 = p$guess$est - 1.96 * p$guess$se,
    guess975 = p$guess$est + 1.96 * p$guess$se,
    slip005 = p$slip$est - 1.96 * p$slip$se,
    slip975 = p$slip$est + 1.96 * p$slip$se
)

# Rearrange to match ggdmc format
nitem <- length(p$guess$est)
cdm_formatted <- data.frame(
    `5%` = c(
        p$guess$est - 1.96 * p$guess$se, # guess lower bounds
        p$slip$est - 1.96 * p$slip$se # slip lower bounds
    ),
    `50%` = c(p$guess$est, p$slip$est), # medians (estimates)
    `97.5%` = c(
        p$guess$est + 1.96 * p$guess$se, # guess upper bounds
        p$slip$est + 1.96 * p$slip$se # slip upper bounds
    ),
    check.names = FALSE
)
rownames(cdm_formatted) <- c(
    paste0("guess", 1:nitem),
    paste0("slip", 1:nitem)
)

cat("\n=== CDM Package Results ===\n")
print(round(cdm_formatted, 3))

cat("\n=== Original cdm_result ===\n")
print(cdm_result)


p$probs$est
p$probs$se
ecpe$attribute.patt
