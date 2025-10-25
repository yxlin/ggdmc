# q(save = "no")
cat("\n\n-------------------- Generate LCDM model 0 --------------------")
rm(list = ls())
pkg <- c("ggdmc", "cdModel", "ggdmcModel")
pkg_ok <- sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc_ecosystem/ggdmc/tests/testthat/"
data_path <- paste0(home_dir, "Group9_gen_cdm/data/lcdm0.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)

Q <- matrix(c(
    1, 0, 0,
    0, 1, 1,
    1, 1, 0
), ncol = 3, byrow = TRUE)
colnames(Q) <- c("A1", "A2", "A3")
rownames(Q) <- c("Item1", "Item2", "Item3")

# Calculate skill probabilities (no correlation, sigma = 0)
K <- ncol(Q)
sigma <- 0
means <- rep(0, K)
Sigma <- cdModel::build_correlation_matrix(K, sigma)
skill_probs <- cdModel::calculate_skill_probabilities(means, Sigma)
skill_probs

model <- BuildModel(
    p_map = list(
        beta1__1 = "1", beta1__A1 = "1",
        beta2__1 = "1", beta2__A2 = "1", beta2__A3 = "1",
        beta2__A2xA3 = "1",
        beta3__1 = "1",
        beta3__A1 = "1",
        beta3__A2 = "1",
        beta3__A1xA2 = "1",
        mean1 = "1", mean2 = "1", mean3 = "1", sigma = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL, # Must enter NULL
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)

# Set DINA rule without multivariate normal
use_mvn <- TRUE
sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "LCDM",
    use_mvn = use_mvn
)
model@pnames
# -------------------- Data Simulation --------------------
# True parameter values
sim_p_vector <- c(
    beta1__1 = .1, beta1__A1 = .2,
    beta2__1 = .3, beta2__A2 = .4,
    beta2__A2xA3 = .6, beta2__A3 = .5,
    beta3__1 = .1,
    beta3__A1 = .2,
    beta3__A1xA2 = .4,
    beta3__A2 = .3,
    mean1 = means[1], mean2 = means[2], mean3 = means[3],
    sigma = sigma
)


# sim_p_vector <- c(
#     guess1 = .10, guess2 = .20, guess3 = .30,
#     pi_00 = skill_probs$probability[1],
#     pi_10 = skill_probs$probability[2],
#     pi_01 = skill_probs$probability[3],
#     slip1 = .01, slip2 = .03, slip3 = .05
# )

N <- 2000
dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = sim_p_vector,
    nschool = 1,
    debug = FALSE,
    seed = 123
)


# StartSampling_subject
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.02, thin = 1,
    is_pblocked = TRUE,
    ncore = 4,
    seed = 9032
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 1,
    seed = 9032
)

fits <- fits1


fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector, sort_name = FALSE)

#               beta1__1 beta1__A1 beta2__1 beta2__A2 beta2__A3 beta2__A2xA3
# True             0.100     0.200    0.300     0.400     0.500        0.600
#    5 Estimate    0.017     0.022    0.069     0.094     0.053        0.056
#   50 Estimate    0.160     0.335    0.390     0.611     0.465        0.578
# 97.5 Estimate    0.363     0.940    0.748     1.024     1.028        1.029
# Median-True      0.060     0.135    0.090     0.211    -0.035       -0.022
#               beta3__1 beta3__A1 beta3__A2 beta3__A1xA2   pi1     pi2    pi3
# True             0.100     0.200     0.300         0.40 0.050  0.1000 0.0600
#    5 Estimate    0.012     0.027     0.036         0.05 0.016  0.0057 0.0086
#   50 Estimate    0.139     0.335     0.403         0.52 0.194  0.0834 0.0981
# 97.5 Estimate    0.403     0.969     0.995         1.02 0.540  0.3509 0.4094
# Median-True      0.039     0.135     0.103         0.12 0.144 -0.0166 0.0381
#                   pi4   pi5     pi6   pi7
# True           0.2000 0.070  0.3000 0.080
#    5 Estimate  0.0041 0.011  0.0073 0.007
#   50 Estimate  0.0597 0.141  0.0820 0.082
# 97.5 Estimate  0.3005 0.502  0.3392 0.373
# Median-True   -0.1403 0.071 -0.2180 0.002




# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# long <- tibble::as_tibble(dat$responses)
# wide <- long2wide_base(long)
