# q(save = "no")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))

dat <- CDM::data.ecpe$data[, paste0("E", 8:10)]
Q <- CDM::data.ecpe$q.matrix[8:10, ]
Q_mat <- as.matrix(Q)
expected_pnames <- cdModel::test_parameter_names(Q_mat)
K <- ncol(Q)
sigma <- 0
means <- rep(0, ncol(Q))
Sigma <- cdModel::build_correlation_matrix(K, sigma)

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
    constants = NULL,
    # constants = c(
    #     pi_000 = skill_probs$probability[1],
    #     pi_100 = skill_probs$probability[2],
    #     pi_010 = skill_probs$probability[3],
    #     pi_110 = skill_probs$probability[4],
    #     pi_001 = skill_probs$probability[5],
    #     pi_101 = skill_probs$probability[6],
    #     pi_011 = skill_probs$probability[7]
    # ),
    match_map = NULL,
    accumulators = Q_mat,
    type = "cdm",
    verbose = TRUE
)

# tibble::as_tibble(CDM::data.ecpe$data)

dat <- CDM::data.ecpe$data[, c("id", paste0("E", 8:10))]
dat_long <- wide2long_base(dat)
dat_long$s <- 1
dat_responses <- dat_long[, 2:5]
# tibble::as_tibble(dat_long)
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
plot_prior(p_prior)

sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = FALSE
)
sim_p_vector <- c(
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
N <- 3000
dat <- simulate(sub_model,
    nsim = N,
    parameter_vector = sim_p_vector,
    nschool = 1,
    debug = FALSE,
    seed = 123
)

# head(dat$responses)
str(dat$responses)

str(dat_responses)
dat_responses$student <- factor(dat_responses$student)
dat_responses$item <- factor(dat_responses$item)
dat_responses$s <- factor(dat_responses$s)

# sub_dmis <- BuildDMI(dat$responses, model,
#     q_matrix = model@cdm_info$q_matrix,
#     rule = "DINA",
#     use_mvn = FALSE
# )


unique(dat_responses$student)

sub_dmis <- BuildDMI(dat_responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = FALSE
)
sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = model@cdm_info$q_matrix,
    rule = "DINA",
    use_mvn = FALSE
)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.00, thin = 1, is_pblocked = FALSE
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 4, is_pblocked = FALSE
)


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
head(ggdmc::summary)
est_theta <- ggdmc::summary(fit)
est_theta$quantiles
est_theta$statistics
# p0 <- ggdmc::plot(fit[[1]], facet_chains = F, pll = T)

p0 <- ggdmc::plot(fits0[[1]], start = 1)
p0 <- ggdmc::plot(fits0[[2]], start = 1)
p0 <- ggdmc::plot(fits0[[3]], start = 1)

p0 <- ggdmc::plot(fits0[[1]], start = 200)
p0 <- ggdmc::plot(fits0[[2]], start = 200)
p0 <- ggdmc::plot(fits0[[3]], start = 30)

p0 <- ggdmc::plot(fits0[[1]], pll = FALSE, den = TRUE, start = fits0[[1]]@nmc * 0.5)

p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
p0 <- ggdmc::plot(fit,
    pll = FALSE, den = TRUE, start = fit@nmc * 0.5,
    hide_legend = FALSE
)



dat <- CDM::data.ecpe$data[, -1]
Q <- CDM::data.ecpe$q.matrix
ecpe <- CDM::din(data = dat, q.matrix = Q)

param <- CDM::IRT.se(ecpe, extended = TRUE)
p <- split(param, param$partype)
cdm_result <- tibble::tibble(
    guess = p$guess$est,
    slip = p$slip$est,
    g_se = p$guess$se,
    s_se = p$slip$se
)

p$probs$est
p$probs$se
