# q(save = "no")
cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
save_path <- paste0(home_dir, "Group9_gen_cdm/data/dina_ecpe.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
# tibble::as_tibble(data.ecpe$data)
# tibble::as_tibble(data.ecpe$data)
Q <- as.matrix(CDM::data.ecpe$q.matrix)
# Q <- matrix(c(
#     1, 0,
#     0, 1,
#     1, 1
# ), ncol = 2, byrow = TRUE)
# str(Q)
# colnames(Q)
# rownames(Q)

# colnames(Q) <- c("Algebra", "Geometry")
# rownames(Q) <- c("Item1", "Item2", "Item3")

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        guess4 = "1", guess5 = "1", guess6 = "1",
        guess7 = "1", guess8 = "1", guess9 = "1",
        guess10 = "1", guess11 = "1", guess12 = "1",
        guess13 = "1", guess14 = "1", guess15 = "1",
        guess16 = "1", guess17 = "1", guess18 = "1",
        guess19 = "1", guess20 = "1", guess21 = "1",
        guess22 = "1", guess23 = "1", guess24 = "1",
        guess25 = "1", guess26 = "1", guess27 = "1",
        guess28 = "1",
        pi1 = "1", pi2 = "1", pi3 = "1",
        pi4 = "1", pi5 = "1", pi6 = "1",
        pi7 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1",
        slip4 = "1", slip5 = "1", slip6 = "1",
        slip7 = "1", slip8 = "1", slip9 = "1",
        slip10 = "1", slip11 = "1", slip12 = "1",
        slip13 = "1", slip14 = "1", slip15 = "1",
        slip16 = "1", slip17 = "1", slip18 = "1",
        slip19 = "1", slip20 = "1", slip21 = "1",
        slip22 = "1", slip23 = "1", slip24 = "1",
        slip25 = "1", slip26 = "1", slip27 = "1",
        slip28 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = Q,
    type = "cdm",
    verbose = TRUE
)


raw_pop_mean <- c(
    guess1 = 0.7053341, guess2 = 0.7238056, guess3 = 0.4381031,
    guess4 = 0.4804197, guess5 = 0.7636764, guess6 = 0.7173405, guess7 = 0.5438179,
    guess8 = 0.8017958, guess9 = 0.5343676, guess10 = 0.4827295,
    guess11 = 0.5562213, guess12 = 0.1946139, guess13 = 0.6330542, guess14 = 0.5166922,
    guess15 = 0.7489095, guess16 = 0.5492260, guess17 = 0.8155547,
    guess18 = 0.7293164, guess19 = 0.4732080, guess20 = 0.2388085, guess21 = 0.6214631,
    guess22 = 0.3217567, guess23 = 0.6370685, guess24 = 0.3134852,
    guess25 = 0.5117201, guess26 = 0.5549666, guess27 = 0.2650583, guess28 = 0.6591149,
    pi1 = 0.311074178, pi2 = 0.006133434, pi3 = 0.040339807, pi4 = 0.049576658,
    pi5 = 0.012468516, pi6 = 0.025885034, pi7 = 0.103326045,
    slip1 = 0.08503517, slip2 = 0.10092311, slip3 = 0.26573233, slip4 = 0.16201356,
    slip5 = 0.04046587, slip6 = 0.06648817, slip7 = 0.08462383, slip8 = 0.03977443,
    slip9 = 0.19913180, slip10 = 0.16276272, slip11 = 0.09893619, slip12 = 0.30515493,
    slip13 = 0.12169667, slip14 = 0.21181583, slip15 = 0.03994119, slip16 = 0.12570135,
    slip17 = 0.05795799, slip18 = 0.08601661, slip19 = 0.15017251, slip20 = 0.29549079,
    slip21 = 0.09655996, slip22 = 0.18779334, slip23 = 0.07527230, slip24 = 0.32193191,
    slip25 = 0.27165688, slip26 = 0.21068721, slip27 = 0.36867236, slip28 = 0.08607046
)

raw_pop_scale <- c(
    guess1 = 0.01256021, guess2 = 0.01557211, guess3 = 0.01380687,
    guess4 = 0.01738231, guess5 = 0.01391434,
    guess6 = 0.01507642,
    guess7 = 0.01399755,
    guess8 = 0.01345476,
    guess9 = 0.01714521,
    guess10 = 0.01443041,
    guess11 = 0.01372773,
    guess12 = 0.01167786,
    guess13 = 0.01365902, guess14 = 0.01427147, guess15 = 0.01441244, guess16 = 0.01391425,
    guess17 = 0.01200980, guess18 = 0.01488188,
    guess19 = 0.01759947,
    guess20 = 0.01232369,
    guess21 = 0.01346624,
    guess22 = 0.01702903,
    guess23 = 0.01691710,
    guess24 = 0.01769723,
    guess25 = 0.01427662,
    guess26 = 0.01703229,
    guess27 = 0.01303236,
    guess28 = 0.01621625,
    pi1 = 0.016615891, pi2 = 0.008651352, pi3 = 0.012668820, pi4 = 0.010517344,
    pi5 = 0.007852105, pi6 = 0.007090650, pi7 = 0.011238043,
    slip1 = 0.009349528, slip2 = 0.008886350, slip3 = 0.013341317, slip4 = 0.009833474,
    slip5 = 0.005295011, slip6 = 0.006745318, slip7 = 0.009047344, slip8 = 0.006163573,
    slip9 = 0.010528718, slip10 = 0.011378382, slip11 = 0.009481677, slip12 = 0.014057252,
    slip13 = 0.009910174, slip14 = 0.012272701, slip15 = 0.005382480, slip16 = 0.010455328,
    slip17 = 0.007078148, slip18 = 0.007365173, slip19 = 0.009597655, slip20 = 0.013927080,
    slip21 = 0.009180179, slip22 = 0.010552648, slip23 = 0.008456843, slip24 = 0.013412625,
    slip25 = 0.013371752, slip26 = 0.010515228, slip27 = 0.014472118, slip28 = 0.007426670
)

pop_mean <- round(raw_pop_mean, 2)
pop_scale <- round(raw_pop_scale, 2)



pop_dist <- ggdmcPrior::BuildPrior(
    p0 = pop_mean,
    p1 = pop_scale,
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dists = rep("tnorm", model@npar),
    log_p = rep(F, model@npar)
)

# plot_prior(pop_dist)

# 1. The probability of occurrence in the population, written as P(alpha_l)
# 2. profile probability, gamma (𝛾) = [P(alpha_1), P(alpha_2), ..., P(alpha_{2^K})].
sub_model <- setCDM(model,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINA"
)

pop_model <- setCDM(model,
    population_distribution = pop_dist,
    q_matrix = model@cdm_info$q_matrix,
    profile_probability = model@cdm_info$profile_probability, rule = "DINA"
)

p_vector <- pop_mean


N_total_student <- 2922
dat <- simulate(sub_model,
    nsim = N_total_student, parameter_vector = p_vector,
    nschool = 1,
    seed = 123,
    debug = FALSE
)

hdat <- simulate(pop_model,
    nsim = N_total_student,
    nschool = 32,
    seed = 123
)

names(pop_dist)
model@pnames == names(pop_dist)

ps <- attr(hdat, "parameters")


wide <- long2wide_base(dat$response)
ecpe2 <- CDM::din(data = wide[, -1], q.matrix = model@cdm_info$q_matrix)


tibble::as_tibble(dat$response)
tibble::as_tibble(wide)
# head(hdat$response)
# table(hdat$response$s)
# table(dat$response$s)

sub_dmis <- BuildDMI(dat$responses, model,
    q_matrix = sub_model@q_matrix,
    profile_probability = sub_model@profile_probability, rule = "DINA"
)

pop_dmis <- BuildDMI(hdat$responses, model,
    q_matrix = sub_model@q_matrix,
    profile_probability = sub_model@profile_probability, rule = "DINA"
)

true_mean <- pop_mean[sort(names(pop_mean))]
true_scale <- pop_scale[sort(names(pop_scale))]
names(true_mean) <- paste0("loc_", names(true_mean))
names(true_scale) <- paste0("sca_", names(true_scale))
true_vector <- c(true_mean, true_scale)


# Fit subject/school theta starting samples -------------------------
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
p_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(1.1, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

sub_priors <- set_priors(p_prior = p_prior)


nmc <- 500
sub_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = model@pnames)
sub_samples <- ggdmc::initialise_theta(sub_theta_input, sub_priors, sub_dmis[[1]],
    verbose = F
)


# Generate pop samples -------------------------
p0 <- runif(model@npar)
names(p0) <- model@pnames
model_likelihood <- ggdmcPrior::BuildPrior(
    p0 = p0,
    # 1. Setting p1 to 10 in the hyper lileihood function, rendering
    # the model hard to converge;
    # 2. TODO: test if increase the sample size will resolve this problem
    # p1 = rep(10, model@npar),
    # 3. We know that the probability is within 0 and 1, bcz this
    # is what CDM defines it item parameters. They are probabilities.
    p1 = rep(1.1, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("tnorm", model@npar),
    log_p = rep(TRUE, model@npar)
)

# Prior log likelihoods
p0 <- rep(0, model@npar)
names(p0) <- model@pnames
l_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(0, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)
s_prior <- ggdmcPrior::BuildPrior(
    p0 = p0,
    p1 = rep(10, model@npar),
    lower = rep(NA, model@npar),
    upper = rep(NA, model@npar),
    dist = rep("unif", model@npar),
    log_p = rep(TRUE, model@npar)
)

pop_priors <- ggdmcPrior::set_priors(
    p_prior = model_likelihood, l_prior = l_prior,
    s_prior = s_prior
)
pop_theta_input <- ggdmc::setThetaInput(nmc = nmc, pnames = pop_priors@pnames)

pop_samples <- ggdmc::initialise_phi(pop_theta_input, pop_priors, pop_dmis,
    seed = 846671, verbose = FALSE
)

save_path

# save(model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
#     sub_dmis, pop_dmis,
#     sub_priors, sub_samples, sub_theta_input,
#     pop_priors, pop_samples, pop_theta_input, Q,
#     file = save_path
# )

save(model, hdat, dat, p_vector, pop_mean, pop_scale, true_vector, ps,
    sub_dmis, pop_dmis,
    sub_priors, sub_samples, sub_theta_input,
    pop_priors, pop_theta_input, Q,
    file = save_path
)
