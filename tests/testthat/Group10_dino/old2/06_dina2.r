# q(save = "no")
cat("\n\n-------------------- Testing model2 parameters -------------\n")
rm(list = ls())
pkg <- c("ggdmc")
sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina2.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)


fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.02, thin = 2,
    is_pblocked = TRUE,
    seed = 9032
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2,
    seed = 9032
)

fits <- fits1


fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)
#               guess1 guess2  guess3 guess4 guess5 pi1.s1 pi1.s2    pi2   pi3
# True           0.100  0.200  0.3000  0.400 0.5000   0.10   0.50  0.250 0.150
#    5 Estimate  0.091  0.181  0.2661  0.351 0.4675   0.20   0.20  0.211 0.210
#   50 Estimate  0.120  0.221  0.2911  0.387 0.5065   0.23   0.23  0.236 0.234
# 97.5 Estimate  0.155  0.266  0.3199  0.425 0.5489   0.27   0.26  0.266 0.264
# Median-True    0.020  0.021 -0.0089 -0.013 0.0065   0.13  -0.27 -0.014 0.084
#                slip1 slip2   slip3   slip4   slip5
# True          0.0100 0.020  0.0300  0.0400  0.0500
#    5 Estimate 0.0045 0.021  0.0044  0.0016  0.0029
#   50 Estimate 0.0264 0.053  0.0277  0.0117  0.0161
# 97.5 Estimate 0.0608 0.090  0.0633  0.0301  0.0368
# Median-True   0.0164 0.033 -0.0023 -0.0283 -0.0339

#     guess1   guess2 guess3 guess4 guess5 pi1.s1 pi1.s2    pi2   pi3
# True           0.100  0.20000  0.300  0.400 0.5000   0.10   0.50 0.2500 0.150
#    5 Estimate  0.105  0.16477  0.292  0.342 0.4782   0.21   0.21 0.2321 0.201
#   50 Estimate  0.137  0.19916  0.315  0.377 0.5087   0.24   0.24 0.2537 0.223
# 97.5 Estimate  0.173  0.23750  0.342  0.418 0.5436   0.27   0.27 0.2818 0.249
# Median-True    0.037 -0.00084  0.015 -0.023 0.0087   0.14  -0.26 0.0037 0.073
#                slip1    slip2   slip3   slip4   slip5
# True          0.0100  0.02000  0.0300  0.0400  0.0500
#    5 Estimate 0.0044  0.00076  0.0016  0.0034  0.0015
#   50 Estimate 0.0271  0.00945  0.0148  0.0165  0.0108
# 97.5 Estimate 0.0617  0.03630  0.0453  0.0360  0.0299
# Median-True   0.0171 -0.01055 -0.0152 -0.0235 -0.0392



# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p2 <- ggdmc::plot(fit, facet_chains = F, pll = F, start = fit@nmc * 0.5)


# long <- tibble::as_tibble(dat$responses)
# wide <- long2wide(long)
# wide
# compare2ecpe <- CDM::din(
#     data = data.frame(wide[, -1]),
#     q.matrix = model@cdm_info$q_matrix
# )

# compared_param <- CDM::IRT.se(compare2ecpe, extended = TRUE)
# compared_p <- split(compared_param, compared_param$partype)
# compared_p




# c_omega1 <- 1 - compared_p$guess$est - compared_p$slip$est # item discrimination
# c_omega2 <- (compared_p$guess$est + (1 - compared_p$slip$est)) / 2 # item easiness

# c_pvalues <- colMeans(wide[, -1], na.rm = TRUE)
# c_pvalues
