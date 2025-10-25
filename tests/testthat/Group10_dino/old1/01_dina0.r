#!/usr/bin/env Rscript
# q(save = "no")
cat("\n\n----------------- Testing model0 --------------------\n")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina0.rda")
save_path <- file.path(home_dir, "Group10_cdm_subjects/data/dina0.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)
# load(save_path)
cat(paste0("Load data file from ", data_path), "\n")

# StartSampling_subject
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.01, thin = 8, is_pblocked = FALSE,
    seed = 9032
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2,
    seed = 9032
)
fits <- fits1

# cat(paste0("Save fit results to ", save_path), "\n")
# save(fits0, fits1, file = save_path)

fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)

p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = 1)
#               guess1  guess2  guess3   pi1     pi2    pi3   slip1  slip2
# True          0.1000  0.2000  0.3000 0.100  0.3000  0.500 0.02000 0.0500
#    5 Estimate 0.0135  0.0164  0.1983 0.132  0.1484  0.168 0.00784 0.0191
#   50 Estimate 0.1392  0.1512  0.2691 0.197  0.2178  0.243 0.06456 0.0901
# 97.5 Estimate 0.3042  0.3272  0.3266 0.273  0.3145  0.342 0.16546 0.1815
# Median-True   0.0392 -0.0488 -0.0309 0.097 -0.0822 -0.257 0.04456 0.0401
#                  slip3
# True          7.00e-02
#    5 Estimate 6.77e-03
#   50 Estimate 7.00e-02
# 97.5 Estimate 1.93e-01
# Median-True   2.77e-06


#               guess1  guess2  guess3   pi1     pi2    pi3   slip1  slip2
# True          0.1000  0.2000  0.3000 0.100  0.3000  0.500 0.02000 0.0500
#    5 Estimate 0.0152  0.0170  0.2059 0.132  0.1462  0.172 0.00689 0.0249
#   50 Estimate 0.1460  0.1605  0.2733 0.200  0.2177  0.248 0.05871 0.0987
# 97.5 Estimate 0.3098  0.3336  0.3282 0.278  0.3151  0.349 0.16282 0.1889
# Median-True   0.0460 -0.0395 -0.0267 0.100 -0.0823 -0.252 0.03871 0.0487
#                 slip3
# True          0.07000
#    5 Estimate 0.00666
#   50 Estimate 0.07139
# 97.5 Estimate 0.19723
# Median-True   0.00139

#
#   partype parindex parameter    est     se   2.5 % 97.5 % item item.name
# 2    slip        2   E1_slip 0.0776 0.0594 -0.0387  0.194    1        E1
# 4    slip        4   E2_slip 0.1136 0.0507  0.0142  0.213    2        E2
# 6    slip        6   E3_slip 0.0651 0.0866 -0.1046  0.235    3        E3
#   skillclass fixed free rule totindex
# 2          0 FALSE TRUE DINA        2
# 4          0 FALSE TRUE DINA        4
# 6          0 FALSE TRUE DINA        6
#
#     partype parindex   parameter  est     se  2.5%  97.5% item item.name
# 7    probs        7 prob_class1 0.194 0.0203 0.155  0.234    0
# 8    probs        8 prob_class2 0.222 0.0196 0.183  0.260    0
# 9    probs        9 prob_class3 0.248 0.0210 0.207  0.289    0
# 10   probs        0 prob_class4 0.336 0.0274 0.283  0.390    0
#    skillclass fixed  free rule totindex
# 7           1 FALSE  TRUE             7
# 8           2 FALSE  TRUE             8
# 9           3 FALSE  TRUE             9
# 10          4 FALSE FALSE            10


#   partype parindex parameter  est    se   2.5%  97.5% item item.name skillclass
# 1   guess        1  E1_guess 0.16 0.079 0.0038   0.31    1        E1          0
# 3   guess        3  E2_guess 0.18 0.075 0.0363   0.33    2        E2          0
# 5   guess        5  E3_guess 0.26 0.045 0.1681   0.34    3        E3          0
#   fixed free rule totindex
# 1 FALSE TRUE DINA        1
# 3 FALSE TRUE DINA        3
# 5 FALSE TRUE DINA        5

#                guess1 guess2 guess3  slip1  slip2   slip3
# True            0.100  0.200  0.300 0.0100 0.0200  0.0300
# 5% Estimate     0.095  0.123  0.246 0.0013 0.0017  0.0013
# 50% Estimate    0.147  0.175  0.275 0.0165 0.0215  0.0166
# 97.5% Estimate  0.212  0.241  0.312 0.0640 0.0759  0.0693
# Median-True     0.047 -0.025 -0.025 0.0065 0.0015 -0.0134



# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = 1)
# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, start = 14)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, pll = F, start = fit@nmc * 0.5)

p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, start = fits[[1]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, start = fits[[2]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, start = fits[[3]]@nmc * 0.5)

p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, start = 10)
p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, start = 1)
p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, start = 1)


# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, den = TRUE, start = fits[[1]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, den = TRUE, start = fits[[2]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, den = TRUE, start = fits[[3]]@nmc * 0.5)


# long <- tibble::as_tibble(dat$responses)
# wide <- long2wide_base(long)
# tibble::as_tibble(wide)

# Q <- matrix(c(
#     1, 0,
#     0, 1,
#     1, 1
# ), ncol = 2, byrow = TRUE)

# colnames(Q) <- c("Algebra", "Geometry")
# rownames(Q) <- c("Item1", "Item2", "Item3")


# compare2ecpe <- CDM::din(data = data.frame(wide[, -1]), q.matrix = Q)
# compared_param <- CDM::IRT.se(compare2ecpe, extended = TRUE)
# compared_p <- split(compared_param, compared_param$partype)
# compared_p$guess
# compared_p$probs
# compared_p$slip

# c_omega1 <- 1 - compared_p$guess$est - compared_p$slip$est # item discrimination
# c_omega2 <- (compared_p$guess$est + (1 - compared_p$slip$est)) / 2 # item easiness

# c_pvalues <- colMeans(wide[, -1], na.rm = TRUE)
