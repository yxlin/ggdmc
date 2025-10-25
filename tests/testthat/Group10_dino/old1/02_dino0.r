# q(save = "no")
cat("\n\n-------------------- Testing DINO model0 --------------------")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dino0.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)

# StartSampling_subject
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 4,
    is_pblocked = FALSE,
    seed = 9032
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 4,
    seed = 9032
)
fits <- fits1


fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 3)
est_theta <- ggdmc::compare(fit, ps = p_vector)

#               guess1  guess2   guess3    pi1     pi2    pi3  slip1  slip2
# True          0.1000 0.20000  0.30000 0.1000  0.3000  0.500 0.0200 0.0500
#    5 Estimate 0.0433 0.11840  0.00999 0.1406  0.1547  0.138 0.0103 0.0109
#   50 Estimate 0.1315 0.20229  0.12066 0.1826  0.2441  0.219 0.1004 0.1049
# 97.5 Estimate 0.2402 0.30363  0.32462 0.2539  0.3459  0.312 0.2298 0.2325
# Median-True   0.0315 0.00229 -0.17934 0.0826 -0.0559 -0.281 0.0804 0.0549
#                 slip3
# True           0.0700
#    5 Estimate  0.0244
#   50 Estimate  0.0577
# 97.5 Estimate  0.0936
# Median-True   -0.0123


#  partype parindex parameter    est     se  2.5 % 97.5 % item item.name
# 2    slip        2   E1_slip 0.101 0.0426 0.0176 0.1847    1        E1
# 4    slip        4   E2_slip 0.111 0.0411 0.0303 0.1915    2        E2
# 6    slip        6   E3_slip 0.052 0.0188 0.0152 0.0887    3        E3
#   skillclass fixed free rule totindex
# 2          0 FALSE TRUE DINO        2
# 4          0 FALSE TRUE DINO        4
# 6          0 FALSE TRUE DINO        6

#  partype parindex   parameter     est     se  2.5 %  97.5% item item.name
# 7    probs        7 prob_class1 0.190 0.0158 0.159  0.221    0
# 8    probs        8 prob_class2 0.243 0.0203 0.203  0.283    0
# 9    probs        9 prob_class3 0.224 0.0191 0.187  0.261    0
# 10   probs        0 prob_class4 0.343 0.0286 0.287  0.400    0
#    skillclass fixed  free rule totindex
# 7           1 FALSE  TRUE             7
# 8           2 FALSE  TRUE             8
# 9           3 FALSE  TRUE             9
# 10          4 FALSE FALSE            10

#    partype parindex parameter  est     se   2.5%   97.5% item item.name
# 1   guess        1  E1_guess 0.140 0.0532  0.0361  0.244    1        E1
# 3   guess        3  E2_guess 0.207 0.0488  0.1113  0.302    2        E2
# 5   guess        5  E3_guess 0.136 0.0924 -0.0449  0.317    3        E3
#   skillclass fixed free rule totindex
# 1          0 FALSE TRUE DINO        1
# 3          0 FALSE TRUE DINO        3
# 5          0 FALSE TRUE DINO        5
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


# compare2ecpe <- CDM::din(data = data.frame(wide[, -1]), q.matrix = Q, rule = "DINO")
# compared_param <- CDM::IRT.se(compare2ecpe, extended = TRUE)
# compared_p <- split(compared_param, compared_param$partype)
# compared_p$guess
# compared_p$probs
# compared_p$slip

#                guess1  guess2 guess3  slip1  slip2   slip3
# True            0.100 0.20000  0.300 0.0100 0.0200  0.0300
# 5% Estimate     0.114 0.15293  0.196 0.0021 0.0034  0.0018
# 50% Estimate    0.161 0.20069  0.271 0.0228 0.0318  0.0131
# 97.5% Estimate  0.219 0.26507  0.353 0.0796 0.0972  0.0346
# Median-True     0.061 0.00069 -0.029 0.0128 0.0118 -0.0169



# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, pll = F, start = fit@nmc * 0.5)

# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, start = fits[[1]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, start = fits[[2]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, start = fits[[3]]@nmc * 0.5)


# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, den = TRUE, start = fits[[1]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, den = TRUE, start = fits[[2]]@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, den = TRUE, start = fits[[3]]@nmc * 0.5)
