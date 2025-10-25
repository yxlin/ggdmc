# q(save = "no")
cat("\n\n-------------------- Generate RRUM model 0 --------------------")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
data_path <- paste0(home_dir, "Group9_gen_cdm/data/rrum0.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)


# StartSampling_subject
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.02, thin = 1,
    # is_pblocked = TRUE,
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

# ps <- p_vector
# ps
options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)

#                 pi1     pi2   pi3 pi_item1 pi_item2 pi_item3  r11   r22   r31
# True           0.10 0.20000 0.300    0.100    0.200    0.300 0.10 0.400 0.300
#    5 Estimate  0.05 0.15467 0.206    0.021    0.029    0.032 0.07 0.078 0.062
#   50 Estimate  0.07 0.20075 0.316    0.209    0.269    0.251 0.52 0.543 0.514
# 97.5 Estimate  0.13 0.35746 0.670    0.616    0.677    0.707 0.98 0.965 0.951
# Median-True   -0.03 0.00075 0.016    0.109    0.069   -0.049 0.42 0.143 0.214
#                  r32
# True           0.600
#    5 Estimate  0.086
#   50 Estimate  0.567
# 97.5 Estimate  0.973
# Median-True   -0.033

# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# long <- tibble::as_tibble(dat$responses)
# wide <- long2wide_base(long)
# mod1 <- CDM::gdina(data = data.frame(wide[, -1]), q.matrix = Q, rule = "RRUM")
# summary(mod1)

# RRUM Parametrization
#      pi r_Algebra r_Geometry
# E1 0.14      0.17
# E2 0.33                 0.25
# E3 0.20      0.89       0.98

# Skill Pattern Probabilities

#   00   10   01   11
# 0.50 0.21 0.26 0.03
