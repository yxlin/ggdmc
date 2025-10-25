# q(save = "no")
cat("\n\n-------------------- Generate model 1 --------------------")
rm(list = ls())
pkg <- c("ggdmc")
sapply(pkg, require, character.only = TRUE)

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dino1.rda")
save_path <- file.path(home_dir, "Group11_hcdm/data/dino1.rda")
load(data_path)
# load(save_path)

# Fit random-effect models -------
fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.01, is_pblocked = TRUE,
    thin = 1L, pop_debug = F, seed = 9032
)
save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    sub_migration_prob = 0.00,
    is_pblocked = TRUE,
    thin = 1L, seed = 9032
)
save(fits0, fits1, file = save_path)

fits2 <- RestartSampling(fits1,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    is_pblocked = FALSE,
    is_hblocked = TRUE,
    thin = 1L, seed = 9032
)
save(fits0, fits1, fits2, file = save_path)


fits3 <- RestartSampling(fits2,
    pop_migration_prob = 0.01,
    sub_migration_prob = 0.00,
    is_pblocked = FALSE,
    is_hblocked = TRUE,
    thin = 1L, seed = 9032
)
save(fits0, fits1, fits2, fits3, file = save_path)


fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

options(digits = 2)
est_phi <- compare(phi, ps = true_vector)

#    loc_guess1 loc_guess2 loc_guess3 loc_guess4 loc_slip1 loc_slip2
# True                0.100      0.200      0.300      0.400    0.0100    0.0200
# 5% Estimate         0.009      0.019      0.032      0.042    0.0079    0.0057
# 50% Estimate        0.129      0.215      0.366      0.543    0.1123    0.0793
# 97.5% Estimate      0.666      0.703      0.857      0.884    0.6220    0.4462
# Median-True         0.029      0.015      0.066      0.143    0.1023    0.0593
#                loc_slip3 loc_slip4 sca_guess1 sca_guess2 sca_guess3 sca_guess4
# True               0.030    0.0400       0.01       0.01       0.01       0.01
# 5% Estimate        0.003    0.0053       0.47       0.31       0.34       0.25
# 50% Estimate       0.041    0.0713       0.76       0.64       0.63       0.49
# 97.5% Estimate     0.279    0.4025       3.39       3.79       3.39       3.73
# Median-True        0.011    0.0313       0.75       0.63       0.62       0.48
#                sca_slip1 sca_slip2 sca_slip3 sca_slip4
# True                0.05      0.03      0.02      0.01
# 5% Estimate         0.47      0.32      0.14      0.31
# 50% Estimate        0.72      0.55      0.34      0.53
# 97.5% Estimate      3.61      3.82      3.72      3.37
# Median-True         0.67      0.52      0.32      0.52


