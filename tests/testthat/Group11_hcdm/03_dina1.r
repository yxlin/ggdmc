# q(save = "no")
cat("\n\n-------------------- Generate model 0 --------------------")
rm(list = ls())

pkg <- c("ggdmc")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
data_path <- paste0(data_dir, "Group9_gen_cdm/data/dina1.rda")
save_path <- file.path(home_dir, "Group11_hcdm/data/dina1.rda")


cat("\nLoad data: ", data_path, "\n")
load(data_path)

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


fits <- fits2
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

options(digits = 2)
est_phi <- compare(phi, ps = true_vector)

#    loc_guess1 loc_guess2 loc_guess3 loc_slip1 loc_slip2 loc_slip3
# True               0.1000      0.200      0.300    0.0100    0.0300    0.0500
# 5% Estimate        0.0066      0.013      0.062    0.0044    0.0041    0.0054
# 50% Estimate       0.0837      0.146      0.386    0.0565    0.0546    0.0682
# 97.5% Estimate     2.6380      0.924      1.075    1.3122    1.7953    1.0408
# Median-True       -0.0163     -0.054      0.086    0.0465    0.0246    0.0182
#                sca_guess1 sca_guess2 sca_guess3 sca_slip1 sca_slip2 sca_slip3
# True                 0.01       0.03       0.05      0.05      0.07      0.09
# 5% Estimate          0.39       0.25       0.13      0.37      0.24      0.44
# 50% Estimate         0.52       0.41       0.18      0.49      0.36      0.58
# 97.5% Estimate       4.31       3.71       3.43      3.58      3.29      4.29
# Median-True          0.51       0.38       0.13      0.44      0.29      0.49


