# q(save = "no")
rm(list = ls())
cat("\n\n-------------------- A model --------------------")

pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
wkdir <- file.path(home_dir, "Group4_hlba/")

data_path <- file.path(home_dir, "Group0_gen_data/data/lba_data4.rda")
save_path <- file.path(home_dir, "Group4_hlba/data/A_model.rda")
cat("\nWorking directory: ", wkdir, "\n")
load(data_path)


fits0 <- StartSampling(pop_dmis, pop_priors,
    sub_migration_prob = 0.01, is_pblocked = TRUE,
    thin = 2L, pop_debug = F, seed = 9032
)

save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    sub_migration_prob = 0.00,
    is_pblocked = TRUE,
    thin = 2L, seed = 9032
)

save(fits0, fits1, file = save_path)

fits2 <- RestartSampling(fits1,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    is_pblocked = FALSE,
    is_hblocked = TRUE,
    thin = 2L, seed = 9032
)

save(fits0, fits1, fits2, file = save_path)

fits <- fits2


thetas <- RebuildPosteriors(fits)
fit <- RebuildHyper(fits)


hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est <- ggdmc::compare(fit, ps = true_vector)
