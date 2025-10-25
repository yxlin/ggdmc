# q(save = "no")
cat("\n\n----------------- Testing model9 --------------------\n")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina9.rda")
save_path <- file.path(home_dir, "Group10_cdm_subjects/data/dina9.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))

cat(paste0("Load data file from ", data_path), "\n")


load(data_path)
# load(save_path)

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 10, is_pblocked = TRUE,
    seed = 9032
)

save(fits0, file = save_path)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2, is_pblocked = TRUE,
    seed = 9032
)

save(fits0, fits1, file = save_path)

# fits2 <- ggdmc:::RestartSampling_subject(fits1,
#     sub_migration_prob = 0.00, thin = 2, is_pblocked = TRUE,
#     seed = 9032
# )
# save(fits0, fits1, fits2, file = save_path)
cat(paste0("Save fit results to ", save_path), "\n")


fits <- fits1
fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")
options(digits = 3)
est_theta <- ggdmc::compare(fit, ps = p_vector)

pdf("Dina9_Rplot.pdf")
p0 <- ggdmc::plot(fit)
p1 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = 1)
dev.off()
