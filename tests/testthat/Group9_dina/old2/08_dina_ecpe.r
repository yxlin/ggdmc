# q(save = "no")
cat("\n\n----------------- Testing model ecpe --------------------\n")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina_ecpe.rda")
save_path <- file.path(home_dir, "Group10_cdm_subjects/data/dina_ecpe.rda")

source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)
load(save_path)

# 44 min for 500
# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     nmc = 10, report_length = 2,
#     sub_migration_prob = 0.03, thin = 1, is_pblocked = TRUE,
#     seed = 9032
# )

# (53.12 * 10 * 5) / 60
# save(fits0, file = save_path)

# (1.85 * 10 * 5) / 60
# 1.54 min
fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 4
)

save(fits0, file = save_path)

# fits1 <- RestartSampling_subject(fits0,
#     nmc = 300, report_length = 50,
#     sub_migration_prob = 0.03, thin = 2, is_pblocked = TRUE,
#     seed = 9032
# )
# save(fits0, fits1, file = save_path)

# fits1 <- RestartSampling_subject(fits0,
#     sub_migration_prob = 0.00, thin = 2,
#     seed = 9032
# )

# save(fits0, fits1, file = save_path)

fits <- fits0
fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)

p0 <- ggdmc::plot(fits0[[1]], start = fits0[[1]]@nmc * 0.5)
p0 <- ggdmc::plot(fits0[[2]], start = fits0[[2]]@nmc * 0.5)
p0 <- ggdmc::plot(fits0[[3]], start = fits0[[3]]@nmc * 0.5)
p1 <- plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# pdf("Rplot.pdf")

p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, start = 1)

# p1 <- ggdmc::plot(fit, facet_chains = F, pll = F, start = fit@nmc * 0.5)
