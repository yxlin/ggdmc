# q(save = "no")
cat("\n\n-------------------- v model --------------------")
rm(list = ls())
pkg <- c("ggdmcDE", "ggdmcModel", "ggdmcPrior", "ggdmcPhi", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")


# fn <- "tests/testthat/Group1/data/ggdmc_data2.rda"
# fg_path <- "tests/testthat/Group5/data/6param.pdf"
# save_path <- "tests/testthat/Group5/data/6param.rda"

fn <- "data/ggdmc_data2.rda"
fg_path <- "data/v_model.pdf"
save_path <- "data/v_model.rda"

load(fn)

model <- ggdmcModel::BuildModel(
    p_map = list(A = "1", B = "1", mean_v = c("POLITICAL_VIEW", "M"), sd_v = "M", st0 = "1", t0 = "1"),
    match_map = list(M = list(s1 = "r1", s2 = "r2")),
    factors = list(S = c("s1", "s2"), POLITICAL_VIEW = c("liberal", "conservative")),
    constants = c(sd_v.false = 1, st0 = 0),
    accumulators = c("r1", "r2"),
    type = "lba"
)


# fits0 <- StartSampling(pop_dmis, pop_priors, sub_migration_prob = 0.06, thin = 8, seed = 9032)
# save(fits0, file = save_path)

# fits1 <- RestartSampling(fits0, sub_migration_prob = 0.02, thin = 4, seed = 9032)
# fits2 <- RestartSampling(fits1, pop_migration_prob = 0.05, sub_migration_prob = 0, seed = 9032)
# save(fits0, fits1, fits2,
#     file = save_path
# )
load(save_path)
fits3 <- RestartSampling(fits2, pop_migration_prob = 0.00, sub_migration_prob = 0.05, thin = 4, seed = 9032)
save(fits0, fits1, fits2, fits3,
    file = save_path
)
fits4 <- RestartSampling(fits3, pop_migration_prob = 0.00, sub_migration_prob = 0.00, thin = 4, seed = 9032)

save(fits0, fits1, fits2, fits3, fits4,
    file = save_path
)

fits <- fits4

fit1 <- fits[[1]]
fit2 <- fits[[2]]
fit3 <- fits[[3]]
phi1 <- fit1$phi
phi2 <- fit2$phi
phi3 <- fit3$phi
thetas1 <- fit1$subject_theta
thetas2 <- fit2$subject_theta
thetas3 <- fit3$subject_theta


phi <- RebuildHyper(fits)
thetas <- ggdmcDE:::RebuildPosteriors(fits)

est_phi <- gdmc:::summary(phi, start = 1, recovery = TRUE, ps = true_vector, verbose = TRUE)
est_theta <- gdmc:::summary(thetas, recovery = TRUE, ps = ps, verbose = TRUE)


pdf(fg_path)
p0 <- gdmc:::plot(thetas, start = 1)
p1 <- gdmc:::plot(phi, den = T, pll = F)
dev.off()

# gdmc:::plot(thetas1, start = 1)
# gdmc:::plot(thetas2, start = 1)

# gdmc:::plot(phi1, den = T, pll = F)
# gdmc:::plot(phi2, den = T, pll = F)
# gdmc:::plot(phi3, den = T, pll = F)
