# q(save = "no")
cat("\n\n-------------------- 6 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmcDE", "ggdmcModel", "ggdmcPrior", "ggdmcPhi", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

fn <- "~/Documents/ggdmcDE/tests/testthat/Group1/data/ggdmc_data1.rda"
fg_path <- "~/Documents/ggdmcDE/tests/testthat/Group5/data/6param.pdf"
save_path <- "~/Documents/ggdmcDE/tests/testthat/Group5/data/6param.rda"

load(fn)
load(save_path)


# fits0 <- StartSampling(pop_dmis, pop_priors, sub_migration_prob = 0.06, thin = 1, seed = 9032)
# save(fits0, file = save_path)

# fits1 <- RestartSampling(fits0, sub_migration_prob = 0.02, thin = 1, seed = 9032)
# save(fits0, fits1, file = save_path)

# fits2 <- RestartSampling(fits1, sub_migration_prob = 0.04, thin = 2, seed = 9032)
# save(fits0, fits1, fits2, file = save_path)

# fits3 <- RestartSampling(fits2, pop_migration_prob = 0.00, sub_migration_prob = 0, thin = 2, seed = 9032)
# save(fits0, fits1, fits2, fits3,
#     file = save_path
# )
fits4 <- RestartSampling(fits3,
    pop_migration_prob = 0.05, sub_migration_prob = 0.00,
    thin = 8, seed = 9032
)
save(fits0, fits1, fits2, fits3, fits4,
    file = save_path
)
# fits5 <- RestartSampling(fits4,
#     pop_migration_prob = 0.05, sub_migration_prob = 0.00,
#     thin = 4, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, fits4, fits5,
#     file = save_path
# )
# fits6 <- RestartSampling(fits5,
#     pop_migration_prob = 0.00, sub_migration_prob = 0,
#     thin = 4, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, fits4, fits5, fits6,
#     file = save_path
# )

# fits7 <- RestartSampling(fits6,
#     pop_migration_prob = 0.00, sub_migration_prob = 0,
#     thin = 8, seed = 9032
# )
# save(fits0, fits1, fits2, fits3, fits4, fits5, fits6, fits7,
#     file = save_path
# )


# fits8 <- RestartSampling(fits7,
#     is_add = TRUE, nmc = 100, thin = 1,
#     pop_migration_prob = 0.00, sub_migration_prob = 0,
#     seed = 9032
# )
# save(fits0, fits1, fits2, fits3, fits4, fits5, fits6, fits7, fits8,
#     file = save_path
# )


fits <- fits4
phi <- ggdmcDE::RebuildHyper(fits)
thetas <- ggdmcDE:::RebuildPosteriors(fits)
options(digits = 2)
est_phi <- ggdmcPhi::compare(phi, ps = true_vector)

# est_phi <- ggdmcPhi::compare(fits[[1]]$phi, ps = true_vector)

# fits[[1]]$phi@theta[1, 1, 1:3]
# fits[[1]]$phi@theta[1, 1, 501:503]


# fits <- fits6

# fit1 <- fits[[1]]
# fit2 <- fits[[2]]
# fit3 <- fits[[3]]
# phi1 <- fit1$phi
# phi2 <- fit2$phi
# phi3 <- fit3$phi
# thetas1 <- fit1$subject_theta
# thetas2 <- fit2$subject_theta
# thetas3 <- fit3$subject_theta


# phi <- RebuildHyper(fits)
# thetas <- ggdmcDE:::RebuildPosteriors(fits)

# est_phi <- gdmc:::summary(phi, start = 1, recovery = TRUE, ps = true_vector, verbose = TRUE)
# est_theta <- gdmc:::summary(thetas, recovery = TRUE, ps = ps, verbose = TRUE)


# pdf(fg_path)
# p0 <- gdmc:::plot(thetas, start = 1)
# p1 <- gdmc:::plot(phi, den = T, pll = F)
# dev.off()

# gdmc:::plot(thetas1, start = 1)
# gdmc:::plot(thetas2, start = 1)
# gdmc:::plot(thetas3, start = 1)

# gdmc:::plot(phi1, den = T, pll = F)
# gdmc:::plot(phi2, den = T, pll = F)
# gdmc:::plot(phi3, den = T, pll = F)
