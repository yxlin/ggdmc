#!/usr/bin/env Rscript
# q(save = "no")
cat("\n\n----------------- Testing model0 --------------------\n")
rm(list = ls())
pkg <- c("ggdmc")
pkg_ok <- sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina12.rda")
save_path <- file.path(home_dir, "Group10_cdm_subjects/data/dina12.rda")
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))
load(data_path)
load(save_path)
cat(paste0("Load data file from ", data_path), "\n")
# sub_dmis[[1]]@use_mvn

# StartSampling_subject
# fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#     nmc = 3, ncore = 1,
#     sub_migration_prob = 0.00, thin = 1, is_pblocked = F,
#     seed = 9032
# )

fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
    sub_migration_prob = 0.05, thin = 2, is_pblocked = TRUE,
    seed = 9032
)

fits1 <- ggdmc:::RestartSampling_subject(fits0,
    sub_migration_prob = 0.00, thin = 2,
    seed = 9032
)

fits <- fits1

cat(paste0("Save fit results to ", save_path), "\n")
save(fits0, fits1, file = save_path)

fit <- RebuildPosterior(fits)
hat <- gelman(fit)
cat("mpsrf = ", hat$mpsrf, "\n")

options(digits = 2)
est_theta <- ggdmc::compare(fit, ps = p_vector)



p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * .5)
p0 <- ggdmc::plot(fit, start = 1)


# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = 1)
# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, start = 14)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, pll = F, start = fit@nmc * 0.5)

p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, start = fits[[1]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, start = fits[[2]]@nmc * 0.5)
p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, start = fits[[3]]@nmc * 0.5)

# p1 <- ggdmc::plot(fits[[1]], facet_chains = F, pll = F, start = 1)
# p1 <- ggdmc::plot(fits[[2]], facet_chains = F, pll = F, start = 1)
# p1 <- ggdmc::plot(fits[[3]], facet_chains = F, pll = F, start = 1)

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
