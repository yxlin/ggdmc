# q(save = "no")
cat("\n\n-------------------- 5 parameters --------------------")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))
cat("\nWorking directory: ", getwd(), "\n")

wkdir <- "~/Documents/ggdmc/tests/testthat/Group1/"
fn <- paste0(wkdir, "data/lba_data0.rda")
load(fn)

theta_input <- ggdmc::setThetaInput(nmc = 2, nchain = 3, pnames = model@pnames, thin = 1)
de_input <- ggdmc::setDEInput(
    sub_migration_prob = 0.00,
    nparameter = as.integer(sub_theta_input@nparameter),
    nchain = as.integer(sub_theta_input@nchain)
)
configs <- ggdmc::set_configs(prior = sub_priors, theta_input = sub_theta_input, de_input = de_input)
cfg <- configs[[1]]

pop_samples <- ggdmc::initialise_phi(pop_theta_input, pop_priors, pop_dmis, seed = 846671, verbose = FALSE)

# hyper_model <- ggdmcModel::BuildModel(
#     p_map = list(A = "1", B = "1", mean_v = "M", sd_v = "1", st0 = "1", t0 = "1"),
#     match_map = list(M = list(s1 = "r1", s2 = "r2")),
#     factors = list(S = c("s1", "s2")),
#     constants = c(sd_v = 1, st0 = 0),
#     accumulators = c("r1", "r2"),
#     type = "hyper"
# )
# hyper_dmi


# pop_samples$phi
fit0 <- run_hyper(cfg, hyper_dmi, pop_samples$phi)
head(run_hyper)
# fit0 <- run_subject(cfg, sub_dmis[[1]], samples, debug = T)
