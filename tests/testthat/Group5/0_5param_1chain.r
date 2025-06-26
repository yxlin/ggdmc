q(save = "no")
cat("\n\n-------------------- Old HB LBA model 0 --------------------")
rm(list = ls())
cat("\nWorking directory: ", getwd(), "\n")
pkg <- c("ggdmc", "ggdmcModel", "ggdmcPrior", "ggdmcLikelihood")
suppressPackageStartupMessages(tmp <- sapply(pkg, require, character.only = TRUE))

wkdir <- "~/Documents/ggdmc/tests/testthat/Group5/"
data_dir <- "~/Documents/ggdmc/tests/testthat/Group1/data/"

data_path <- paste0(data_dir, "lba_data0.rda")
helper_path <- paste0(wkdir, "gplot.r")
# save_path <- paste0(wkdir, "5param_1chain.rda")
# fg_path <- paste0(wkdir, "5param_1chain.pdf")

source(helper_path)
load(data_path)

theta_input <- ggdmcPhi::setThetaInput(nmc = 200L, nchain = 15L, pnames = pop_priors@pnames)

de_input <- ggdmcDE::setDEInput(
    nparameter = as.integer(theta_input@nparameter),
    nchain = as.integer(theta_input@nchain),
    is_hblocked = T, is_pblocked = F,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    pop_debug = T, sub_debug = F
)


configs <- ggdmcDE::set_configs(prior = pop_priors, theta_input = theta_input, de_input = de_input)

samples <- ggdmc::initialise_phi(theta_input, pop_priors, pop_dmis)
cfg <- configs[[1]]

fit0 <- ggdmcDE::run(cfg, pop_dmis, samples)
phi <- fit0$phi

plot(phi, start = phi@nmc * 0.5)
