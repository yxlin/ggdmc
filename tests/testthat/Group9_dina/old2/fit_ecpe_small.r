q(save = "no")
rm(list = ls())
pkg <- c("ggdmc", "ggdmcPrior", "ggdmcModel", "cdModel")
suppressPackageStartupMessages(pkg_ok <- sapply(pkg, require, character.only = TRUE))
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat/"
source(file.path(home_dir, "Group9_gen_cdm/00_helpers.r"))

# nitem <- 10
# dat <- CDM::data.ecpe$data[, paste0("E", seq_len(nitem))]
# Q <- CDM::data.ecpe$q.matrix[seq_len(nitem), ]

dat <- CDM::data.ecpe$data[, paste0("E", 8:10)]
Q <- CDM::data.ecpe$q.matrix[8:10, ]


ecpe <- CDM::din(data = dat, q.matrix = Q)
param <- CDM::IRT.se(ecpe, extended = TRUE)

p <- split(param, param$partype)
cdm_result <- tibble::tibble(
    guess = p$guess$est,
    slip = p$slip$est,
    g_se = p$guess$se,
    s_se = p$slip$se
)
cdm_result
p$probs



# Get the mapping between class numbers and skill profiles
profile_mapping <- data.frame(
    class_num = 1:nrow(ecpe$attribute.patt),
    profile = rownames(ecpe$attribute.patt),
    ecpe$attribute.patt.splitted, # The binary skill indicators
    class_prob = ecpe$attribute.patt$class.prob
)
print(profile_mapping)
