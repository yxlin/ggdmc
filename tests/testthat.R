Sys.setenv("R_TESTS" = "")
## Workaround for the error,
## "cannot open file 'startup.Rs': No such file or directory" in Windows 10

library(testthat)
library(ggdmc)
library(ggplot2)
library(data.table)

cat("Start running testthat\n")
cat("in the directory:", getwd(), "\n")

cat("\n========================== Group 0 tests ==========================\n")

test_file(path = "testthat/Group0/test_BuildModel.R")
test_file(path = "testthat/Group0/test_rprior.R")
test_file(path = "testthat/Group0/test_GetParameterMatrix.R")
test_file(path = "testthat/Group0/test_likelihood.R")
test_file(path = "testthat/Group0/test_p_df.R")
test_file(path = "testthat/Group0/test_make_level_array.R")
test_file(path = "testthat/Group0/test_tnorm.R")
test_file(path = "testthat/Group0/test_BuildPrior.R")
test_file(path = "testthat/Group0/test_GetNsim.R")

cat("\n========================== Group 1 tests ==========================\n")

test_file(path = "testthat/Group1/test_checks.R")
test_file(path = "testthat/Group1/test_analysis_one.R")
test_file(path = "testthat/Group1/test_analysis_hyper.R")
test_file(path = "testthat/Group1/test_LBA1S.R")
test_file(path = "testthat/Group1/test_DDM1S.R")

cat("\n========================== Group 2 tests ==========================\n")

# test_file(path = "testthat/Group2/test_LBA8S.R")
# test_file(path = "testthat/Group2/test_DDM8S.R")
# test_file(path = "testthat/Group2/test_HLBA.R")
# test_file(path = "testthat/Group2/test_HDDM.R")

cat("\n====================== Testing plot functions ======================\n")

test_file(path = "testthat/plots/test_plotsubchain.R")
test_file(path = "testthat/plots/test_postpred.R")
test_file(path = "testthat/plots/test_plot_tnorm.R")