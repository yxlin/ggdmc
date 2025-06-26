q(save = "no")
devtools::test()
devtools::build()
devtools::document()
devtools::install(quick = TRUE, upgrade = "never")
devtools::install()


remove.packages("ggdmc") # Remove existing installation
unlink("src/*.o") # Remove object files
unlink("src/*.so") # Remove shared objects
unlink("src/RcppExports.*") # Remove generated Rcpp files
devtools::document()
devtools::install()
