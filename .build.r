q(save = "no")
tools::showNonASCIIfile("DESCRIPTION")

usethis::use_gpl_license(version = 2) # or use_gpl_license(version = 2)
devtools::check(manual = TRUE, cran = TRUE) # Generates a PDF manual if needed
devtools::check()

devtools::test()
devtools::build()
devtools::install()
devtools::document()
devtools::install(quick = TRUE, upgrade = "never")

Sys.time()

rhub::rhub_setup()
rhub::rhub_platforms()

rhub::rhub_doctor()
windows <- rhub::rhub_check(
    gh_url = "https://github.com/yxlin/ggdmc",
    platforms = "windows"
)
macos <- rhub::rhub_check(
    gh_url = "https://github.com/yxlin/ggdmc",
    platforms = "macos"
)

clang_asan <- rhub::rhub_check(
    gh_url = "https://github.com/yxlin/ggdmc",
    platforms = "clang-asan"
)

valgrind <- rhub::rhub_check(
    gh_url = "https://github.com/yxlin/ggdmc",
    platforms = "valgrind"
)


macos_arm64 <- rhub::rhub_check(
    gh_url = "https://github.com/yxlin/ggdmc",
    platforms = "macos-arm64"
)


remove.packages("ggdmc") # Remove existing installation
unlink("src/*.o") # Remove object files
unlink("src/*.so") # Remove shared objects
unlink("src/RcppExports.*") # Remove generated Rcpp files
devtools::document()
devtools::install()
