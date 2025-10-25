# q(save = "no")
cat("\n\n-------------------- Fitting hierarchical model2 --------------------")
rm(list = ls())
pkg <- c("ggdmc")
sapply(pkg, require, character.only = TRUE)
home_dir <- "/media/yslin/Tui/01_Projects/ggdmc/tests/testthat"
# home_dir <- "~/Documents/ggdmc/tests/testthat"
data_path <- file.path(home_dir, "Group9_gen_cdm/data/dina2.rda")
save_path <- file.path(home_dir, "Group11_hcdm/data/dina2.rda")


# long: columns = student (fct/int/chr), item (fct/int/chr), C (0/1), s (ignore)
long2wide <- function(df) {
    require(dplyr)
    require(tidyr)
    require(stringr)
    wide <- df %>%
        # keep only what we need
        select(student, item, C) %>%
        # make sure ids and items are clean, ordered, and item names look like E1..E28
        mutate(
            id   = as.integer(as.character(student)),
            item = paste0("E", as.integer(as.character(item))),
            C    = as.integer(C)
        ) %>%
        select(id, item, C) %>%
        # if there could be duplicates per (id,item), pick one; change to mean/max if needed
        distinct(id, item, .keep_all = TRUE) %>%
        # ensure a full rectangular layout (so missing item administrations show up as NA)
        complete(id, item) %>%
        # pivot to wide
        pivot_wider(names_from = item, values_from = C, values_fill = NA_integer_) %>%
        arrange(id)

    wide
}


cat("\nLoad data: ", data_path, "\n")
load(data_path)

cat("\nLoad fits data: ", data_path, "\n")
load(save_path)



# Fit random-effect models -------
# fits0 <- StartSampling(pop_dmis, pop_priors,
#     sub_migration_prob = 0.05, is_pblocked = TRUE,
#     thin = 2L, pop_debug = F, seed = 9032
# )

# save(fits0, file = save_path)

fits1 <- RestartSampling(fits0,
    sub_migration_prob = 0.00,
    is_pblocked = TRUE,
    thin = 1L, seed = 9032
)

save(fits0, fits1, file = save_path)


fits2 <- RestartSampling(fits1,
    pop_migration_prob = 0.02,
    sub_migration_prob = 0.02,
    is_pblocked = FALSE,
    is_hblocked = TRUE,
    thin = 1L, seed = 9032
)

save(fits0, fits1, fits2, file = save_path)


fits3 <- RestartSampling(fits2,
    pop_migration_prob = 0.02,
    sub_migration_prob = 0.02,
    is_pblocked = TRUE,
    is_hblocked = TRUE,
    thin = 1L, seed = 9032
)

save(fits0, fits1, fits2, fits3, file = save_path)

fits <- fits3
phi <- RebuildHyper(fits)
thetas <- RebuildPosteriors(fits)

hat <- gelman(phi)
cat("mpsrf = ", hat$mpsrf, "\n")
hats <- lapply(thetas, gelman)
sort(unlist(lapply(hats, function(x) x$mpsrf)), decreasing = TRUE)


options(digits = 2)
est_phi <- compare(phi, ps = true_vector)


# DT1 <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = fits[[1]]$phi@nmc * 0.85)
# DT2 <- ggdmc::prepare_thetas_data(fits[[2]]$subject_theta, start = fits[[2]]$phi@nmc * 0.85)
# DT3 <- ggdmc::prepare_thetas_data(fits[[3]]$subject_theta, start = fits[[3]]$phi@nmc * 0.85)
# p1 <- plot_thetas(DT1)
# p1 <- plot_thetas(DT2)
# p1 <- plot_thetas(DT3)
# p2 <- plot(fits[[1]]$phi, facet_chains = F, start = fits[[1]]$phi@nmc * 0.5)
# p2 <- plot(fits[[2]]$phi, facet_chains = F, start = fits[[2]]$phi@nmc * 0.5)
# p2 <- plot(fits[[3]]$phi, facet_chains = F, start = fits[[3]]$phi@nmc * 0.5)

# p3 <- plot(phi, facet_chains = F, start = phi@nmc * 0.5)
# p4 <- plot(phi, den = TRUE, pll = F, start = phi@nmc * 0.5)

# options(digits = 2)
# result <- compare_many(thetas, ps = ps)
# result <- compare_many(thetas, ps = ps, verbose = TRUE)


# Fixed-effect fits
#                guess1 guess2 guess3  pi1  pi2    pi3  slip1  slip2  slip3
# True           0.1000  0.200  0.300 0.10 0.20  0.300 0.0100 0.0200 0.0300
# 5% Estimate    0.0108  0.011  0.204 0.15 0.17  0.190 0.0048 0.0081 0.0057
# 50% Estimate   0.1068  0.116  0.253 0.20 0.23  0.251 0.0408 0.0577 0.0585
# 97.5% Estimate 0.2507  0.266  0.299 0.26 0.31  0.330 0.1230 0.1387 0.1628
# Median-True    0.0068 -0.084 -0.047 0.10 0.03 -0.049 0.0308 0.0377 0.0285

# $guess
#   partype parindex parameter  est    se   2.5 % 97.5 % item item.name
# 1   guess        1  E1_guess 0.13 0.071 -0.0128   0.26    1        E1
# 3   guess        3  E2_guess 0.14 0.069  0.0045   0.28    2        E2
# 5   guess        5  E3_guess 0.25 0.037  0.1753   0.32    3        E3
#   skillclass fixed free rule totindex
# 1          0 FALSE TRUE DINA        1
# 3          0 FALSE TRUE DINA        3
# 5          0 FALSE TRUE DINA        5
#
#   partype parindex parameter   est    se  2.5 % 97.5 % item item.name
# 2    slip        2   E1_slip 0.046 0.055 -0.062   0.15    1        E1
# 4    slip        4   E2_slip 0.064 0.050 -0.033   0.16    2        E2
# 6    slip        6   E3_slip 0.045 0.081 -0.114   0.20    3        E3
#   skillclass fixed free rule totindex
# 2          0 FALSE TRUE DINA        2
# 4          0 FALSE TRUE DINA        4
# 6          0 FALSE TRUE DINA        6
#
# $probs
#    partype parindex   parameter  est    se 2.5 % 97.5 % item item.name
# 7    probs        7 prob_class1 0.21 0.018  0.17   0.24    0
# 8    probs        8 prob_class2 0.23 0.018  0.20   0.27    0
# 9    probs        9 prob_class3 0.25 0.019  0.22   0.29    0
# 10   probs        0 prob_class4 0.31 0.027  0.26   0.36    0
#    skillclass fixed  free rule totindex
# 7           1 FALSE  TRUE             7
# 8           2 FALSE  TRUE             8
# 9           3 FALSE  TRUE             9
# 10          4 FALSE FALSE            10

# p0 <- ggdmc::plot(fit, pll = FALSE, den = TRUE, start = fit@nmc * 0.5)
# p1 <- ggdmc::plot(fit, facet_chains = F, start = fit@nmc * 0.5)
# p2 <- ggdmc::plot(fit, facet_chains = F, pll = F, start = fit@nmc * 0.5)


# long <- tibble::as_tibble(dat$responses)
# wide <- long2wide(long)
# wide
