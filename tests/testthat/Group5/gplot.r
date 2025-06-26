gplot_thetas <- function(x, start = 1L, end = NULL, subchain = FALSE,
                         chains = NA, hide_legend = TRUE) {
    require(ggplot2)
    f1 <- ggplot(x, aes(x = Iteration, y = value, color = Chain)) +
        geom_line() +
        facet_wrap(~s, scales = "free") +
        theme_bw(base_size = 14) +
        ylab("Log-posterior likelihood")
    if (hide_legend) {
        f1 <- f1 + theme(legend.position = "none")
    }
    print(f1)
    return(invisible(f1))
}
