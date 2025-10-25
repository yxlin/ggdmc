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

#' Construct Data Table from the Theta Estiamtes
#'
#' Extracts and formats theta parameter samples or log-posterior values from
#' a model object for downstream analysis or plotting.
#'
#' @param x A model object containing MCMC samples (requires specific slots:
#'        @nmc, @nchain, @npar, @summed_log_prior, @log_likelihoods, @theta, @pnames)
#' @param start First iteration to include (default: 1)
#' @param end Last iteration to include (default: NULL uses all samples)
#' @param pll Logical indicating whether to prepare log-posterior values (TRUE)
#'        or theta parameters (FALSE) (default: TRUE)
#'
#' @return A data.table with columns:
#' \itemize{
#'   \item Iteration: The MCMC iteration number
#'   \item Chain: The chain identifier (factor)
#'   \item Parameter: The parameter name (factor)
#'   \item value: The sampled values or log-posterior values
#' }
#'
#' @details
#' This function prepares MCMC output in a long format suitable for:
#' \itemize{
#'   \item Trace plots and other diagnostics
#'   \item Posterior summary statistics
#'   \item Model comparison
#' }
#'
#' When pll=TRUE, it computes log-posterior = log_prior + log_likelihood.
#' When pll=FALSE, it extracts all theta parameters.
#'
#' @examples
#' \dontrun{
#' fits <- fits1
#' fit <- ggdmcDE::RebuildHyper(fits)
#' fit_thetas <- ggdmcDE:::RebuildPosteriors(fits)
#' DT <- prepare_theta_data(fit)
#' }
#' @export
#' @importFrom data.table rbindlist
prepare_theta_data <- function(x, start = 1L, end = NULL, pll = TRUE) {
    if (is.null(end)) {
        end <- x@nmc
    }
    if (end <= start) {
        stop("End must be greater than start")
    }

    iter <- start:end

    if (pll) {
        log_posterior <- x@summed_log_prior[, start:end] + x@log_likelihoods[, start:end]
        rownames(log_posterior) <- seq_len(x@nchain)

        v <- lapply(seq_len(x@nchain), function(k) {
            dd <- data.frame(
                Iteration = iter,
                Chain = k,
                Parameter = "log-posterior",
                value = sapply(log_posterior[k, ], c)
            )
        })

        DT <- data.table::rbindlist(v)
    } else {
        ## npar x nchain x nmc
        v <- lapply(seq_len(x@nchain), function(k) {
            dd <- data.frame(
                Iteration = rep(iter, x@npar),
                Chain = k,
                Parameter = rep(x@pnames, each = length(iter)),
                value = sapply(t(x@theta[, k, start:end]), c)
            )
        })

        DT <- data.table::rbindlist(v)
    }
    # DT
    DT$Parameter <- factor(DT$Parameter)
    DT$Chain <- factor(DT$Chain)
    return(DT)
}
