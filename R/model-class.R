################## ggdmcPhi class ---------------------------------------------
#' An S4 class to store MCMC sampling parameters.
#'
#' This class holds parameters that control the Markov Chain Monte Carlo (MCMC)
#' sampling process.
#'
#' @slot nmc A integer specifying the total number of MCMC iterations
#' (including burn-in).
#' @slot nchain A integer indicating the number of MCMC chains to run.
#' @slot thin A integer indicating the thinning interval.  This
#' determines how often samples are kept (e.g., a `thin` of 10 keeps every
#' 10th sample).  Thinning can help reduce autocorrelation in the MCMC samples.
#' @slot nparameter A integer storing the number of free parameter.
#' @slot pnames A string vector storing the name of the free parameter.
#' @slot report_length A integer specifying the interval at which
#' progress reports are printed during the MCMC run. When running parallel cores
#' on Windows, the print function may not work.
#' @slot max_init_attempts A integer specifying the maximum number of
#' attempts to initialise an MC sampler via the \code{rprior()} random number
#' generation.
#' @slot is_print A logical value specifying whether to print out the progress.
#' upon previous samples.
#'
#' @details
#' This class is used to encapsulate the key parameters that govern the
#' MCMC sampling process.  These parameters are essential for controlling
#' the length of the MCMC run, how often to keep samples, and the
#' storage and reporting of the sampling process.
#'
#' @section Objects from the Class:
#' An \code{theta_input} instance can be created using
#' the `setThetaInput` or the conventional `new("theta_input", ...)`
#' constructor, where `...` are named arguments corresponding to
#' the slots.
#' @export
setClass("theta_input",
    slots = c(
        nmc = "integer",
        nchain = "integer",
        thin = "integer",
        nparameter = "integer",
        pnames = "ANY",
        report_length = "integer",
        max_init_attempts = "integer",
        is_print = "logical"
    ),
    prototype = list(
        nmc = 200L,
        nchain = 15L,
        thin = 1L,
        nparameter = 5L,
        pnames = NULL,
        report_length = 100L,
        max_init_attempts = 10000L,
        is_print = TRUE
    )
)

#' Create a theta_input Object for MCMC Configuration
#'
#' Constructs a \code{theta_input} instance storing the parameter that controls
#' Markov Chain Monte Carlo (MCMC) sampling.
#'
#' @param nmc Integer specifying total number of MCMC iterations (including burn-in).
#'   Default: 500L
#' @param nchain Integer specifying number of MCMC chains facilitating the
#'   DE-MC sampling. If NULL (default), it will be set to 3 times the number of
#'   parameters.
#' @param thin Integer specifying thinning interval (keep every nth sample).
#'   Default: 1L (no thinning)
#' @param pnames Character vector naming the free parameters. Default: NULL
#' @param report_length Integer controlling how often progress is reported.
#'   Default: 100L (report every 100 iterations)
#' @param max_init_attempts Integer specifying maximum attempts to initialise
#'   valid starting values. Default: 1000L
#' @param is_print Logical controlling whether to print sampling progress.
#'   Default: TRUE
#'
#' @return A \code{theta_input} S4 object containing MCMC configuration
#'   parameters with the following slots:
#'   \itemize{
#'     \item \code{nmc}: Total iterations
#'     \item \code{nchain}: Number of chains
#'     \item \code{thin}: Thinning interval
#'     \item \code{nparameter}: Number of parameters (derived from \code{pnames})
#'     \item \code{pnames}: Parameter names
#'     \item \code{report_length}: Progress report frequency
#'     \item \code{max_init_attempts}: Initialisation attempts
#'     \item \code{is_print}: Print progress flag
#'   }
#'
#' @details
#' This constructor:
#' \itemize{
#'   \item Validates all input parameters
#'   \item Automatically sets \code{nchain} to 3 × number of parameters if NULL
#'   \item Derives \code{nparameter} from length of \code{pnames}
#'   \item Ensures all numeric parameters are converted to integers
#' }
#'
#' @examples
#' # A minimal LBA model with 5 free parameters
#' nmc <- 1000L
#' npar <- 5L
#' nchain <- npar * 3
#' thin <- 1L
#' pnames <- c("A", "B", "mean_v.false", "mean_v.true", "t0")
#'
#' sub_theta_input <- setThetaInput(
#'     nmc = nmc, nchain = as.integer(nchain), pnames = pnames,
#'     thin = thin
#' )
#'
#' # theta_input has its own print method
#' print(sub_theta_input)
#'
#' @importFrom methods new
#' @export
setThetaInput <- function(nmc = 500L, nchain = NULL, thin = 1L, pnames = NULL, report_length = 100L, max_init_attempts = 1000L, is_print = TRUE) {
    if (is.null(pnames)) {
        stop("Must provide a string vector, representing the ordered parameter name")
    }

    nparameter <- length(pnames)

    if (is.null(nchain)) {
        nchain <- 3 * nparameter
        message("Using ", nchain, " chains is to optimise the model")
    }


    new("theta_input",
        nmc = as.integer(nmc),
        nchain = as.integer(nchain),
        thin = as.integer(thin),
        nparameter = as.integer(nparameter),
        pnames = pnames,
        report_length = as.integer(report_length),
        max_init_attempts = as.integer(max_init_attempts),
        is_print = is_print
    )
}

#' Print Method for theta_input Objects
#'
#' Displays a human-readable summary of a \code{theta_input} object's contents.
#'
#' @param x A \code{theta_input} object to be printed
#'
#' @return Invisibly returns the input object \code{x}. Called for its side effect
#'   of printing to the console.
#'
#' @details
#' This method provides a concise console representation of \code{theta_input} objects,
#' showing all key configuration parameters in a standardised format:
#' \itemize{
#'   \item Basic MCMC dimensions (iterations, chains, thinning)
#'   \item Parameter information (count and names)
#'   \item Runtime control settings (reporting, initialization attempts)
#'   \item Operational flags (printing, accumulation)
#' }
#'
#' The output is designed to give users immediate visibility into the MCMC
#' configuration while avoiding overwhelming technical details.
#'
#' @examples
#' # Create and print a theta_input object
#' theta <- setThetaInput(
#'     nmc = 1000L,
#'     nchain = 3L,
#'     pnames = c("alpha", "beta", "sigma")
#' )
#' print(theta)
#'
#' # typing the object name to use the default print method
#' theta
#' @export
#' @name theta_input-print
#' @aliases print,theta_input-method
#' @docType methods
setMethod("print", "theta_input", function(x) {
    cat("ThetaInput Object Summary\n")
    cat("===================\n\n")

    # Basic configuration
    cat("Configuration Parameters:\n")
    cat("-----------------------\n")
    cat(sprintf("  MCMC Chains: %d\n", x@nchain))
    cat(sprintf("  Iterations per chain: %d\n", x@nmc))
    cat(sprintf("  Thinning interval: %d\n", x@thin))
    cat(sprintf("  Parameter numbers: %d\n", x@nparameter))
    cat(sprintf("  Report length: %d\n", x@report_length))
    cat(sprintf("  Max initialization attempts: %d\n", x@max_init_attempts))
    cat("\n")

    # Improved parameter names display
    cat("Model Parameters:\n")
    cat("----------------\n")
    if (length(x@pnames) == 0) {
        cat("  No parameter names specified.\n")
    } else if (length(x@pnames) < 10) {
        for (i in seq_along(x@pnames)) {
            cat(sprintf("  Parameter %d: %s\n", i, x@pnames[i]))
        }
    } else {
        message("Print the first 9 parameters")
        for (i in seq_len(9)) {
            cat(sprintf("  Parameter %d: %s\n", i, x@pnames[i]))
        }
    }
    cat("\n")

    # Display status flags
    cat("Status Flags:\n")
    cat("------------\n")
    cat(sprintf("  Print progress: %s\n", ifelse(x@is_print, "Yes", "No")))
    cat("\n")
    cat("Object Class: theta (S4)\n")

    invisible(x)
})


## Posterior ---------------------------------------------------------
#' An S4 class to represent an object storing posterior samples.
#'
#' @slot theta posterior samples.
#' @slot summed_log_prior summed log prior likelihoods.
#' @slot log_likelihoods logged likelihoods
#' @slot start the index of starting sample
#' @slot npar number of free parameters
#' @slot pnames free parameter names
#' @slot nmc number of Monte Carlo samples
#' @slot thin thinning length
#' @slot nchain number of Markov chains facilitating the DE-MCMC sampling
#' @exportClass posterior
#' @export
setClass("posterior",
    slot = c(
        theta = "ANY",
        summed_log_prior = "ANY",
        log_likelihoods = "ANY",
        start = "integer",
        npar = "integer",
        pnames = "character",
        nmc = "integer",
        thin = "integer",
        nchain = "integer"
    ),
    prototype = list(
        theta = NULL, summed_log_prior = NULL, log_likelihoods = NULL,
        start = 1L, npar = 5L, pnames = NULL, nmc = 200L, thin = 1L, nchain = 3L
    )
)

#' @importFrom stats ar residuals sd lm
spectrum0_ar <- function(x) {
    x <- as.matrix(x)
    v0 <- order <- numeric(ncol(x))
    names(v0) <- names(order) <- colnames(x)
    z <- 1:nrow(x)

    for (i in 1:ncol(x)) {
        lm.out <- stats::lm(x[, i] ~ z)
        if (identical(all.equal(stats::sd(stats::residuals(lm.out)), 0), TRUE)) {
            v0[i] <- 0
            order[i] <- 0
        } else {
            ar.out <- stats::ar(x[, i], aic = TRUE)
            v0[i] <- ar.out$var.pred / (1 - sum(ar.out$ar))^2
            order[i] <- ar.out$order
        }
    }
    return(list(spec = v0, order = order))
}

safespec0 <- function(x) {
    result <- try(spectrum0_ar(x)$spec)
    if (inherits(result, "try-error")) result <- NA
    result
}

#' Summarise Posterior Distribution
#'
#' Computes summary statistics for samples from a posterior distribution object.
#'
#' @param object A \code{posterior} class object containing MCMC samples
#' @param start First iteration to include in summary (default: 1)
#' @param end Last iteration to include in summary (default: NULL uses all samples)
#' @param probability Vector of probabilities for quantiles (default: c(0.05, 0.5, 0.975))
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{statistics}: Matrix of summary statistics (mean, sd, etc.)
#'   \item \code{quantiles}: Matrix of requested quantiles
#'   \item \code{start}: First iteration used
#'   \item \code{end}: Last iteration used
#' }
#'
#' @details
#' This method provides a statistical summary of the `theta` in the posterior `samples`:
#' \itemize{
#'   \item Calculates standard statistics (mean, median, SD) for all parameters
#'   \item Computes user-specified quantiles
#'   \item Handles burn-in through the \code{start} and \code{end} parameters
#' }
#'
#' The default quantiles (2.5%, 50%, 97.5%) provide the median and 95% credible interval.
#'
#' @examples
#' \dontrun{
#'
#' # Assuming you have set up DMI, priors etc.
#' # First make sure that the subject level has converged
#' fits0 <- StartSampling(pop_dmis, pop_priors,
#'     sub_migration_prob = 0.05,
#'     thin = 8L, pop_debug = F, seed = 9032
#' )
#'
#' # Then turn on the migration sampler at the population level
#' fits1 <- RestartSampling(fits0,
#'     pop_migration_prob = 0.02,
#'     sub_migration_prob = 0.00,
#'     thin = 4L, seed = 9032
#' )
#'
#' # Turn down the migration sampler, so that we may not fall
#' # into local minimal
#' fits2 <- RestartSampling(fits1,
#'     pop_migration_prob = 0.01,
#'     sub_migration_prob = 0.00,
#'     thin = 2L, seed = 9032
#' )
#'
#' fits <- fits2
#' phi <- RebuildHyper(fits)
#' thetas <- RebuildPosteriors(fits)
#' est_phi <- summary(phi)
#'
#' # str(est_phi)
#' # List of 2
#' #  $ statistics: num [1:10, 1:4] 1.024 0.246 0.15 2.283 0.406 ...
#' #   ..- attr(*, "dimnames")=List of 2
#' #   .. ..$ : chr [1:10] "loc_a" "loc_sz" "loc_t0" "loc_v" ...
#' #   .. ..$ : chr [1:4] "Mean" "SD" "Naive SE" "Time-series SE"
#' #  $ quantiles : num [1:10, 1:3] 0.9614 0.0266 0.134 0.2257 0.2275 ...
#' #   ..- attr(*, "dimnames")=List of 2
#' #   .. ..$ : chr [1:10] "loc_a" "loc_sz" "loc_t0" "loc_v" ...
#' #   .. ..$ : chr [1:3] "5%" "50%" "97.5%"
#' #  - attr(*, "start")= int 1
#' #  - attr(*, "end")= int 15000
#' #  - attr(*, "thin")= int 4
#' #  - attr(*, "nchain")= int 3
#' 
#' post_summary <- summary(phi, start = 501) # Discard first 500 as burn-in
#'
#' # Custom quantiles
#' detailed_summary <- summary(phi, probability = seq(0.1, 0.9, by = 0.1))
#' 
#' Extract specific elements
#' posterior_means <- post_summary$statistics[, "Mean"]
#' credible_intervals <- post_summary$quantiles[, c("5%", "97.5%")]
#' }
#'
#' @docType methods
#' @rdname summary.posterior-method
#' @aliases summary,posterior-method
#' @importFrom matrixStats rowMeans2 rowVars colMeans2 rowQuantiles
#' @export
setMethod(
    "summary", "posterior",
    function(object, start = 1L, end = NULL, probability = c(0.05, 0.5, 0.975)) {
        if (is.null(end)) {
            end <- object@nmc
        }
        if (length(probability) == 1) {
            stop("prob must be more than one element")
        }

        nchain <- object@nchain
        iter <- start:end
        niter <- length(iter)
        npar <- object@npar

        statnames <- c("Mean", "SD", "Naive SE", "Time-series SE")
        varstats <- matrix(
            nrow = npar, ncol = length(statnames),
            dimnames = list(object@pnames, statnames)
        )

        xtsvar <- matrix(nrow = nchain, ncol = object@npar)
        x0 <- NULL

        for (k in seq_len(nchain))
        {
            for (i in seq_len(npar))
            {
                xtsvar[k, i] <- safespec0(object@theta[i, k, ]) ## npar x nchain x nmc
            }
            x0 <- cbind(x0, object@theta[, k, ])
        }

        xmean <- matrixStats::rowMeans2(x0)
        xvar <- matrixStats::rowVars(x0)
        xtsvar2 <- matrixStats::colMeans2(xtsvar)
        varquant <- matrixStats::rowQuantiles(x0, probs = probability)

        varstats[, 1] <- xmean
        varstats[, 2] <- sqrt(xvar)
        varstats[, 3] <- sqrt(xvar / niter * nchain)
        varstats[, 4] <- sqrt(xtsvar2 / niter * nchain)
        varquant <- drop(varquant)
        varstats <- drop(varstats)
        rownames(varquant) <- object@pnames


        out <- list(statistics = varstats, quantiles = varquant)
        attr(out, "start") <- start
        attr(out, "end") <- end
        attr(out, "thin") <- object@thin
        attr(out, "nchain") <- nchain
        return(out)
    }
)

#' Compare True Parameters to Estimated Quantiles
#'
#' Compares true parameters to posterior quantiles from a fitted model object,
#' calculating estimation accuracy and bias.
#'
#' @param object A fitted instance from the a posterior class
#' @param start First iteration to include in comparison (default: 1)
#' @param end Last iteration to include in comparison (default: NULL uses all samples)
#' @param ps Named vector of true parameter values to compare against estimates
#' @param probability Vector of quantiles to compute (default: c(0.05, 0.5, 0.975))
#' @param verbose Logical indicating whether to print results (default: TRUE)
#'
#' @return A matrix with rows containing:
#' \itemize{
#'   \item True parameter values
#'   \item Estimated quantiles (labeled with percentage)
#'   \item Bias (Median - True)
#' }
#'
#' @details
#' This function:
#' \itemize{
#'   \item Extracts quantile summaries from the posterior samples
#'   \item Validates parameter names match between true values and estimates
#'   \item Computes estimation bias using the median
#'   \item Returns (and optionally prints) a comparison matrix
#' }
#'
#' @examples
#' \dontrun{
#' model <- ggdmcModel::BuildModel(
#'     p_map = list(A = "1", B = "1", mean_v = "M", sd_v = "1", st0 = "1", t0 = "1"),
#'     match_map = list(M = list(s1 = "r1", s2 = "r2")),
#'     factors = list(S = c("s1", "s2")),
#'     constants = c(sd_v = 1, st0 = 0),
#'     accumulators = c("r1", "r2"),
#'     type = "lba"
#' )
#'
#' fits0 <- ggdmcDE::StartSampling_subject(sub_dmis[[1]], sub_priors,
#'     sub_migration_prob = 0.00, thin = 1L, seed = 9032
#' )
#' fit0 <- ggdmcDE::RebuildPosterior(fits0)
#' options(digits = 2)
#' est_phi <- compare(fit0, ps = p_vector)
#' }
#'
#' @export
compare <- function(object, start = 1L, end = NULL, ps = NULL, probability = c(0.05, 0.5, 0.975), verbose = TRUE) {
    if (is.null(ps)) {
        stop("Must provide the true parameter")
    }
    # object <- phi
    # probability <- c(.05, .5, .975)
    # start <- 1
    # end <- phi@nmc * 0.5

    qs <- summary(object, start, end, probability)
    es <- qs$quantiles
    pnames_from_samples <- dimnames(qs)[[1]]
    pnames_from_true <- names(ps)

    if (!is.null(ps) && (!all(pnames_from_samples %in% pnames_from_true))) {
        stop("The name in the true vector (ps) do not match the parameter names stored in samples (ie object)")
    }
    # Create labels for the quantiles based on probability argument
    prob_labels <- paste0(100 * probability, "%")

    # Initialize output matrix
    out <- rbind(
        "True" = ps[order(names(ps))]
    )

    # Add each quantile estimate
    for (i in seq_along(probability)) {
        est <- qs$quantiles[names(ps), i]
        oest <- est[order(names(est))]
        out <- rbind(out, oest)
        rownames(out)[nrow(out)] <- paste0(prob_labels[i], " Estimate")
    }
    # Calculate bias using the median (50% quantile)
    median_idx <- which.min(abs(probability - 0.5)) # Find closest to median
    if (length(median_idx) == 0) median_idx <- ceiling(length(probability) / 2) # fallback
    est <- qs$quantiles[names(ps), median_idx]
    oest <- est[order(names(est))]
    bias <- oest - out["True", ]
    out <- rbind(out, "Median-True" = bias)

    if (verbose) {
        print(out)
    }

    return(out)
}


#' Summarise Multiple Posterior Objects
#'
#' @param objects A list of posterior objects to summarise
#' @param start First iteration to include in summary
#' @param end Last iteration to include in summary (NULL for all)
#' @param verbose Logical, whether to print summarised results. Default: FALSE
#'
#' @return Either a list of summary objects (if verbose=FALSE) or a matrix of 
#' means
#' @examples
#' \dontrun{
#' # Assuming you have set up DMI, priors etc.
#' # First make sure that the subject level has converged
#' fits0 <- StartSampling(pop_dmis, pop_priors,
#'     sub_migration_prob = 0.05,
#'     thin = 8L, pop_debug = F, seed = 9032
#' )
#'
#' # Then turn on the migration sampler at the population level
#' fits1 <- RestartSampling(fits0,
#'     pop_migration_prob = 0.02,
#'     sub_migration_prob = 0.00,
#'     thin = 4L, seed = 9032
#' )
#'
#' # Turn down the migration sampler, so that we may not fall
#' # into local minimal
#' fits2 <- RestartSampling(fits1,
#'     pop_migration_prob = 0.01,
#'     sub_migration_prob = 0.00,
#'     thin = 2L, seed = 9032
#' )
#'
#' fits <- fits2
#' phi <- RebuildHyper(fits)
#' thetas <- RebuildPosteriors(fits)
#' result <- summary_many(thetas)
#' 
#' # Return all details for every participant
#' result <- summary_many(thetas, verbose = TRUE)
#' }
#' @export
summary_many <- function(objects, start = 1L, end = NULL, 
                         verbose = FALSE) {
    # objects <- thetas
    #    start = 1
    #    end = NULL
    probability = c(0.05, 0.5, 0.975)
    #  verbose = TRUE


    # Validate input
    if (!is.list(objects)) {
        stop("objects must be a list of posterior objects")
    }

    # Generate summaries for all objects
    posterior_summaries <- lapply(objects, function(obj) {
        summary(object = obj, start = start, end = end, probability = probability)
    })

    # Return full summaries if not in verbose mode
    if (verbose) {
        return(posterior_summaries)
    }

    # Extract means from each summary for verbose output
    posterior_means <- lapply(posterior_summaries, function(summary) {
        summary$quantiles[, "50%"] # Extract mean estimates
    })

    # Combine means into matrix
    median_matrix <- do.call(rbind, posterior_means)

    # Add overall mean as last row
    combined_median <- rbind(
        median_matrix,
        "Mean" = colMeans(median_matrix)
    )

    # Set row names
    if (is.null(names(objects))) {
        rownames(combined_median) <- c(seq_along(objects), "Mean")
    } else {
        rownames(combined_median) <- c(names(objects), "Mean")
    }

    return(combined_median)
}


#' Compare Multiple Posterior Distributions Against True Values
#'
#' Computes and compares summary statistics across multiple posterior instances,
#' against user-provided true parameter values.
#'
#' @param objects A list of posterior objects from parameter optimisation
#' @param start First iteration to include in comparison (default: 1)
#' @param end Last iteration to include in comparison (default: NULL uses all samples)
#' @param ps Named matrix of true parameter values. Each row represents a participant
#' @param verbose Logical indicating whether to print results (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{summary_stats}: Array of summary statistics across all objects
#'   \item \code{quantiles}: Array of requested quantiles across all objects
#'   \item \code{bias}: Matrix of biases (median - true values) if ps is provided
#'   \item \code{comparison}: Summary comparison table if ps is provided
#' }
#'
#' @details
#' This function:
#' \itemize{
#'   \item Computes summary statistics across multiple posterior objects
#'   \item Calculates specified quantiles for each parameter in each object
#'   \item Optionally compares estimates to true values (when ps is provided)
#' }
#'
#' The function computes:
#' \itemize{
#'   \item Bias (difference between median estimate and true value)
#'   \item Comparison tables summarising accuracy
#' }
#'
#' @examples
#' \dontrun{
#' options(digits = 2)
#' est_thetas <- compare_many(fit_thetas, ps = ps)
#' #           A      B mean_v.false mean_v.true     t0
#' #  Mean 0.106  0.761       0.4326       2.524 0.2755
#' #  True 0.367  0.488       0.2402       2.502 0.3085
#' #  Diff 0.261 -0.274      -0.1924      -0.022 0.0330
#' #  Sd   0.039  0.075       0.1425       0.230 0.0491
#' #  True 0.106  0.093       0.1512       0.167 0.0556
#' #  Diff 0.068  0.018       0.0087      -0.063 0.0065
#' 
#' fits <- fits1
#' phi <- RebuildHyper(fits)
#' thetas <- RebuildPosteriors(fits)
#' result <- compare_many(thetas, ps = ps, verbose = TRUE)
#'               a     sz      t0       v        z
#' # Mean  1.00635  0.284  0.1463  2.6565  0.40525
#' # True  1.00235  0.251  0.1441  2.6892  0.38429
#' # Diff -0.00400 -0.033 -0.0022  0.0327 -0.02096
#' # Sd    0.05936  0.071  0.0155  0.5395  0.00851
#' # True  0.05884  0.011  0.0186  0.5343  0.00920
#' # Diff -0.00052 -0.060  0.0031 -0.0052  0.00069
#' }
#'
#' @export
#' @importFrom matrixStats colSds colMeans2
compare_many <- function(objects, start = 1L, end = NULL, ps = NULL,
                         verbose = TRUE) {
    if (is.null(ps)) {
        stop("Must provide the true parameter")
    }
    est <- summary_many(objects, start, end, FALSE)


    mean.est <- matrixStats::colMeans2(est)
    mean.ps <- matrixStats::colMeans2(ps)
    sd.est <- matrixStats::colSds(est)
    sd.ps <- matrixStats::colSds(ps)

    loc <- rbind(mean.est, mean.ps, mean.ps - mean.est)
    sca <- rbind(sd.est, sd.ps, sd.ps - sd.est)
    out <- rbind(loc, sca)

    rownames(out) <- c("Mean", "True", "Diff", "Sd", "True", "Diff")
    colnames(out) <- objects[[1]]@pnames


    if (verbose) {
        print(out)
    }
    return(list(summary = out, estimate = est))
}

#' Construct Data Table from the Theta Estiamtes
#'
#' Extracts and formats theta parameter samples or log-posterior values
#' from a model object for downstream analysis or plotting.
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


#' Plot Posterior Distributions or Traces using Lattice
#'
#' Visualise MCMC chains using lattice graphics, showing either trace plots of sampling
#' evolution or density plots of parameter distributions.
#'
#' @param x A \code{posterior} object containing MCMC samples
#' @param y Ignored (for consistency with generic plot method)
#' @param start First iteration to include (default: 1)
#' @param end Last iteration to include (default: NULL uses all samples)
#' @param pll Logical indicating whether to plot posterior log-likelihood (TRUE)
#'        or parameters (FALSE) (default: TRUE)
#' @param den Logical indicating whether to plot densities (TRUE) or traces (FALSE)
#'        when pll=FALSE (default: FALSE)
#' @param subchain Logical indicating whether to subset chains (default: FALSE)
#' @param chains Numeric vector specifying which chains to include when subchain=TRUE
#'        (default: NA randomly selects 3 chains)
#' @param facet_chains Logical indicating whether to facet chains vertically when
#'        pll=TRUE (default: TRUE)
#' @param hide_legend Logical indicating whether to hide legend (default: TRUE)
#' @param ... Additional arguments passed to lattice plotting functions
#'
#' @return Invisibly returns the lattice plot object. Primarily called for side effects.
#'
#' @details
#' This method provides three plotting modes using lattice graphics:
#' \enumerate{
#'   \item \strong{Log-likelihood traces} (pll=TRUE): Shows sampling evolution of
#'         log-likelihood values with optional chain faceting
#'   \item \strong{Parameter trace plots} (pll=FALSE, den=FALSE): Shows sampling
#'         evolution for each parameter with chains colored differently
#'   \item \strong{Parameter densities} (pll=FALSE, den=TRUE): Shows kernel density
#'         estimates for each parameter's posterior distribution
#' }
#'
#' The function automatically handles:
#' \itemize{
#'   \item Chain coloring using rainbow palette
#'   \item Free y-axis scaling across facets
#'   \item Intelligent legend placement
#'   \item Subsampling of chains when requested
#' }
#'
#' @examples
#' \dontrun{
#' # Plot posterior log-likelihood traces
#' plot(posterior_obj, pll = TRUE)
#'
#' # Plot parameter traces (first 500 iterations as burn-in)
#' plot(posterior_obj, pll = FALSE, start = 501)
#'
#' # Plot parameter densities for 3 random chains
#' plot(posterior_obj, pll = FALSE, den = TRUE, subchain = TRUE)
#' }
#'
#' @export
#' @docType methods
#' @rdname plot.posterior-method
#' @aliases plot,posterior-method
#'
#' @importFrom lattice xyplot densityplot lattice.options
#' @importFrom grDevices rainbow
setMethod(
    "plot",
    signature(x = "posterior"),
    function(x, y, start = 1L, end = NULL, pll = TRUE, den = FALSE,
             subchain = FALSE, chains = NA, facet_chains = FALSE,
             hide_legend = TRUE, ...) {
        if (!requireNamespace("lattice", quietly = TRUE)) {
            stop("Please install lattice: install.packages('lattice')")
        }

        # Prepare data
        message("Construct a data table instance for plotting...")
        DT <- prepare_theta_data(x, start, end, pll)

        # Handle subchain selection
        if (subchain) {
            if (any(is.na(chains))) {
                nsubchain <- 3
                message("Randomly selecting ", nsubchain, " chains")
                chains <- sample(unique(DT$Chain), nsubchain)
            }
            cat("Plotting chains: ", chains, "\n")
            DT <- DT[DT$Chain %in% chains, ]
        }

        # Convert Chain to factor for proper coloring
        DT$Chain <- factor(DT$Chain)

        # Define common lattice parameters
        lattice::lattice.options(
            default.theme = list(
                plot.symbol = list(col = "black"),
                superpose.line = list(col = rainbow(length(unique(DT$Chain)))),
                strip.background = list(col = "lightgray")
            )
        )

        if (pll) {
            # Posterior log-likelihood plot
            if (facet_chains) {
                p <- lattice::xyplot(
                    value ~ Iteration | Chain,
                    data = DT,
                    type = "l",
                    layout = c(1, length(unique(DT$Chain))),
                    xlab = "Iteration",
                    ylab = "Posterior log-likelihood",
                    scales = list(y = list(relation = "free")),
                    as.table = TRUE,
                    par.settings = list(
                        superpose.line = list(col = rainbow(length(unique(DT$Chain))))
                    )
                )
            } else {
                p <- lattice::xyplot(
                    value ~ Iteration,
                    groups = DT$Chain,
                    data = DT,
                    type = "l",
                    xlab = "Iteration",
                    ylab = "Posterior log-likelihood",
                    auto.key = if (!hide_legend) list(space = "right") else FALSE,
                    par.settings = list(
                        superpose.line = list(col = rainbow(length(unique(DT$Chain))))
                    )
                )
            }
        } else if (den) {
            # Density plot
            p <- lattice::densityplot(
                ~ value | Parameter,
                groups = DT$Chain,
                data = DT,
                xlab = "Parameter value",
                ylab = "Density",
                plot.points = FALSE,
                auto.key = if (!hide_legend) list(space = "right") else FALSE,
                scales = list(relation = "free"),
                par.settings = list(
                    superpose.line = list(col = rainbow(length(unique(DT$Chain))))
                )
            )
        } else {
            # Trace plot
            p <- lattice::xyplot(
                value ~ Iteration | Parameter,
                groups = DT$Chain,
                data = DT,
                type = "l",
                xlab = "Iteration",
                ylab = "",
                scales = list(y = list(relation = "free")),
                auto.key = if (!hide_legend) list(space = "right") else FALSE,
                par.settings = list(
                    superpose.line = list(col = rainbow(length(unique(DT$Chain))))
                )
            )
        }

        print(p)
        return(invisible(p))
    }
)

#' Prepare Theta Parameters Data from MCMC Samples
#'
#' Extracts and formats theta parameter samples from MCMC output for analysis or visualization,
#' with options for subchain extraction and chain selection.
#'
#' @param x A model object containing MCMC samples (requires slots: @theta, @pnames,
#'        @nchain, @nmc, and optionally @subchain_theta if subchain=TRUE)
#' @param start First iteration to include (default: 1)
#' @param end Last iteration to include (default: NULL uses all available samples)
#' @param subchain Logical indicating whether to extract subchain samples (default: FALSE)
#' @param chains Numeric vector specifying which chains to include (default: NA includes all chains)
#'
#' @return A data.table in long format with columns:
#' \itemize{
#'   \item Iteration: The MCMC iteration number
#'   \item Chain: The chain identifier (factor)
#'   \item Parameter: The parameter name (factor)
#'   \item value: The sampled parameter values
#' }
#'
#' @details
#' This function provides flexible extraction of theta parameters:
#' \itemize{
#'   \item Handles both main chains and subchains (when available)
#'   \item Allows selection of specific chains
#'   \item Supports iteration range selection
#'   \item Returns data in tidy format for easy plotting with ggplot2
#' }
#'
#' When subchain=TRUE, the function will attempt to access @subchain_theta slot.
#' If chains are specified, only those chains will be included in the output.
#'
#' @examples
#' \dontrun{
#' # Basic usage with main chains
#' theta_data <- prepare_thetas_data(fit, start = 500)
#'
#' # Extract subchains instead of main chains
#' subchain_data <- prepare_thetas_data(fit, subchain = TRUE)
#'
#' # Specific chains only
#' chain1_data <- prepare_thetas_data(fit, chains = 1)
#'
#' # Using with ggplot2
#' library(ggplot2)
#' ggplot(
#'     prepare_thetas_data(fit),
#'     aes(x = Iteration, y = value, color = Chain)
#' ) +
#'     geom_line(alpha = 0.7) +
#'     facet_wrap(~Parameter, scales = "free_y") +
#'     labs(title = "Theta Parameter Traces")
#' }
#'
#' @export
#' @importFrom data.table rbindlist
prepare_thetas_data <- function(x, start = 1L, end = NULL, subchain = FALSE, chains = NA) {
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("data.table package required but not installed")
    }

    if (!is.list(x)) stop("x must be a list")
    if (length(x) == 0) stop("x cannot be empty")

    nsubject <- length(x)

    # Parallel processing for large datasets
    if (nsubject > 50 && requireNamespace("parallel", quietly = TRUE)) {
        DT_list <- parallel::mclapply(x, prepare_theta_data,
            start = start, end = end, pll = TRUE
        )
    } else {
        DT_list <- lapply(x, prepare_theta_data, start = start, end = end, pll = TRUE)
    }

    if (subchain) {
        if (any(is.na(chains))) {
            nsubchain <- 3
            message("Randomly select ", nsubchain, "chains")
            # Get chains from first subject as sample
            chains <- sample(unique(DT_list[[1]]$Chain), nsubchain)
        }
        message("Plot chains: ", paste(chains, collapse = ", "))
        DT_list <- lapply(DT_list, function(dt) {
            dt[dt$Chain %in% chains, ]
            # dt[Chain %in% chains]
        })
    }

    # More efficient subject naming
    subjectnames <- names(x) %||% seq_len(nsubject)

    # Efficient rbind with proper data.table environment
    x0 <- lapply(seq_len(nsubject), function(i) {
        dt <- data.table::copy(DT_list[[i]]) # Avoid modifying by reference
        dt$s <- factor(subjectnames[i])
        dt
    })


    DT <- data.table::rbindlist(x0)

    return(DT)
}

#' Plot Theta Parameters using Lattice Graphics
#'
#' Visualises theta parameter traces using lattice graphics with faceting by parameter.
#' Provides flexible subject selection options.
#'
#' @param x A data frame/table containing columns: Iteration, value, Chain, and 
#' s (subject label)
#' @param start First iteration to include (default: 1)
#' @param end Last iteration to include (default: NULL uses all samples)
#' @param subchain Logical indicating whether to subset chains (default: FALSE)
#' @param chains Numeric vector specifying which chains to include when subchain=TRUE
#'        (default: NA includes all chains)
#' @param hide_legend Logical indicating whether to hide legend (default: TRUE)
#' @param subjects Number of subjects to plot or vector of specific subjects:
#'        - NULL (default): plots all subjects, if less than 5. Plot the first 5, when
#'                          the subject number exceeeds 5.
#'        - Integer N: randomly samples N subjects
#'        - Character vector: plots specified subjects
#' @param seed Optional random seed for reproducible sampling (default: NULL)
#'
#' @return Invisibly returns the lattice plot object. Primarily called for side effects.
#'
#' @examples
#' \dontrun{
#' # Prepare data
#' fits <- fits1
#' phi <- RebuildHyper(fits)
#' thetas <- RebuildPosteriors(fits)
#'
#' DT <- ggdmc::prepare_thetas_data(fits[[1]]$subject_theta, start = 250)
#' DT <- ggdmc::prepare_thetas_data(thetas, start = 5000)
#' p1 <- plot_thetas(DT)
#' p1 <- plot_thetas(DT, start = 300, end = 400)
#' p1 <- plot_thetas(DT, start = 300, end = 400, subjects = 5)
#' p1 <- plot_thetas(DT, start = 300, end = 400, subjects = as.character(1:10))
#' p1 <- plot_thetas(DT, start = 300, end = 400, max_subjects = 8)
#' }
#'
#' @export
#' @importFrom lattice xyplot
#' @importFrom grDevices rainbow
plot_thetas <- function(x, start = 1L, end = NULL, subchain = FALSE,
                       chains = NA, hide_legend = TRUE, subjects = NULL,
                       seed = NULL, max_subjects = 50) {
    if (!requireNamespace("lattice", quietly = TRUE)) {
        stop("Please install lattice: install.packages('lattice')")
    }

    # Filter data if needed
    if (!is.null(end)) {
        x <- x[x$Iteration >= start & x$Iteration <= end, ]
    } else {
        x <- x[x$Iteration >= start, ]
    }

    # Handle subchain selection
    if (subchain) {
        if (any(is.na(chains))) {
            chains <- unique(x$Chain)
        }
        x <- x[x$Chain %in% chains, ]
    }

    # Handle subject selection
    all_subjects <- sort(unique(x$s))
    n_total <- length(all_subjects)

    if (!is.null(subjects)) {
        if (is.numeric(subjects)) {
            # Random sampling case
            if (!is.null(seed)) set.seed(seed)
            selected <- sample(all_subjects, min(subjects, n_total))
            message("Plotting ", length(selected), " randomly selected subjects")
            x <- x[x$s %in% selected, ]
        } else if (is.character(subjects)) {
            # Specific subjects case
            missing <- setdiff(subjects, all_subjects)
            if (length(missing) > 0) {
                warning("Some subjects not found: ", paste(missing, collapse = ", "))
            }
            selected <- intersect(subjects, all_subjects)
            if (length(selected) == 0) stop("No valid subjects specified")
            message("Plotting ", length(selected), " specified subjects")
            x <- x[x$s %in% selected, ]
        }
    } else if (!is.null(max_subjects) && n_total > max_subjects) {
        # Plot first max_subjects if specified and there are many subjects
        selected <- all_subjects[1:max_subjects]
        message(
            n_total, " subjects detected. Plotting first ", max_subjects, ". ",
            "Use 'subjects' parameter to customise or set max_subjects = NULL to plot all."
        )
        x <- x[x$s %in% selected, ]
    } else {
        # Default case: plot all subjects
        message("Plotting all ", n_total, " subjects")
    }

    # Determine number of subjects to plot
    n_subjects <- length(unique(x$s))
    
    # Calculate layout dimensions for grid arrangement
    if (n_subjects > 1) {
        ncol <- min(ceiling(sqrt(n_subjects)), 4)  # Max 4 columns for readability
        nrow <- ceiling(n_subjects / ncol)
        layout_dim <- c(ncol, nrow)
    } else {
        layout_dim <- c(1, 1)  # Single plot layout
    }

    # Convert Chain to factor for proper coloring
    x$Chain <- factor(x$Chain)

    # Create plot
    p <- lattice::xyplot(
        value ~ Iteration | s,
        groups = x$Chain,
        data = x,
        type = "l",
        xlab = "Iteration",
        ylab = "Log-posterior likelihood",
        scales = list(y = list(relation = "free")),
        auto.key = if (!hide_legend) {
            list(space = "right", title = "Chain", cex.title = 1)
        } else {
            FALSE
        },
        par.settings = list(
            superpose.line = list(col = rainbow(length(unique(x$Chain)))),
            strip.background = list(col = "lightgray")
        ),
        layout = layout_dim,  # Use the calculated layout dimensions
        as.table = TRUE
    )

    print(p)
    invisible(p)
}


################## ggdmcDE class ------------------------

#' Differential Evolution (DE) MCMC Input Configuration
#'
#' Provides configuration objects and methods for Differential Evolution
#' Markov Chain Monte Carlo (DE-MCMC) sampling, including migration
#' probabilities and tuning parameters.
#'
#' @slot pop_migration_prob Numeric. Population-level chain migration
#' probability (default = 0.00)
#' @slot sub_migration_prob Numeric. Subject-level chain migration
#' probability (default = 0.00)
#' @slot gamma_precursor Numeric. Scaling factor for proposal
#' distribution (default = 2.38)
#' @slot rp Numeric. Random perturbation factor for proposals (default = 0.001)
#' @slot is_hblocked Logical. Whether hyperparameters use block
#' updating (default = FALSE)
#' @slot is_pblocked Logical. Whether subject parameters use block
#' updating (default = FALSE)
#' @slot nparameter Integer. Number of parameters in model (default = 5L)
#' @slot nchain Integer. Number of MCMC chains (default = 3L)
#' @slot pop_debug Logical. Whether to print hyper-level debugging information
#' @slot sub_debug Logical. Whether to print subject-level debugging information
#'
#' @section Functions:
#' \describe{
#'   \item{`setClass("de_input")`}{Defines the DE-MCMC configuration class}
#'   \item{`setDEInput()`}{Constructor function for creating `de_input` objects}
#'   \item{`print` method}{Displays a readable summary of DE-MCMC configuration}
#' }
#'
#' @param pop_migration_prob Population-level migration probability (0-1)
#' @param sub_migration_prob Subject-level migration probability (0-1)
#' @param gamma_precursor Scaling factor for proposal generation (typically 2.38)
#' @param rp Random perturbation size for proposals
#' @param is_hblocked Use block updating for hyperparameters?
#' @param is_pblocked Use block updating for subject parameters?
#' @param nparameter Number of parameters in the model
#' @param nchain Number of parallel chains
#' @param pop_debug Logical. If \code{TRUE}, enables verbose diagnostic
#' output at the population level for debugging purposes (default = FALSE).
#' @param sub_debug Logical. If \code{TRUE}, enables verbose diagnostic
#' output at the subject level for debugging purposes (default = FALSE).
#' @param x A `de_input` object (for print method)
#'
#' @return
#' \describe{
#'   \item{`setDEInput()`}{Returns a `de_input` class object}
#'   \item{`print` method}{Invisibly returns the input object}
#' }
#'
#' @details
#' The `de_input` class encapsulates key configuration parameters for DE-MCMC sampling:
#' \itemize{
#'   \item Chain migration probabilities control between-chain mixing
#'   \item Gamma precursor affects proposal jump sizes
#'   \item Random perturbation (rp) maintains chain diversity
#'   \item Block updating options can improve efficiency for correlated parameters
#' }
#'
#' @examples
#' # Create a default configuration
#' de_config <- setDEInput()
#'
#' # Custom configuration with higher migration and block updating
#' de_custom <- setDEInput(
#'     pop_migration_prob = 0.1,
#'     sub_migration_prob = 0.05,
#'     is_hblocked = TRUE,
#'     nparameter = 7L
#' )
#'
#' # Print configuration
#' print(de_custom)
#'
#' @seealso
#' \code{\link{StartSampling}} for functions using these configurations,
#' \code{\link{setClass}} for S4 class definitions
#'
#' @rdname de_input
#' @aliases de_input-class setDEInput print,de_input-method
#' @export
setClass("de_input",
    slots = c(
        pop_migration_prob = "numeric",
        sub_migration_prob = "numeric",
        gamma_precursor = "numeric",
        rp = "numeric",
        is_hblocked = "logical",
        is_pblocked = "logical",
        nparameter = "integer",
        nchain = "integer",
        pop_debug = "logical",
        sub_debug = "logical"
    ),
    prototype = list(
        pop_migration_prob = 0.00,
        sub_migration_prob = 0.00,
        gamma_precursor = 2.38,
        rp = 0.001,
        is_hblocked = FALSE,
        is_pblocked = FALSE,
        nparameter = 5L,
        nchain = 3L,
        pop_debug = FALSE,
        sub_debug = FALSE
    )
)


#' @rdname de_input
#' @export
setDEInput <- function(
    pop_migration_prob = 0.00, sub_migration_prob = 0.00,
    gamma_precursor = 2.38, rp = 0.001, is_hblocked = FALSE, is_pblocked = FALSE,
    nparameter = 5L, nchain = 3L, pop_debug = FALSE, sub_debug = FALSE) {
    # Validate inputs
    if (any(pop_migration_prob < 0 | sub_migration_prob > 1)) {
        stop("migration_prob must be between 0 and 1.")
    }

    if (gamma_precursor <= 0) {
        stop("gamma_precursor must be a positive number.")
    }
    if (rp <= 0) {
        stop("rp must be a positive number.")
    }
    if (!is.logical(is_hblocked) | !is.logical(is_pblocked)) {
        stop("is_hblocked / is_pblocked must be a logical value (TRUE or FALSE).")
    }

    if (!is.integer(nparameter)) {
        message("You enter a numeric nparameter, but we need it to be an integer.")
        nparameter <- as.integer(nparameter)
    }

    if (!is.integer(nchain)) {
        message("You enter a numeric nchain, but we need it to be an integer.")
        nchain <- as.integer(nchain)
    }

    # Create and return the object
    new("de_input",
        pop_migration_prob = pop_migration_prob,
        sub_migration_prob = sub_migration_prob,
        gamma_precursor = gamma_precursor,
        rp = rp,
        is_hblocked = is_hblocked,
        is_pblocked = is_pblocked,
        nparameter = nparameter,
        nchain = nchain,
        pop_debug = pop_debug,
        sub_debug = sub_debug
    )
}

#' @rdname de_input
#' @export
setMethod(
    "print", "de_input",
    function(x) {
        cat('An object of class "de_input"\n')
        cat("-----------------------------\n")
        cat(sprintf(
            "Migration probabilities (pop & sub): %.2f, %.2f\n", x@pop_migration_prob,
            x@sub_migration_prob
        ))
        cat(sprintf("Tuning parameters (gamma precursor & rp): %.2f, %.3f\n", x@gamma_precursor, x@rp))
        cat(sprintf("Using population block update? %s\n", ifelse(x@is_hblocked, "Yes", "No")))
        cat(sprintf("Using subject block update? %s\n", ifelse(x@is_pblocked, "Yes", "No")))

        cat(sprintf("population level debug? %s\n", ifelse(x@pop_debug, "Yes", "No")))
        cat(sprintf("subject level debug? %s\n", ifelse(x@sub_debug, "Yes", "No")))
        cat(sprintf("The number of parameters and chains: %d, %d\n", x@nparameter, x@nchain))
    }
)

#' MCMC Configuration Class and Constructor
#'
#' The `config` class collects complete MCMC sampling configurations, including
#' priors, parameter specifications, and Differential Evolution (DE) tuning
#' parameters. The `set_configs()` function provides a convenient constructor
#' for creating configuration objects.
#'
#' @slot prior An object of class `prior` from ggdmcPrior package containing
#' prior distributions
#' @slot theta_input An object of class `theta_input` containing parameter
#' specifications
#' @slot de_input An object of class `de_input` containing DE-MCMC tuning
#' parameters
#' @slot seed An integer vector of random seeds for worker chains
#' @slot main_seed An integer for the main random seed
#' @slot core_id An integer vector mapping cores to seeds
#'
#' @param prior A `prior` object specifying the model's prior distributions
#' @param theta_input A `theta_input` object specifying parameter configurations
#' @param de_input A `de_input` object specifying DE-MCMC tuning parameters
#' @param ncore Integer specifying number of cores to use (default = 3L)
#' @param seed Optional random seed (if NULL, seeds are generated automatically)
#'
#' @return
#' \describe{
#'   \item{`config` class}{An S4 object containing complete MCMC configuration}
#'   \item{`set_configs()`}{Returns a validated `config` object}
#' }
#'
#' @details
#' The `config` class integrates all components needed for MCMC sampling:
#' \itemize{
#'   \item **Prior distributions** (from ggdmcPrior)
#'   \item **Parameter specifications**
#'   \item **DE-MCMC tuning parameters** (jump sizes: rp and gamma_precursor; 
#' migration probabilities, etc.)
#'   \item **Random number generation** with proper seed management
#' }
#'
#' The constructor `set_configs()`:
#' \itemize{
#'   \item Automatically generates seeds when not specified
#'   \item Ensures worker seeds differ from main seed
#'   \item Validates core assignments
#'   \item Checks consistency between components
#' }
#' 
#' Each core runs a true independent `chain`, differing from the 
#' `chains`/`chromosome` used to help the DE-MC sampling to work.
#'
#' @section Validation Rules:
#' The class includes these validity checks:
#' \itemize{
#'   \item Worker seeds must differ from main seed
#'   \item Core IDs must be positive integers
#'   \item Components must be properly classed objects
#' }
#'
#' @examples
#' \dontrun{
#' # To create a configuration profile, we first set up a de_input
#' 
#' de_input <- ggdmc::setDEInput(
#'     sub_migration_prob = 0.00,
#'     nparameter = as.integer(sub_theta_input@nparameter),
#'     nchain = as.integer(sub_theta_input@nchain)
#' )
#' 
#' # Then a theta_input
#' theta_input <- ggdmc::setThetaInput(nmc = 2, nchain = 3, pnames = model@pnames, thin = 1)
#' 
#' # Finally we can use set_configs to create  configuration profile
#' configs <- ggdmc::set_configs(prior = sub_priors, theta_input = theta_input, 
#' de_input = de_input)
#' 
#' # The default setting is to set up three configurations for three independent chains
#' # You have to ensure that you have three different starting `samples`.
#' cfg <- configs[[1]]
#' }
#'
#' @rdname mcmc_config
#' @aliases config-class set_configs
#' @importClassesFrom ggdmcPrior prior
#' @export
setClass("config",
    slots = c(
        prior = "prior",
        theta_input = "theta_input",
        de_input = "de_input",
        seed = "integer",
        main_seed = "integer",
        core_id = "integer"
    ),
    validity = function(object) {
        if (object@seed == object@main_seed) {
            return("Worker seed should differ from main seed")
        }
        if (object@core_id < 1) {
            return("core_id must be positive")
        }
        TRUE
    }
)


#' @rdname mcmc_config
#' @export
set_configs <- function(prior = NULL, theta_input = NULL, de_input = NULL, ncore = 3L, seed = NULL) {
    if (is.null(prior)) {
        stop("Must supply a joint prior distribution.")
    }
    if (is.null(theta_input)) {
        stop("Must supply a theta_input.")
    }
    if (is.null(de_input)) {
        stop("Must supply a de_input.")
    }
    if (is.null(prior@h_prior)) {
        if (prior@nparameter != theta_input@nparameter) {
            stop("The number of parameters in `prior` are not the same as that in the theta_input ")
        }
    }


    # Generate seeds based on number of cores
    if (is.null(seed)) {
        main_seed <- as.integer(Sys.time()) %% 1000000L
    } else {
        main_seed <- seed
    }

    set.seed(main_seed)
    seeds <- sample.int(1e6, ncore)

    # Print header information
    message("Sampling settings:")
    message("---------------------")
    message("Main seed: ", main_seed)
    message("Number of cores: ", ncore)
    message("\nSeeds for each instance:")

    # Print seeds in a more readable format
    for (i in seq_len(ncore)) {
        message("  Instance ", i, ": ", seeds[i])
    }

    # Create list of configs, each with a different seed
    config_list <- lapply(seq_len(ncore), function(i) {
        new("config",
            prior = prior,
            theta_input = theta_input,
            de_input = de_input,
            seed = as.integer(seeds[i]), # Each config gets its own seed
            main_seed = as.integer(main_seed),
            core_id = as.integer(i)
        )
    })

    return(config_list)
}


########### gelman method ------------------------------------
#' @title Brook-Gelman Diagnostic (R-hat)
#' @description Calculates potential scale reduction factors for MCMC chains.
#' @name gelman
#' @aliases gelman-methods
NULL

## Helper Functions -----------------------------------------------------------
#' @title Calculate within-chain statistics
#' @description Internal function to compute within-chain variances (W)
#' @param theta Array of MCMC samples (iterations × chains × parameters)
#' @return List with W matrix and S2 array
#' @keywords internal
#' @importFrom stats var
.calculate_within_chain_stats <- function(theta) {
    nchain <- dim(theta)[2]
    npar <- dim(theta)[1]

    # Calculate variance for each chain
    chain_variances <- sapply(1:nchain, function(k) var(t(theta[, k, ])))

    list(
        W = matrix(rowMeans(chain_variances), nrow = npar), # Mean within-chain variance
        S2 = array(chain_variances, dim = c(npar, npar, nchain)) # Chain-specific variances
    )
}

#' @title Calculate between-chain statistics
#' @description Internal function to compute between-chain variances
#' @param xbar Matrix of chain means (parameters × chains)
#' @param niter Number of iterations used
#' @return Between-chain covariance matrix
#' @keywords internal
#' @importFrom stats var
.calculate_between_chain_stats <- function(xbar, niter) {
    niter * var(t(xbar))
}

#' @title Compute PSRF components
#' @description Internal function for R-hat calculations
#' @inheritParams gelman
#' @param W Within-chain covariance matrix
#' @param B Between-chain covariance matrix
#' @param niter Number of iterations
#' @param nchain Number of chains
#' @return List with R-hat components
#' @keywords internal
#' @importFrom stats qf var
.compute_psrf_components <- function(W, B, niter, nchain, conf) {
    w <- diag(W)
    b <- diag(B)

    # Variance components
    var_w <- apply(W, 1, var) / nchain
    var_b <- (2 * b^2) / (nchain - 1)

    # Degrees of freedom adjustments
    V <- (niter - 1) / niter * w + (1 + 1 / nchain) / niter * b
    var_V <- ((niter - 1)^2 * var_w + (1 + 1 / nchain)^2 * var_b) / niter^2
    df_V <- (2 * V^2) / var_V
    df_adj <- (df_V + 3) / (df_V + 1)

    # R-hat calculations
    R2_fixed <- (niter - 1) / niter
    R2_random <- (1 + 1 / nchain) * (1 / niter) * (b / w)

    list(
        R2_estimate = R2_fixed + R2_random,
        R2_upper = R2_fixed + qf(0.5 * (1 + conf), nchain - 1, (2 * w^2) / var_w) * R2_random,
        df_adj = df_adj
    )
}

#' @title Get multivariate PSRF
#' @description Internal function for multivariate R-hat
#' @inheritParams gelman
#' @param npar Number of parameters
#' @param W Within-chain covariance matrix
#' @param B Between-chain covariance matrix
#' @return Multivariate PSRF or NULL
#' @keywords internal
.get_multivariate_psrf <- function(npar, niter, W, B, multivariate) {
    if (!multivariate || npar < 2) {
        return(NULL)
    }

    CW <- chol(W)
    emax <- eigen(
        backsolve(CW, t(backsolve(CW, B, transpose = TRUE)),
            transpose = TRUE
        ),
        symmetric = TRUE,
        only.values = TRUE
    )$values[1]

    mpsrf <- sqrt((1 - 1 / niter) + (1 + 1 / npar) * emax / niter)
    return(mpsrf)
}


## Internal Hyper Function -----------------------------------------------------


#' @title Brook-Gelman diagnosis
#' @description Internal function for hierarchical models
#' @inheritParams gelman
#' @return List with psrf and mpsrf
#' @keywords internal
#'
#'
.gelman_diag <- function(x, start = 1, end = NA, conf = 0.95, digits = 2) {
    # Common parameter setup

    # start <- 1
    # end <- NA
    # conf <- 0.95
    # digits <- 2
    # x <- fit
    if (is.na(end)) {
        end <- slot(x, "nmc")
    }
    if (slot(x, "nchain") < 3) stop("Minimum three chains needed")

    iter <- start:end
    niter <- length(iter)
    nchain <- slot(x, "nchain")
    pnames <- slot(x, "pnames")

    # Regular posterior case
    data_list <- list(theta = slot(x, "theta")[, , iter, drop = FALSE])

    # Process each component (location/scale or just theta)
    results <- lapply(data_list, function(data) {
        # Calculate statistics
        wc_stats <- .calculate_within_chain_stats(data)
        xbar <- apply(data, c(1, 2), mean)
        B <- .calculate_between_chain_stats(xbar, niter)

        # Compute diagnostics
        psrf_components <- .compute_psrf_components(wc_stats$W, B, niter, nchain, conf)
        mpsrf <- .get_multivariate_psrf(dim(data)[1], niter, wc_stats$W, B, TRUE)

        # Prepare output
        psrf <- cbind(
            sqrt(psrf_components$df_adj * psrf_components$R2_estimate),
            sqrt(psrf_components$df_adj * psrf_components$R2_upper)
        )
        dimnames(psrf) <- list(pnames, c("Point est.", "Upper C.I."))

        list(psrf = psrf, mpsrf = mpsrf)
    })
    
    return(results$theta)
}

## S4 Method Definitions ------------------------------------------------------

#' Calculate Brook-Gelman diagnostic (R-hat) for MCMC convergence
#'
#' @param x A \code{posterior} object or list of posterior samples
#' @param start First iteration to use
#' @param end Last iteration to use (NA for all)
#' @param conf Confidence level for upper CI (default = 0.95)
#' @param digits Number of digits to print (default = 2)
#' @return A list containing psrf and mpsrf (if requested)
#' @rdname gelman
#' @examples
#' \dontrun{
#' fits <- fits1
#' phi <- RebuildHyper(fits)
#' thetas <- RebuildPosteriors(fits)
#' hat <- gelman(phi)
#' hat <- gelman(fits[[1]]$phi)
#' hat <- gelman(fits[[2]]$phi)
#' hat <- gelman(fits[[3]]$phi)
#' }
#' @export
setGeneric(
    "gelman",
    function(x, start = 1, end = NA, conf = 0.95,
             digits = 2) {
        standardGeneric("gelman")
    }
)

#' @rdname gelman
#' @export
setMethod(
    "gelman", "posterior",
    function(x, start = 1, end = NA, conf = 0.95, 
             digits = 2) {
        out <- tryCatch(
            {
                .gelman_diag(x, start, end, conf, digits)
            },
            error = function(e) {
                message("Error: chains not converged?\n")
                NULL
            }
        )
        return(out)
    }
)



