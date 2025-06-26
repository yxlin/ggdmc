.get_likelihood_matrix <- function(x, is_prior = FALSE, use_phi = FALSE) {
    if (length(x) == 0) stop("Input list is empty")

    extract_matrix <- function(obj) {
        source <- if (use_phi) obj$phi else obj
        mat <- if (is_prior) source@summed_log_prior else source@log_likelihoods
        t(mat)
    }

    mats <- lapply(x, extract_matrix)
    do.call(rbind, lapply(mats, as.vector))
}


.GetParamArray <- function(fits_list, param_type = c("theta", "phi")) {
    if (length(fits_list) == 0) {
        stop("Input list is empty")
    }

    # Validate parameter type
    param_type <- match.arg(param_type)

    # Extract the first element to get reference dimensions and names
    if (param_type == "phi") {
        fit1 <- fits_list[[1]]$phi
    } else {
        fit1 <- fits_list[[1]]
    }

    npar <- fit1@npar
    nchain <- fit1@nchain
    nmc <- fit1@nmc
    pnames <- fit1@pnames

    # Get number of chains (elements in the input list)
    new_nchain <- length(fits_list)

    # Initialise array to store combined results
    combined_array <- array(dim = c(npar, new_nchain, nchain * nmc))

    # Fill the array
    for (i in seq_along(fits_list)) {
        # Extract the appropriate array based on parameter type
        if (param_type == "phi") {
            theta_array <- fits_list[[i]]$phi@theta
        } else {
            theta_array <- fits_list[[i]]@theta
        }

        # Reshape and store
        dim(theta_array) <- c(npar, nchain * nmc)
        combined_array[, i, ] <- theta_array
    }

    return(combined_array)
}


#' @importFrom methods new
#' @rdname rebuild_posteriors
#' @export
RebuildPosterior <- function(x) {
    theta_array <- .GetParamArray(x, param_type = "theta")

    summed_log_prior <- .get_likelihood_matrix(x, is_prior = TRUE, use_phi = FALSE)
    log_likelihood <- .get_likelihood_matrix(x, is_prior = FALSE, use_phi = FALSE)

    x1 <- x[[1]]
    new_nmc <- x1@nmc * x1@nchain

    new("posterior",
        theta = theta_array,
        summed_log_prior = summed_log_prior,
        log_likelihoods = log_likelihood,
        start = as.integer(x1@start),
        npar = as.integer(x1@npar),
        pnames = x1@pnames,
        nmc = as.integer(new_nmc),
        thin = as.integer(x1@thin),
        nchain = as.integer(length(x))
    )
}

#' @importFrom methods new
#' @rdname rebuild_posteriors
#' @export
RebuildHyper <- function(x) {
    phi_array <- .GetParamArray(x, param_type = "phi")

    new_hyper_summed_log_prior <- .get_likelihood_matrix(x, TRUE, TRUE)
    new_hyper_log_likelihood <- .get_likelihood_matrix(x, FALSE, TRUE)

    x1 <- x[[1]]
    new_nmc <- x1$phi@nmc * x1$phi@nchain

    new(
        "posterior",
        theta = phi_array,
        summed_log_prior = new_hyper_summed_log_prior,
        log_likelihoods = new_hyper_log_likelihood,
        start = as.integer(x1$phi@start),
        npar = as.integer(x1$phi@npar),
        pnames = x1$phi@pnames,
        nmc = as.integer(new_nmc),
        thin = as.integer(x1$phi@thin),
        nchain = as.integer(length(x))
    )
}

#' Rebuild Posterior Objects from MCMC Chains
#'
#' These functions reconstruct posterior objects by
#' combining multiple DE-MC chains.
#'
#' @param x A list of model fit objects (typically of class `posterior`)
#'   containing MCMC chains to be combined.
#'
#' @return
#' \describe{
#'   \item{RebuildPosterior}{Returns a single `posterior` object combining
#' all DE-MC chains from the theta object}
#'   \item{RebuildHyper}{Returns a single `posterior` object from hierarchical
#' model fits}
#'   \item{RebuildPosteriors}{Returns a list of `posterior` objects (one per subject)}
#' }
#'
#' @details
#' These functions are used to:
#' \itemize{
#'   \item Combine multiple DE-MC chains (e.g., from genetic samplers)
#'   \item Reconstruct proper posterior objects after manipulation
#'   \item Maintain chain information while aggregating results
#' }
#'
#' The functions handle:
#' \itemize{
#'   \item Combining theta/phi parameters across DE-MC chains
#'   \item Aggregating log-prior and log-likelihood matrices
#'   \item Properly setting dimensions (nmc, nchain, etc.)
#' }
#'
#' @section Notes:
#' \itemize{
#'   \item Input objects should have consistent dimensions across chains
#'   \item The `thin` parameter is taken from the first chain
#'   \item For hierarchical models, use `RebuildPosteriors` for subject-level
#' parameters and `RebuildHyper` for group-level parameters
#' }
#'
#' @examples
#' \dontrun{
#' # After running multiple MCMC chains:
#' fits <- fits1
#' phi <- RebuildHyper(fits)
#' thetas <- RebuildPosteriors(fits)
#' }
#'
#' @importFrom methods new
#' @rdname rebuild_posteriors
#' @aliases RebuildPosterior RebuildHyper RebuildPosteriors
#' @export
RebuildPosteriors <- function(x) {
    # Input validation
    if (!is.list(x) || length(x) == 0) {
        stop("Input must be a non-empty list of model fits")
    }

    # Get basic information from first chain and first subject
    first_chain <- x[[1]]
    first_subject <- first_chain$subject_theta[[1]]
    nsubject <- length(first_chain$subject_theta)
    nchain <- length(x)

    # Pre-allocate the output list
    posteriors <- vector("list", nsubject)

    # Process each subject
    for (i in seq_len(nsubject)) {
        # Extract all chains for this subject
        subject_chains <- lapply(x, function(chain) chain$subject_theta[[i]])

        # Combine thetas from all chains
        theta_array <- .GetParamArray(subject_chains, param_type = "theta")

        # Get priors and likelihoods
        summed_log_prior <- .get_likelihood_matrix(subject_chains, TRUE, FALSE)
        log_likelihood <- .get_likelihood_matrix(subject_chains, FALSE, FALSE)

        # Calculate new nmc
        new_nmc <- first_subject@nmc * first_subject@nchain

        # Create new posterior object
        posteriors[[i]] <- new(
            "posterior",
            theta = theta_array,
            summed_log_prior = summed_log_prior,
            log_likelihoods = log_likelihood,
            start = as.integer(first_subject@start),
            npar = as.integer(first_subject@npar),
            pnames = first_subject@pnames,
            nmc = as.integer(new_nmc),
            thin = as.integer(first_subject@thin),
            nchain = as.integer(nchain)
        )
    }

    return(posteriors)
}
