######### Helper functions ---------------------------------------------------

.sumlog <- function(x, debug) {
    sum(vapply(x, function(xi) {
        if (any(xi <= 0, na.rm = TRUE)) {
            if (debug) {
                warning("Invalid xi <= 0 or NA/NaN detected. Replacing with Machine eps.")
            }
            xi <- pmax(xi, .Machine$double.eps, na.rm = TRUE)
        }
        sum(log(xi))
    }, numeric(1)))
}

.update_prior_p0_p1 <- function(prior_list, param_vec) {
    n_params <- length(prior_list)

    if (length(param_vec) != 2 * n_params) {
        stop("param_vec must be twice as long as the number of parameters in prior_list")
    }

    p0_vals <- param_vec[1:n_params]
    p1_vals <- param_vec[(n_params + 1):(2 * n_params)]

    param_names <- names(prior_list)

    for (i in seq_along(param_names)) {
        prior_list[[param_names[i]]]$p0 <- p0_vals[i]
        prior_list[[param_names[i]]]$p1 <- p1_vals[i]
    }

    return(prior_list)
}

#' @importFrom methods slot
#' @importFrom ggdmcPrior dprior
.sumloghlike <- function(p_prior, phi_parameters, subject_theta, chain_idx, nmc_idx) {
    new_p_prior <- .update_prior_p0_p1(p_prior, phi_parameters)
    nsubject <- length(subject_theta)
    out <- numeric(nsubject)

    for (subject_idx in seq_len(nsubject)) {
        theta <- slot(subject_theta[[subject_idx]], "theta")
        theta_parameter <- theta[, chain_idx, nmc_idx]
        out[subject_idx] <- sum(dprior(new_p_prior, theta_parameter))
    }
    sum(out)
}

#' Set Theta Storage for Model
#'
#' Creates and configures a theta container for use in model fitting.
#'
#' @param theta_input An S4 object containing parameters specifying
#' the dimension and other information for the storage.
#'
#' @return An S4 object of class `posterior` containing:
#' \itemize{
#'   \item The configured theta parameters
#'   \item Additional metadata needed for posterior calculations
#' }
#'
#' @details
#' This function takes the configuration parameter needed to create a posterior
#'
#' @examples
#' \dontrun{
#' result <- set_up_new_samples(theta_input)
#' }
#' @importFrom methods new
#' @export
set_up_new_samples <- function(theta_input) {
    nparameter <- theta_input@nparameter
    nchain <- theta_input@nchain
    nmc <- theta_input@nmc

    out <- new("posterior",
        theta = array(NaN, dim = c(nparameter, nchain, nmc)),
        summed_log_prior = matrix(-Inf, nrow = nchain, ncol = nmc),
        log_likelihoods = matrix(-Inf, nrow = nchain, ncol = nmc),
        start = 1L,
        npar = nparameter,
        pnames = theta_input@pnames,
        nmc = nmc,
        thin = theta_input@thin,
        nchain = nchain
    )
}


#' Initialise Theta Parameters for MCMC Sampling
#'
#' Creates initial parameter values (theta) for MCMC sampling by combining prior distributions
#' with likelihood information from the data model instance (DMI).
#'
#' @param theta_input An S4 object of class `ThetaInput` containing parameter specifications
#' @param priors a `prior` object from ggdmcPrior carrying prior distributions.
#' @param dmi Data model instance containing data and experiment design specifications.
#' @param nmc_idx Iteration index for MCMC chain initialisation (default: 1)
#' @param seed Optional random seed for reproducible initialisation (default: NULL)
#' @param verbose Logical flag for printing initialisation details (default: FALSE)
#'
#' @return An S4 object of class `Theta` containing:
#' \itemize{
#'   \item Initial parameter values drawn from priors and adjusted by likelihood
#'   \item Metadata about initialisation process
#'   \item References to the priors and DMI used
#' }
#'
#' @details
#' This initialisation function:
#' \itemize{
#'   \item Draws initial values from specified prior distributions using `ggdmcPrior::rprior`
#'   \item Optionally adjusts values based on likelihood computed via
#'         ggdmcLikelihood::compute_subject_likelihood
#'   \item Handles both fixed and free parameters according to theta_input specifications
#'   \item Ensures parameters remain within valid bounds during initialisation
#' }
#'
#' When `verbose=TRUE`, prints internal diagnostic information about:
#' \itemize{
#'   \item Which parameters were initialised from priors
#'   \item Any likelihood-based adjustments made parameter bounds checking results
#' }
#'
#' @examples
#' \dontrun{
#' # Basic initialisation
#' theta <- initialise_theta(theta_input, priors, dmi)
#'
#' # Reproducible initialisation with seed
#' theta <- initialise_theta(theta_input, priors, dmi, seed = 123)
#'
#' # Verbose initialisation showing details
#' theta <- initialise_theta(theta_input, priors, dmi, verbose = TRUE)
#' }
#'
#' @export
#' @importFrom ggdmcPrior rprior dprior
#' @importFrom ggdmcLikelihood compute_subject_likelihood
initialise_theta <- function(
    theta_input, priors, dmi, nmc_idx = 1L,
    seed = NULL, verbose = FALSE) {
    if (is.list(dmi) && isS4(dmi[[1]])) {
        message("Use the first instance in the dmi list.")
        dmi <- dmi[[1]]
    }

    # Generate seeds based on number of cores
    if (is.null(seed)) {
        seed <- as.integer(Sys.time()) %% 1000000L
    }
    if (verbose) {
        message("Set seed to: ", seed)
        set.seed(seed)
    }

    if (is.list(priors)) {
        stop("Did you enter a S4 'priors' or a p_prior list?")
    }

    theta <- set_up_new_samples(theta_input)
    nchain <- theta_input@nchain
    max_attempts <- theta_input@max_init_attempts

    for (chain_idx in seq_len(nchain)) {
        valid_found <- FALSE

        for (attempt in seq_len(max_attempts)) {
            tmp_parameter <- ggdmcPrior::rprior(priors@p_prior, n = 1)

            # Calculate prior and likelihood
            slp <- sum(ggdmcPrior::dprior(priors@p_prior, tmp_parameter[, 1]))
            sll <- tryCatch(
                {
                    likelihoods <- ggdmcLikelihood::compute_subject_likelihood(dmi, tmp_parameter[, 1])
                    sum(sapply(likelihoods, .sumlog, verbose))
                },
                error = function(e) {
                    if (verbose) message("Likelihood error: ", conditionMessage(e))
                    -Inf
                }
            )


            if (is.finite(slp) && is.finite(sll)) {
                theta@summed_log_prior[chain_idx, nmc_idx] <- slp
                theta@log_likelihoods[chain_idx, nmc_idx] <- sll
                theta@theta[, chain_idx, nmc_idx] <- tmp_parameter[, 1]
                valid_found <- TRUE

                if (verbose) {
                    message(sprintf("Chain %d: Valid theta found on attempt %d", chain_idx, attempt))
                }
                break
            }
        }
        if (!valid_found) {
            stop(sprintf("Chain %d: Failed to find valid theta after %d attempts", chain_idx, max_attempts))
        }
    }

    return(theta)
}

#' Initialize Hierarchical Parameters (Phi) for MCMC Sampling
#'
#' Creates initial values for hierarchical parameters (phi) by combining population-level
#' priors with subject-specific parameter estimates.
#'
#' @param theta_input An S4 object of class `ThetaInput` containing individual-level parameter specifications
#' @param priors List of prior distribution specifications for population-level parameters
#' @param dmis List of data model interface objects (one per subject)
#' @param nmc_idx Iteration index for MCMC chain initialisation (default: 1)
#' @param seed Optional random seed for reproducible initialisation (default: NULL)
#' @param verbose Logical flag for printing initialisation details (default: FALSE)
#'
#' @return An S4 object of class `Phi` containing:
#' \itemize{
#'   \item Population-level parameter values drawn from hierarchical priors
#'   \item Individual-level parameter variations
#'   \item Metadata about the initialisation process
#' }
#'
#' @details
#' This initialisation function:
#' \itemize{
#'   \item Draws population-level parameters from hierarchical priors using `ggdmcPrior::rprior`
#'   \item Generates subject-specific variations around population values
#'   \item Ensures parameters respect the hierarchical structure specified in theta_input
#'   \item Performs bounds checking on all generated parameters
#' }
#'
#' When `verbose=TRUE`, prints:
#' \itemize{
#'   \item Population parameter initialisation values
#'   \item Summary statistics of subject-level variations
#'   \item Any parameter bounds adjustments made
#' }
#'
#' @examples
#' \dontrun{
#' # Basic initialisation
#' phi <- initialise_phi(theta_input, priors, dmis)
#'
#' # Reproducible initialisation with seed
#' phi <- initialise_phi(theta_input, priors, dmis, seed = 123)
#'
#' # Verbose initialisation showing details
#' phi <- initialise_phi(theta_input, priors, dmis, verbose = TRUE)
#' }
#'
#' @export
#' @importFrom ggdmcPrior rprior dprior
#' @export
initialise_phi <- function(
    theta_input, priors, dmis, nmc_idx = 1L,
    seed = NULL, verbose = FALSE) {
    if (is.null(seed)) {
        seed <- as.integer(Sys.time()) %% 1000000L
    }
    if (verbose) {
        message("Set seed to: ", seed)
        set.seed(seed)
    }
    if (is.null(priors@h_prior)) {
        stop("hyper prior is NULL")
    }

    p_prior <- priors@p_prior
    phi_prior <- priors@h_prior
    nsubject <- as.integer(length(dmis))
    subject_theta <- vector("list", nsubject)

    phi_input <- theta_input
    theta_input@pnames <- names(priors@p_prior)
    theta_input@nparameter <- as.integer(length(names(priors@p_prior)))

    if (verbose) cat("Processing subjects:\n")
    for (subject_idx in seq_len(nsubject)) {
        if (verbose) cat("  Subject", subject_idx, "\n")
        subject_seed <- seed + subject_idx

        subject_theta[[subject_idx]] <- initialise_theta(
            theta_input,
            priors, dmis[[subject_idx]], nmc_idx, subject_seed, verbose
        )
    }

    names(subject_theta) <- names(dmis)
    phi <- set_up_new_samples(phi_input)
    nchain <- phi_input@nchain
    max_attempts <- phi_input@max_init_attempts

    for (chain_idx in seq_len(nchain)) {
        valid_found <- FALSE

        for (attempt in seq_len(max_attempts)) {
            tmp_parameter <- rprior(phi_prior, n = 1)
            slp <- sum(dprior(phi_prior, tmp_parameter[, 1]))

            sll <- tryCatch(
                {
                    .sumloghlike(p_prior, tmp_parameter[, 1], subject_theta, chain_idx, nmc_idx)
                },
                error = function(e) {
                    if (verbose) message(sprintf("Chain %d, attempt %d: Likelihood error: %s", chain_idx, attempt, e$message))
                    -Inf
                }
            )

            if (is.finite(slp) && is.finite(sll)) {
                phi@summed_log_prior[chain_idx, nmc_idx] <- slp
                phi@log_likelihoods[chain_idx, nmc_idx] <- sll
                phi@theta[, chain_idx, nmc_idx] <- tmp_parameter[, 1]


                valid_found <- TRUE
                if (verbose) {
                    message(sprintf("Chain %d: Valid phi found on attempt %d", chain_idx, attempt))
                }
                break
            }
        }

        if (!valid_found) {
            stop(sprintf("Chain %d: Failed to find valid phi after %d attempts", chain_idx, max_attempts))
        }
    }

    return(list(phi = phi, subject_theta = subject_theta))
}
