.get_nchain <- function(priors, nchain = NULL) {
    nparameter <- priors@nparameter
    pnames <- priors@pnames

    if (is.null(nchain)) {
        nchain <- nparameter * 3.0
        message("Use ", nchain, " chains in an instance to optimise the model")
    }
    return(nchain)
}


parallel_lapply <- function(
    ncore,
    config_list,
    fun, # The function to run (run_hyper or run)
    ..., # Additional arguments needed by fun
    samples_list = NULL, # Optional samples list
    hyper_dmi = NULL, # Optional hyper_dmi
    dmis = NULL, # Optional dmis
    dmi = NULL # Optional dmi
    ) {
    seq_list <- seq_len(ncore)
    dots <- list(...)

    if (.Platform$OS.type == "unix" && ncore > 1) {
        # Unix-like systems (Linux, Mac) - use mclapply
        message("Using mclapply (fork-based parallelisation)")

        out <- parallel::mclapply(seq_list, function(i) {
            set.seed(config_list[[i]]@seed)

            tryCatch(
                {
                    # Dynamically call the function with appropriate args
                    if (identical(fun, run_hyper)) {
                        # run_hyper(config_list[[i]], hyper_dmi, samples_list[[i]]$phi)

                        fun(config_list[[i]], hyper_dmi, samples_list[[i]])
                    } else if (identical(fun, run)) {
                        # run(config_list[[i]], dmis, samples_list[[i]])
                        fun(config_list[[i]], dmis, samples_list[[i]])
                    } else if (identical(fun, run_subject)) {
                        fun(config_list[[i]], dmi, samples_list[[i]])
                    } else {
                        # Generic case for other functions
                        do.call(fun, c(list(config_list[[i]]), dots))
                    }
                },
                error = function(e) {
                    message("Chain ", i, " failed: ", e$message)
                    NULL
                }
            )
        }, mc.cores = ncore)
    } else if (ncore > 1) {
        # Windows - use parLapply with PSOCK cluster
        message("Using parLapply (PSOCK cluster)")
        cl <- parallel::makeCluster(ncore)
        on.exit(parallel::stopCluster(cl))

        # Export necessary objects
        parallel::clusterExport(cl,
            varlist = c(
                "config_list", "fun", "samples_list",
                "hyper_dmi", "dmis", "dmi", "dots"
            ),
            envir = environment()
        )


        # Set seeds on workers
        parallel::clusterApply(cl, seq_list, function(i) {
            set.seed(config_list[[i]]@seed)
        })

        out <- parallel::parLapply(cl, seq_list, function(i) {
            tryCatch(
                {
                    if (identical(fun, run_hyper)) {
                        fun(config_list[[i]], hyper_dmi, samples_list[[i]]$phi)
                    } else if (identical(fun, run)) {
                        fun(config_list[[i]], dmis, samples_list[[i]])
                    } else if (identical(fun, run_subject)) {
                        fun(config_list[[i]], dmi, samples_list[[i]])
                    } else {
                        do.call(fun, c(list(config_list[[i]]), dots))
                    }
                },
                error = function(e) {
                    message("Chain ", i, " failed: ", e$message)
                    NULL
                }
            )
        })
    } else {
        # Sequential execution
        message("Running sequentially (ncore = 1)")
        out <- lapply(seq_list, function(i) {
            set.seed(config_list[[i]]@seed)
            tryCatch(
                {
                    if (identical(fun, run_hyper)) {
                        fun(config_list[[i]], hyper_dmi, samples_list[[i]]$phi)
                    } else if (identical(fun, run)) {
                        fun(config_list[[i]], dmis, samples_list[[i]])
                    } else if (identical(fun, run_subject)) {
                        fun(config_list[[i]], dmi, samples_list[[i]])
                    } else {
                        do.call(fun, c(list(config_list[[i]]), dots))
                    }
                },
                error = function(e) {
                    message("Chain ", i, " failed: ", e$message)
                    NULL
                }
            )
        })
    }
    return(out)
}


#' Initialise or Continue MCMC Sampling
#'
#' These functions handle initialisation and continuation of Markov Chain Monte
#' Carlo (MCMC) sampling for hierarchical models, individual subjects, and
#' (merely) hyperparameter estimation (i.e., standard models).
#'
#' @param dmis For hierarchical models: A list of data model instances (DMIs)
#' for each subject
#' @param dmi For single-subject models: A DMI instance for one subject
#' @param hyper_dmi For standard models: a DMI for hyperparameters using
#' `type`, `hyper`)
#' @param priors A list of prior distributions
#' @param samples_list Either:
#'   \itemize{
#'     \item For Start* functions: Optional the user-initislised samples.
#'     \item For Restart* functions: Required previous samples object to
#' continue sampling
#'   }
#' @param nmc Number of MCMC iterations to run (default = 500)
#' @param nchain Number of independent chains (default = NULL for automatic
#' determination; not used in RestartSampling_subject)
#' @param thin Thinning interval (default = 1)
#' @param report_length Frequency of progress reports in iterations (default = 100)
#' @param max_init_attempts Maximum attempts to find valid initial values (default = 1000)
#' @param is_print Whether to print progress messages (default = TRUE)
#' @param pop_migration_prob Population-level chain migration probability (default = 0)
#' @param sub_migration_prob Subject-level chain migration probability (default = 0)
#' @param gamma_precursor Scaling factor for proposal distribution (default = 2.38)
#' @param rp Random perturbation factor to add noise to DE-MC proposal. The large the value,
#' the big the proposal will jump (default = 0.001)
#' @param is_hblocked Whether to use block updating for hyperparameters (default = FALSE)
#' @param is_pblocked Whether to use block updating for subject parameters (default = FALSE)
#' @param ncore Number of CPU cores for parallel computation (default = 3)
#' @param seed Random seed for reproducibility (default = NULL for random seed)
#' @param pop_debug Logical. Whether to print population level debugging information.
#' @param sub_debug Logical. Whether to print subject level debugging information.
#' @return An object containing MCMC samples, with class:
#' \itemize{
#'   \item `hierarchical_samples` for hierarchical versions
#'   \item `subject_samples` for single-subject versions
#'   \item `hyper_samples` for hyperparameter versions
#' }
#'
#' @details
#' \subsection{Function Types}{
#' \itemize{
#'   \item \strong{Hierarchical Models}:
#'   \itemize{
#'     \item `StartSampling`: initialise full hierarchical sampling
#'     \item `RestartSampling`: Continue hierarchical sampling
#'   }
#'   \item \strong{Single-Subject Models}:
#'   \itemize{
#'     \item `StartSampling_subject`: initialise single-subject sampling
#'     \item `RestartSampling_subject`: Continue single-subject sampling
#'   }
#'   \item \strong{Hyperparameter Models}:
#'   \itemize{
#'     \item `StartSampling_hyper`: initialise hyperparameter-only sampling
#'     \item `RestartSampling_hyper`: Continue hyperparameter-only sampling
#'   }
#' }}
#'
#' \subsection{Key Features}{
#' \itemize{
#'   \item Adaptive MCMC with optional chain migration
#'   \item Parallel chain execution
#'   \item Progress reporting and convergence monitoring
#'   \item Three-level modeling (hyperparameter-only, population, and individual
#'         subjects)
#' }}
#'
#' @section Migration Parameters:
#' The migration sampler helps to improve chain mixing at different levels:
#' \itemize{
#'   \item `pop_migration_prob`: For population-level parameters
#'   \item `sub_migration_prob`: For subject-level parameters and
#' hyper-only parameters.
#' }. The sampler essentially compares the self chain against its neighbouring
#' chains, and updates the self chain when the neighboring chain results
#' in higher likelihood. The user must be aware that such a strategy may
#' result in an optimisation process falls into a local maximum (likelihood).
#'
#' @examples
#' \dontrun{
#' # Hierarchical models
#' hier_fits0 <- StartSampling(pop_dmis, pop_priors,
#'     sub_migration_prob = 0.05,
#'     thin = 8L, pop_debug = F, seed = 9032
#' )
#'
#' hier_fits1 <- RestartSampling(hier_fits0,
#'     pop_migration_prob = 0.02,
#'     sub_migration_prob = 0.00,
#'     thin = 4L, seed = 9032
#' )
#'
#' # Single-subject models
#' subj_fits0 <- StartSampling_subject(sub_dmis[[1]], sub_priors,
#'     sub_migration_prob = 0.02,
#'     thin = 2, seed = 9032
#' )
#'
#' subj_fits1 <- RestartSampling_subject(subj_fits0,
#'     sub_migration_prob = 0.00, thin = 4, seed = 9032
#' )
#'
#' # Hyperparameter models
#' hyper_fits0 <- StartSampling_hyper(hyper_dmi, pop_dmis, pop_priors,
#'     sub_migration_prob = 0.05, thin = 4
#' )
#'
#' hyper_fits1 <- RestartSampling_hyper(fits0, sub_migration_prob = 0.00, thin = 2)
#' }
#'
#' @seealso
#' \code{\link{RebuildPosterior}} for combining chains
#'
#' @rdname sampling_functions
#' @aliases StartSampling RestartSampling StartSampling_subject RestartSampling_subject
#'          StartSampling_hyper RestartSampling_hyper
#' @export
StartSampling <- function(
    dmis, priors, samples_list = NULL,
    nmc = 500L, nchain = NULL, thin = 1L,
    report_length = 100L, max_init_attempts = 1000L,
    is_print = TRUE,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    gamma_precursor = 2.38,
    rp = 0.001, is_hblocked = FALSE, is_pblocked = FALSE,
    pop_debug = FALSE, sub_debug = FALSE,
    ncore = 3L, seed = NULL) {
    nparameter <- priors@nparameter
    pnames <- priors@pnames
    nchain <- .get_nchain(priors, nchain)

    theta_input <- setThetaInput(
        nmc = nmc, nchain = nchain, thin = thin, pnames = pnames,
        report_length = report_length, max_init_attempts = max_init_attempts,
        is_print = is_print
    )

    de_input <- setDEInput(
        pop_migration_prob = pop_migration_prob, sub_migration_prob = sub_migration_prob, gamma_precursor = gamma_precursor,
        rp = rp, is_hblocked = is_hblocked, is_pblocked = is_pblocked,
        nparameter = as.integer(nparameter), nchain = as.integer(nchain),
        pop_debug = pop_debug, sub_debug = sub_debug
    )

    config_list <- set_configs(
        prior = priors, theta_input = theta_input, de_input = de_input,
        ncore = ncore, seed = seed
    )

    if (is.null(samples_list)) {
        message("Initialise ", ncore, " independent sets of new phi samples")
        samples_list <- lapply(seq_len(ncore), function(i) {
            initialise_phi(theta_input, priors, dmis, seed = config_list[[i]]@seed)
        })
    }

    message("\nLaunching ", ncore, " cores/chains (mclapply)")
    t0 <- Sys.time()
    out <- parallel_lapply(
        ncore = ncore, config_list = config_list,
        fun = run,
        samples_list = samples_list, dmis = dmis
    )
    t1 <- Sys.time()

    proc_time <- difftime(t1, t0, units = "secs")[[1]]
    message("Processing time: ", round(proc_time, 2), " secs.")

    attr(out, "configs") <- config_list
    attr(out, "dmis") <- dmis
    attr(out, "priors") <- priors
    attr(out, "start_samples") <- samples_list
    return(out)
}

#' @rdname sampling_functions
#' @export
RestartSampling <- function(
    samples_list, nmc = 500L, nchain = NULL, thin = 1L,
    report_length = 100L, max_init_attempts = 1000L,
    is_print = TRUE,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00,
    gamma_precursor = 2.38,
    rp = 0.001, is_hblocked = FALSE, is_pblocked = FALSE,
    pop_debug = FALSE, sub_debug = FALSE,
    seed = NULL) {
    if (missing(samples_list)) {
        stop("Must provide a previously fitted samples")
    }

    # Get the original configs and other attributes from the samples
    old_configs <- attr(samples_list, "configs")
    dmis <- attr(samples_list, "dmis")
    priors <- attr(samples_list, "priors")
    ncore <- as.integer(length(old_configs))

    # Get old parameters
    old_params <- list(
        nmc = as.integer(old_configs[[1]]@theta_input@nmc),
        nchain = as.integer(old_configs[[1]]@theta_input@nchain),
        thin = as.integer(old_configs[[1]]@theta_input@thin),
        report_length = as.integer(old_configs[[1]]@theta_input@report_length),
        max_init_attempts = as.integer(old_configs[[1]]@theta_input@max_init_attempts),
        is_print = old_configs[[1]]@theta_input@is_print,
        pop_migration_prob = old_configs[[1]]@de_input@pop_migration_prob,
        sub_migration_prob = old_configs[[1]]@de_input@sub_migration_prob,
        gamma_precursor = old_configs[[1]]@de_input@gamma_precursor,
        rp = old_configs[[1]]@de_input@rp,
        is_hblocked = old_configs[[1]]@de_input@is_hblocked,
        is_pblocked = old_configs[[1]]@de_input@is_pblocked,
        pop_debug = old_configs[[1]]@de_input@pop_debug,
        sub_debug = old_configs[[1]]@de_input@sub_debug,
        ncore = ncore
    )

    if (is.null(nchain)) {
        nchain <- as.integer(3 * priors@nparameter)
    }

    new_params <- list(
        nmc = as.integer(nmc),
        nchain = as.integer(nchain),
        thin = as.integer(thin),
        report_length = as.integer(report_length),
        max_init_attempts = as.integer(max_init_attempts),
        is_print = is_print,
        pop_migration_prob = pop_migration_prob,
        sub_migration_prob = sub_migration_prob,
        gamma_precursor = gamma_precursor,
        rp = rp,
        is_hblocked = is_hblocked,
        is_pblocked = is_pblocked,
        pop_debug = pop_debug,
        sub_debug = sub_debug,
        ncore = ncore
    )

    # Function to compare and report differences
    compare_configs <- function(old, new) {
        changed <- FALSE
        all_names <- union(names(old), names(new))

        for (name in all_names) {
            if (!name %in% names(old)) {
                message(sprintf("New parameter added: %s = %s", name, new[[name]]))
                changed <- TRUE
            } else if (!name %in% names(new)) {
                message(sprintf("Parameter removed: %s (was %s)", name, old[[name]]))
                changed <- TRUE
            } else if (!identical(old[[name]], new[[name]])) {
                message(sprintf(
                    "Parameter changed: %s (old: %s, new: %s)",
                    name, old[[name]], new[[name]]
                ))
                changed <- TRUE
            }
        }
        return(changed)
    }

    # Compare and report differences
    message("Comparing configuration parameters:")
    need_new_configs <- compare_configs(old_params, new_params)

    if (need_new_configs) {
        message("\nConfiguration parameters changed - creating new configs")


        theta_input <- setThetaInput(
            nmc = nmc, nchain = nchain, thin = thin, pnames = priors@pnames,
            report_length = report_length, max_init_attempts = max_init_attempts,
            is_print = is_print
        )

        de_input <- setDEInput(
            pop_migration_prob = pop_migration_prob,
            sub_migration_prob = sub_migration_prob,
            gamma_precursor = gamma_precursor, rp = rp,
            is_hblocked = is_hblocked, is_pblocked = is_pblocked,
            nparameter = as.integer(priors@nparameter), nchain = as.integer(nchain),
            pop_debug = pop_debug, sub_debug = sub_debug
        )
        config_list <- set_configs(
            prior = priors, theta_input = theta_input, de_input = de_input,
            ncore = ncore, seed = seed
        )
    } else {
        message("\nAll configuration parameters match - using existing configs")
        config_list <- old_configs
    }

    message("\nSampling ", ncore, " independent chains")

    t0 <- Sys.time()
    out <- parallel_lapply(
        ncore = ncore, config_list = config_list,
        fun = run, samples_list = samples_list, dmis = dmis
    )
    t1 <- Sys.time()

    proc_time <- difftime(t1, t0, units = "secs")[[1]]
    message("Processing time: ", round(proc_time, 2), " secs.")

    # Preserve attributes from original samples
    attr(out, "configs") <- config_list
    attr(out, "dmis") <- dmis
    attr(out, "priors") <- priors

    return(out)
}


#' @rdname sampling_functions
#' @export
StartSampling_subject <- function(
    dmi, priors, samples_list = NULL,
    nmc = 500L, nchain = NULL, thin = 1L,
    report_length = 100L, max_init_attempts = 1000L,
    is_print = TRUE,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00, gamma_precursor = 2.38,
    rp = 0.001,
    is_hblocked = FALSE, is_pblocked = FALSE,
    ncore = 3L, seed = NULL) {
    nparameter <- priors@nparameter
    pnames <- priors@pnames
    nchain <- .get_nchain(priors, nchain)

    theta_input <- setThetaInput(
        nmc = nmc, nchain = nchain, thin = thin, pnames = pnames,
        report_length = report_length, max_init_attempts = max_init_attempts,
        is_print = is_print
    )

    de_input <- setDEInput(
        pop_migration_prob = pop_migration_prob,
        sub_migration_prob = sub_migration_prob,
        gamma_precursor = gamma_precursor, rp = rp,
        is_hblocked = is_hblocked, is_pblocked = is_pblocked,
        nparameter = as.integer(priors@nparameter), nchain = as.integer(nchain)
    )
    config_list <- set_configs(
        prior = priors, theta_input = theta_input, de_input = de_input,
        ncore = ncore, seed = seed
    )

    if (is.null(samples_list)) {
        message("Initialise ", ncore, " independent sets of new theta samples")

        samples_list <- lapply(seq_len(ncore), function(i) {
            initialise_theta(theta_input, priors, dmi, seed = config_list[[i]]@seed)
        })
    }


    message("\nLaunching ", ncore, " independent chains")
    t0 <- Sys.time()
    out <- parallel_lapply(ncore = ncore, config_list = config_list, fun = run_subject, samples_list = samples_list, dmi = dmi)
    t1 <- Sys.time()

    proc_time <- difftime(t1, t0, units = "secs")[[1]]
    message("Processing time (in R): ", round(proc_time, 2), " secs.")

    attr(out, "configs") <- config_list
    attr(out, "dmi") <- dmi
    attr(out, "priors") <- priors
    attr(out, "start_samples") <- samples_list


    return(out)
}


#' @rdname sampling_functions
#' @export
RestartSampling_subject <- function(
    samples_list,
    nmc = 500L, thin = 1L,
    report_length = 100L, max_init_attempts = 1000L,
    is_print = TRUE,
    pop_migration_prob = 0.00,
    sub_migration_prob = 0.00, gamma_precursor = 2.38,
    rp = 0.001,
    is_hblocked = FALSE, is_pblocked = FALSE,
    seed = NULL) {
    # Get the original configs and other attributes from the samples

    old_configs <- attr(samples_list, "configs")
    dmi <- attr(samples_list, "dmi")
    priors <- attr(samples_list, "priors")
    ncore <- length(old_configs)

    # Check if any configuration parameters have changed
    old_params <- list(
        nmc = old_configs[[1]]@theta_input@nmc,
        nchain = old_configs[[1]]@theta_input@nchain,
        thin = old_configs[[1]]@theta_input@thin,
        report_length = old_configs[[1]]@theta_input@report_length,
        max_init_attempts = old_configs[[1]]@theta_input@max_init_attempts,
        is_print = old_configs[[1]]@theta_input@is_print,
        pop_migration_prob = old_configs[[1]]@de_input@pop_migration_prob,
        sub_migration_prob = old_configs[[1]]@de_input@sub_migration_prob,
        gamma_precursor = old_configs[[1]]@de_input@gamma_precursor,
        rp = old_configs[[1]]@de_input@rp,
        is_hblocked = old_configs[[1]]@de_input@is_hblocked,
        is_pblocked = old_configs[[1]]@de_input@is_pblocked,
        ncore = ncore
    )

    new_params <- list(
        nmc = nmc,
        nchain = old_configs[[1]]@theta_input@nchain,
        thin = thin,
        report_length = report_length,
        max_init_attempts = max_init_attempts,
        is_print = is_print,
        pop_migration_prob = pop_migration_prob,
        sub_migration_prob = sub_migration_prob,
        gamma_precursor = gamma_precursor,
        rp = rp,
        is_hblocked = is_hblocked,
        is_pblocked = is_pblocked,
        ncore = ncore
    )

    # Compare parameters to see if we need new configs
    need_new_configs <- !isTRUE(all.equal(old_params, new_params))
    # need_new_configs

    if (need_new_configs) {
        message("Configuration parameters changed - creating new configs")
        nchain <- old_configs[[1]]@theta_input@nchain

        nparameter <- priors@nparameter
        pnames <- priors@pnames


        theta_input <- setThetaInput(
            nmc = nmc, nchain = nchain, thin = thin, pnames = pnames,
            report_length = report_length, max_init_attempts = max_init_attempts,
            is_print = is_print
        )

        de_input <- setDEInput(
            pop_migration_prob = pop_migration_prob,
            sub_migration_prob = sub_migration_prob,
            gamma_precursor = gamma_precursor, rp = rp, is_hblocked = is_hblocked,
            is_pblocked = is_pblocked,
            nparameter = as.integer(nparameter), nchain = as.integer(nchain)
        )

        config_list <- set_configs(
            prior = priors, theta_input = theta_input, de_input = de_input,
            ncore = ncore, seed = seed
        )
    } else {
        message("Using existing configs")
        config_list <- old_configs
    }

    message("\nSampling ", ncore, " independent chains")

    t0 <- Sys.time()
    out <- parallel_lapply(
        ncore = ncore, config_list = config_list,
        fun = run_subject,
        samples_list = samples_list, dmi = dmi
    )

    t1 <- Sys.time()

    proc_time <- difftime(t1, t0, units = "secs")[[1]]
    message("Processing time (in R): ", round(proc_time, 2), " secs.")

    # Preserve attributes from original samples
    attr(out, "configs") <- config_list
    attr(out, "dmi") <- dmi
    attr(out, "priors") <- priors

    return(out)
}


#' @rdname sampling_functions
#' @export
StartSampling_hyper <- function(
    hyper_dmi, dmis, priors, samples_list = NULL,
    nmc = 500L, nchain = NULL, thin = 1L,
    report_length = 100L, max_init_attempts = 1000L,
    is_print = TRUE, pop_migration_prob = 0.00,
    sub_migration_prob = 0.00, gamma_precursor = 2.38,
    rp = 0.001, is_hblocked = FALSE, is_pblocked = FALSE,
    ncore = 3L, seed = NULL) {
    nparameter <- priors@nparameter # 10
    pnames <- priors@pnames
    nchain <- .get_nchain(priors, nchain)

    if (is.null(priors@h_prior)) {
        stop("h_prior has not found")
    }

    theta_input <- setThetaInput(
        nmc = nmc, nchain = nchain, thin = thin,
        pnames = pnames, report_length = report_length, max_init_attempts = max_init_attempts,
        is_print = is_print
    )

    de_input <- setDEInput(
        pop_migration_prob = pop_migration_prob,
        sub_migration_prob = sub_migration_prob,
        gamma_precursor = gamma_precursor, rp = rp, is_hblocked = is_hblocked,
        is_pblocked = is_pblocked,
        # different nparameter
        nparameter = as.integer(nparameter), nchain = as.integer(nchain)
    )
    config_list <- set_configs(
        prior = priors, theta_input = theta_input, de_input = de_input,
        ncore = ncore, seed = seed
    )

    if (is.null(samples_list)) {
        message("Initialise ", ncore, " independent sets of new phi samples")
        samples_list <- lapply(seq_len(ncore), function(i) {
            tmp_list <- initialise_phi(theta_input, priors, dmis, seed = config_list[[i]]@seed)
            tmp_list$phi
        })
    }

    message("\nLaunching ", ncore, " independent chains")

    t0 <- Sys.time()
    out <- parallel_lapply(
        ncore = ncore, config_list = config_list,
        fun = run_hyper, hyper_dmi = hyper_dmi, samples_list = samples_list
    )
    t1 <- Sys.time()

    proc_time <- difftime(t1, t0, units = "secs")[[1]]
    message("Processing time: ", round(proc_time, 2), " secs.")

    attr(out, "configs") <- config_list
    attr(out, "hyper_dmi") <- hyper_dmi
    attr(out, "dmis") <- dmis
    attr(out, "priors") <- priors
    attr(out, "start_samples") <- samples_list
    return(out)
}

#' @rdname sampling_functions
#' @export
RestartSampling_hyper <- function(
    samples_list,
    nmc = 500L, nchain = NULL, thin = 1L,
    report_length = 100L, max_init_attempts = 1000L,
    is_print = TRUE, pop_migration_prob = 0.00,
    sub_migration_prob = 0.00, gamma_precursor = 2.38,
    rp = 0.001, is_hblocked = FALSE, is_pblocked = FALSE, seed = NULL) {
    if (missing(samples_list)) {
        stop("Must provide a previously fitted samples")
    }


    old_configs <- attr(samples_list, "configs")
    hyper_dmi <- attr(samples_list, "hyper_dmi")
    dmis <- attr(samples_list, "dmis")
    priors <- attr(samples_list, "priors")
    ncore <- as.integer(length(old_configs))

    # Check if any configuration parameters have changed
    old_params <- list(
        nmc = old_configs[[1]]@theta_input@nmc,
        nchain = old_configs[[1]]@theta_input@nchain,
        thin = old_configs[[1]]@theta_input@thin,
        report_length = old_configs[[1]]@theta_input@report_length,
        max_init_attempts = old_configs[[1]]@theta_input@max_init_attempts,
        is_print = old_configs[[1]]@theta_input@is_print,
        pop_migration_prob = old_configs[[1]]@de_input@pop_migration_prob,
        sub_migration_prob = old_configs[[1]]@de_input@sub_migration_prob,
        gamma_precursor = old_configs[[1]]@de_input@gamma_precursor,
        rp = old_configs[[1]]@de_input@rp,
        is_hblocked = old_configs[[1]]@de_input@is_hblocked,
        is_pblocked = old_configs[[1]]@de_input@is_pblocked,
        ncore = ncore
    )


    if (is.null(nchain)) {
        nchain <- 3 * priors@nparameter
    }

    new_params <- list(
        nmc = nmc,
        nchain = nchain,
        thin = thin,
        report_length = report_length,
        max_init_attempts = max_init_attempts,
        is_print = is_print,
        pop_migration_prob = pop_migration_prob,
        sub_migration_prob = sub_migration_prob,
        gamma_precursor = gamma_precursor,
        rp = rp,
        is_hblocked = is_hblocked,
        is_pblocked = is_pblocked,
        ncore = ncore
    )
    # Function to compare and report differences
    compare_configs <- function(old, new) {
        changed <- FALSE
        all_names <- union(names(old), names(new))

        for (name in all_names) {
            if (!name %in% names(old)) {
                message(sprintf("New parameter added: %s = %s", name, new[[name]]))
                changed <- TRUE
            } else if (!name %in% names(new)) {
                message(sprintf("Parameter removed: %s (was %s)", name, old[[name]]))
                changed <- TRUE
            } else if (!identical(old[[name]], new[[name]])) {
                message(sprintf(
                    "Parameter changed: %s (old: %s, new: %s)",
                    name, old[[name]], new[[name]]
                ))
                changed <- TRUE
            }
        }
        return(changed)
    }

    # Compare and report differences
    message("Comparing configuration parameters:")
    need_new_configs <- compare_configs(old_params, new_params)

    if (need_new_configs) {
        message("Configuration parameters changed - creating new configs")


        theta_input <- setThetaInput(
            nmc = nmc, nchain = nchain, thin = thin, pnames = priors@pnames,
            report_length = report_length, max_init_attempts = max_init_attempts,
            is_print = is_print
        )

        de_input <- setDEInput(
            pop_migration_prob = pop_migration_prob,
            sub_migration_prob = sub_migration_prob,
            gamma_precursor = gamma_precursor, rp = rp,
            is_hblocked = is_hblocked, is_pblocked = is_pblocked,
            nparameter = as.integer(priors@nparameter), nchain = as.integer(nchain)
        )

        config_list <- set_configs(
            prior = priors, theta_input = theta_input, de_input = de_input,
            ncore = ncore, seed = seed
        )
    } else {
        message("\nAll configuration parameters match - using existing configs")
        config_list <- old_configs
    }

    message("\nSampling ", ncore, " independent chains")


    t0 <- Sys.time()
    out <- parallel_lapply(
        ncore = ncore, config_list = config_list,
        fun = run_hyper, hyper_dmi = hyper_dmi, samples_list = samples_list
    )

    t1 <- Sys.time()

    proc_time <- difftime(t1, t0, units = "secs")[[1]]
    message("Processing time: ", round(proc_time, 2), " secs.")

    # Preserve attributes from original samples
    attr(out, "configs") <- config_list
    attr(out, "hyper_dmi") <- hyper_dmi
    attr(out, "dmis") <- dmis
    attr(out, "priors") <- priors

    return(out)
}
