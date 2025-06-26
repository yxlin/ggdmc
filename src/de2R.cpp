#include "de.h"
#include "type_casting.h"

/* ------------- Run subjects-------------------------- */
//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::S4 run_subject(const Rcpp::S4 &config_r, const Rcpp::S4 &dmi,
                     const Rcpp::S4 &samples)
{
    Rcpp::S4 priors = config_r.slot("prior");
    auto p_ptr = prior::new_prior(priors.slot("p_prior"));
    auto l_ptr = new_likelihood(dmi);

    auto de_input = new_DEInput(config_r.slot("de_input"));
    auto th_input = new_ThetaInput(config_r.slot("theta_input"));
    auto theta = new_SampleInput(samples);
    auto t_ptr = std::make_shared<theta_class>(th_input, theta);
    auto de_ptr = std::make_shared<de_class>(de_input);

    de_ptr->run_chains(t_ptr, p_ptr, l_ptr, de_input.sub_debug);
    return new_posterior(t_ptr, th_input);
}

/* ------------- Run hyper-------------------------- */

//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::S4 run_hyper(const Rcpp::S4 &config_r, const Rcpp::S4 &dmi,
                   const Rcpp::S4 &samples)
{
    Rcpp::S4 priors = config_r.slot("prior");
    auto p_ptr = prior::new_prior(priors.slot("p_prior"));
    auto h_ptr = prior::new_prior(priors.slot("h_prior"));
    auto l_ptr = new_likelihood(dmi, p_ptr); // hyper_dmi

    auto de_input = new_DEInput(config_r.slot("de_input"));
    auto th_input = new_ThetaInput(config_r.slot("theta_input"));
    auto theta = new_SampleInput(samples);

    auto t_ptr = std::make_shared<theta_class>(th_input, theta);
    auto de_ptr = std::make_shared<de_class>(de_input);

    de_ptr->run_chains(t_ptr, h_ptr, l_ptr, de_input.sub_debug);
    return new_posterior(t_ptr, th_input);
}

/* ------------- Run HB---------------------------------------------- */

//' Run MCMC Sampling for Cognitive Models
//'
//' Executes MCMC-based posterior sampling for LBA or decision-diffusion
//' models. The model can be structured as either one or two level(s),
//' including subject-only, group-level, and hierarchical
//' multi-subject models.
//'
//' @description
//' These functions manage the MCMC sampling process for Linear
//' Ballistic Accumulator (LBA) models via C++ backends. The sampling is
//' configured using an S4 `config` object, which carries sampler settings
//' and model specifications. The `dmi` (data-model instance) contains
//' trial-level data and likelihood structure, and `samples` represents an
//' initial or resumed parameter state.
//'
//' @param config_r An S4 object of class \code{config} specifying model
//' settings, priors, tuning parameters, etc.
//' @param dmi For \code{run_subject()} and \code{run_hyper()}, a single S4
//' object containing a data-model instance for one subject or one set of
//' pressumed true theta values (in a R matrix format). For \code{run()}, a
//' list of \code{dmi} objects (one per subject).
//' @param dmis For \code{run()}, a list of \code{dmi} objects (one per
//'             subject).
//' @param samples For \code{run_subject()} and \code{run_hyper()}, an
//' S4 \code{samples} object containing starting parameter values.
//' For \code{run()}, a list of such \code{samples} objects for multiple
//' subjects and a phi `samples` for the hyper-level.
//'
//' @return
//' \itemize{
//'   \item \code{run_subject()} returns an S4 \code{posterior} object
//' for a single subject.
//'   \item \code{run_hyper()} returns an S4 \code{posterior} object
//'          for hyperparameters (group-level parameters).
//'   \item \code{run()} returns a named \code{list} of \code{posterior}
//'           objects, one per subject.
//' }
//'
//' @details
//' These functions are typically called from higher-level model fitting
//' workflows. Internally, they coordinate log-likelihood evaluation,
//' parameter proposal, acceptance checks, and diagnostics.
//'
//' \code{run_subject()} and \code{run_hyper()} differ in that the former
//' operates on trial-level likelihoods for an individual subject, while the
//' latter works with summary statistics and priors for group-level parameters.
//'
//' \code{run()} enables full hierarchical modeling across multiple subjects
//' by looping through each subject's `dmi` and `samples`,
//' applying \code{run_subject()}.
//'
//' @examples
//' \dontrun{
//' # Setup (assuming valid config_r, dmi, and samples objects exist)
//' posterior_subject <- run_subject(config_r, dmi, samples)
//'
//' # Group-level inference
//' posterior_hyper <- run_hyper(config_r, dmi, samples)
//'
//' # Hierarchical model
//' posteriors <- run(config_r, dmis = list(dmi1, dmi2), samples = list(s1, s2))
//' }
//'
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List run(const Rcpp::S4 &config_r, const Rcpp::List &dmis,
               const Rcpp::List &samples)
{
    Rcpp::S4 priors = config_r.slot("prior");
    auto hyper_likelihood = prior::new_prior(priors.slot("p_prior"));
    auto hyper_prior = prior::new_prior(priors.slot("h_prior"));
    auto l_ptrs = new_likelihoods(dmis); // likelihood of the LBA/DDM etc

    auto de_input = new_DEInput(config_r.slot("de_input"));
    auto phi_input = new_ThetaInput(config_r.slot("theta_input"));

    ThetaInput th_input = phi_input;
    th_input.pnames = l_ptrs[0]->m_model->m_free_parameter_names;
    th_input.nparameter = l_ptrs[0]->m_model->m_n_free_parameter;

    Rcpp::List subj_theta_r = samples["subject_theta"];
    unsigned int n_subject = l_ptrs.size();

    /* Initialise the storage */
    std::vector<std::shared_ptr<theta_class>> subj_theta(n_subject);
    for (size_t i = 0; i < n_subject; ++i)
    {
        auto theta = new_SampleInput(subj_theta_r[i]);
        subj_theta[i] = std::make_shared<theta_class>(th_input, theta);
    }

    auto phi = new_SampleInput(samples["phi"]); // phi

    auto phi_ptr = std::make_shared<theta_class>(phi_input, phi);

    auto de_ptr = std::make_shared<de_class>(de_input, n_subject);

    de_ptr->run_hchains(phi_ptr, l_ptrs, subj_theta, hyper_likelihood,
                        hyper_prior, de_input.pop_debug, de_input.sub_debug);

    Rcpp::List theta_out(n_subject);
    for (size_t subj_idx = 0; subj_idx < n_subject; ++subj_idx)
    {

        theta_out[subj_idx] = new_posterior(subj_theta[subj_idx], th_input);
    }

    Rcpp::S4 phi_out = new_posterior(phi_ptr, phi_input);

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("phi") = phi_out, Rcpp::Named("subject_theta") = theta_out);

    return out;
}
