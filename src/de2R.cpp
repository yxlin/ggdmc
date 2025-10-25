#include "de.h"
#include "type_casting.h"

/* ------------- Helper to extract timing stats ----------------------- */
Rcpp::List get_timing_stats(const de_class::TimingStats &stats)
{
    return Rcpp::List::create(
        Rcpp::Named("total_time") = stats.total_time,
        Rcpp::Named("crossover_time") = stats.crossover_time,
        Rcpp::Named("migration_time") = stats.migration_time,
        Rcpp::Named("likelihood_time") = stats.likelihood_time,
        Rcpp::Named("n_crossover") = stats.n_crossover,
        Rcpp::Named("n_migration") = stats.n_migration);
}

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

//' @rdname run
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

/* ------------- Fast versions with timing stats ---------------------- */

//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::List run_subject_fast(const Rcpp::S4 &config_r, const Rcpp::S4 &dmi,
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

    de_ptr->run_chains_fast(t_ptr, p_ptr, l_ptr, de_input.sub_debug);

    Rcpp::S4 posterior = new_posterior(t_ptr, th_input);
    Rcpp::List timing = get_timing_stats(de_ptr->m_timing_stats);

    return Rcpp::List::create(Rcpp::Named("posterior") = posterior,
                              Rcpp::Named("timing") = timing);
}

//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::List run_hyper_fast(const Rcpp::S4 &config_r, const Rcpp::S4 &dmi,
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

    de_ptr->run_chains_fast(t_ptr, h_ptr, l_ptr, de_input.sub_debug);

    Rcpp::S4 posterior = new_posterior(t_ptr, th_input);
    Rcpp::List timing = get_timing_stats(de_ptr->m_timing_stats);

    return Rcpp::List::create(Rcpp::Named("posterior") = posterior,
                              Rcpp::Named("timing") = timing);
}

/* ------------- Fast versions WITHOUT timing overhead --------------- */

//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::S4 run_subject_fast_notiming(const Rcpp::S4 &config_r,
                                   const Rcpp::S4 &dmi, const Rcpp::S4 &samples)
{
    Rcpp::S4 priors = config_r.slot("prior");
    auto p_ptr = prior::new_prior(priors.slot("p_prior"));
    auto l_ptr = new_likelihood(dmi);

    auto de_input = new_DEInput(config_r.slot("de_input"));
    auto th_input = new_ThetaInput(config_r.slot("theta_input"));
    auto theta = new_SampleInput(samples);
    auto t_ptr = std::make_shared<theta_class>(th_input, theta);
    auto de_ptr = std::make_shared<de_class>(de_input);

    de_ptr->run_chains_fast_notiming(t_ptr, p_ptr, l_ptr, de_input.sub_debug);

    return new_posterior(t_ptr, th_input);
}

//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::S4 run_hyper_fast_notiming(const Rcpp::S4 &config_r, const Rcpp::S4 &dmi,
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

    de_ptr->run_chains_fast_notiming(t_ptr, h_ptr, l_ptr, de_input.sub_debug);

    return new_posterior(t_ptr, th_input);
}

/* ------------- HB Fast versions with cached likelihood ------------- */

//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::List run_fast(const Rcpp::S4 &config_r, const Rcpp::List &dmis,
                    const Rcpp::List &samples)
{
    Rcpp::S4 priors = config_r.slot("prior");
    auto hyper_likelihood = prior::new_prior(priors.slot("p_prior"));
    auto hyper_prior = prior::new_prior(priors.slot("h_prior"));
    auto l_ptrs = new_likelihoods(dmis);

    auto de_input = new_DEInput(config_r.slot("de_input"));
    auto phi_input = new_ThetaInput(config_r.slot("theta_input"));

    ThetaInput th_input = phi_input;
    th_input.pnames = l_ptrs[0]->m_model->m_free_parameter_names;
    th_input.nparameter = l_ptrs[0]->m_model->m_n_free_parameter;

    Rcpp::List subj_theta_r = samples["subject_theta"];
    unsigned int n_subject = l_ptrs.size();

    std::vector<std::shared_ptr<theta_class>> subj_theta(n_subject);
    for (size_t i = 0; i < n_subject; ++i)
    {
        auto theta = new_SampleInput(subj_theta_r[i]);
        subj_theta[i] = std::make_shared<theta_class>(th_input, theta);
    }

    auto phi = new_SampleInput(samples["phi"]);
    auto phi_ptr = std::make_shared<theta_class>(phi_input, phi);
    auto de_ptr = std::make_shared<de_class>(de_input, n_subject);

    de_ptr->run_hchains_fast(phi_ptr, l_ptrs, subj_theta, hyper_likelihood,
                             hyper_prior, de_input.pop_debug,
                             de_input.sub_debug);

    Rcpp::List theta_out(n_subject);
    for (size_t subj_idx = 0; subj_idx < n_subject; ++subj_idx)
    {
        theta_out[subj_idx] = new_posterior(subj_theta[subj_idx], th_input);
    }

    Rcpp::S4 phi_out = new_posterior(phi_ptr, phi_input);
    Rcpp::List timing = get_timing_stats(de_ptr->m_timing_stats);

    Rcpp::List out = Rcpp::List::create(
        Rcpp::Named("phi") = phi_out, Rcpp::Named("subject_theta") = theta_out,
        Rcpp::Named("timing") = timing);

    return out;
}

//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::List run_fast_notiming(const Rcpp::S4 &config_r, const Rcpp::List &dmis,
                             const Rcpp::List &samples)
{
    Rcpp::S4 priors = config_r.slot("prior");
    auto hyper_likelihood = prior::new_prior(priors.slot("p_prior"));
    auto hyper_prior = prior::new_prior(priors.slot("h_prior"));
    auto l_ptrs = new_likelihoods(dmis);

    auto de_input = new_DEInput(config_r.slot("de_input"));
    auto phi_input = new_ThetaInput(config_r.slot("theta_input"));

    ThetaInput th_input = phi_input;
    th_input.pnames = l_ptrs[0]->m_model->m_free_parameter_names;
    th_input.nparameter = l_ptrs[0]->m_model->m_n_free_parameter;

    Rcpp::List subj_theta_r = samples["subject_theta"];
    unsigned int n_subject = l_ptrs.size();

    std::vector<std::shared_ptr<theta_class>> subj_theta(n_subject);
    for (size_t i = 0; i < n_subject; ++i)
    {
        auto theta = new_SampleInput(subj_theta_r[i]);
        subj_theta[i] = std::make_shared<theta_class>(th_input, theta);
    }

    auto phi = new_SampleInput(samples["phi"]);
    auto phi_ptr = std::make_shared<theta_class>(phi_input, phi);
    auto de_ptr = std::make_shared<de_class>(de_input, n_subject);

    de_ptr->run_hchains_fast_notiming(phi_ptr, l_ptrs, subj_theta,
                                      hyper_likelihood, hyper_prior,
                                      de_input.pop_debug, de_input.sub_debug);

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

/* ==================== Step 3: Type Conversion Optimization ==================
 * These wrappers use run_chains_typeconv() which calls sumloglike_fast()
 * Eliminates arma::vec -> std::vector<double> conversions in likelihood
 * evaluation
 * =================================================================================*/

//' @rdname run
//' @export
// [[Rcpp::export]]
Rcpp::S4 run_subject_typeconv(const Rcpp::S4 &config_r, const Rcpp::S4 &dmi,
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

    de_ptr->run_chains_typeconv(t_ptr, p_ptr, l_ptr, de_input.sub_debug);
    return new_posterior(t_ptr, th_input);
}
