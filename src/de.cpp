#include "de.h"

/* ----------------Constructor----------------*/
de_class::de_class(const DEInput &de_input)
{
    yank_de_input(de_input);
    if (m_nchain <= m_min_nsubchain)
    {
        throw std::runtime_error("Require three or more chains.");
    }
    /*--------tuning parameters------------------------------*/
    m_gamma = m_gamma_precursor / std::sqrt(2.0 * m_nparameter);
    m_chains = arma::linspace<arma::uvec>(0, m_nchain - 1, m_nchain);
}

de_class::de_class(const DEInput &de_input, unsigned int nsubject)
{
    yank_de_input(de_input);
    if (m_nchain <= m_min_nsubchain)
    {
        throw std::runtime_error("Require three or more chains.");
    }
    /*--------tuning parameters------------------------------*/
    m_gamma = m_gamma_precursor / std::sqrt(2.0 * m_nparameter);
    m_chains = arma::linspace<arma::uvec>(0, m_nchain - 1, m_nchain);

    m_nsubject = nsubject;
}

de_class::~de_class()
{
}

/* -------------------Private functions -------------------------*/
void de_class::yank_de_input(const DEInput &de_input)
{
    m_pop_migration_prob = de_input.pop_migration_prob;
    m_sub_migration_prob = de_input.sub_migration_prob;

    m_gamma_precursor = de_input.gamma_precursor;
    m_rp = de_input.rp;

    m_is_hblocked = de_input.is_hblocked;
    m_is_pblocked = de_input.is_pblocked;

    m_nparameter = de_input.nparameter;
    m_half_nparameter = m_nparameter / 2;
    m_nchain = de_input.nchain;
}

/* -------------------Generic functions -------------------------*/
// Exclude the kth chain, shuffle the remaining, and then select nsubchain
// chains
arma::uvec de_class::get_chains(unsigned int k, unsigned int nsubchain)
{
    arma::uvec chains_copy = m_chains;
    chains_copy.shed_row(k);
    arma::uvec rchains = arma::shuffle(chains_copy);
    return rchains.rows(0, nsubchain - 1);
}

arma::uvec de_class::get_subchains()
{
    // Calculate the number of subchains to pick by generating random proportion
    double proportion = Rf_runif(0.0, 1.0);
    m_nsubchain = static_cast<unsigned int>(std::ceil(m_nchain * proportion));

    // Ensure nsubchain is at least min_nsubchain and at most m_nchain
    m_nsubchain = std::max(m_nsubchain, m_min_nsubchain);
    m_nsubchain = std::min(m_nsubchain, m_nchain);

    // Then shuffle all chains.
    arma::uvec rchains = arma::shuffle(m_chains);
    // Select the first nsubchain chains (no need to sort unless required)
    arma::uvec out = rchains.head(m_nsubchain);
    return arma::sort(out);
    // return rchains.rows(0, m_nsubchain - 1);
}

/* -------------------De Samplers -------------------------*/
void de_class::update_theta(ThetaPtr t_ptr, size_t chain_idx, bool debug)
{
    if (std::isnan(m_mh_ratio))
    {
        if (debug)
            Rcpp::Rcout << "rejected (mh nan)\n";
    }
    else if (Rf_runif(0, 1) < m_mh_ratio)
    {
        if (debug)
        {
            Rcpp::Rcout << "  [C" << chain_idx << "]: ";
            Rcpp::Rcout << m_tmp_lp << ", " << m_tmp_ll << ", ("
                        << m_tmp_log_pos << " vs " << m_cur_log_pos << "), \t["
                        << m_mh_ratio << "]-> ";
            Rcpp::Rcout << "accepted\n";
        }

        t_ptr->m_used_theta.col(chain_idx) = m_tmp_theta;
        t_ptr->m_used_lp(chain_idx) = m_tmp_lp;
        t_ptr->m_used_ll(chain_idx) = m_tmp_ll;
    }
    else
    {
        if (debug)
            Rcpp::Rcout << "rejected\n";
    }
}

/* -------------Theta and hyper only Samplers -------------*/
void de_class::crossover(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr, bool debug,
                         size_t para_idx)
{
    if (debug)
    {
        Rcpp::Rcout << "-----Start theta crossover step: "
                    << "Proposed log prior, like, & post [MH ratio]\n";
    }
    for (size_t i = 0; i < m_nchain; ++i)
    {
        m_cur_log_pos = t_ptr->m_used_ll(i) + t_ptr->m_used_lp(i);
        m_tmp_theta = t_ptr->m_used_theta.col(i);
        m_subchains = get_chains(i, 2);
        const arma::vec &theta0 = t_ptr->m_used_theta.col(m_subchains(0));
        const arma::vec &theta1 = t_ptr->m_used_theta.col(m_subchains(1));

        /*-----------------Blocking switch----------------*/
        if (para_idx != SIZE_T_MAX)
        {
            m_tmp_theta(para_idx) +=
                Rf_runif(-m_rp, m_rp) +
                m_gamma * (theta0(para_idx) - theta1(para_idx));
        }
        else
        {
            for (size_t j = 0; j < m_nparameter; ++j)
            {
                m_tmp_theta(j) +=
                    Rf_runif(-m_rp, m_rp) + m_gamma * (theta0(j) - theta1(j));
            }
        }

        /* ---------------------- M-H algorithm----------------------*/
        m_tmp_lp = p_ptr->sumlogprior(m_tmp_theta);
        m_tmp_ll = l_ptr->sumloglike(m_tmp_theta, debug);
        m_tmp_log_pos = m_tmp_lp + m_tmp_ll;
        m_mh_ratio = std::exp(m_tmp_log_pos - m_cur_log_pos);

        update_theta(t_ptr, i, debug);
    }
    if (debug)
    {
        Rcpp::Rcout << "-----End theta crossover\n\n";
    }
}

void de_class::migration(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr, bool debug,
                         size_t para_idx)
{
    if (debug)
    {
        Rcpp::Rcout << "-----Start theta migration step: ";
    }
    m_subchains = get_subchains(); // eg, 0, 1, 3, 4, 8;
    unsigned int next_chain;

    for (size_t i = 0; i < m_nsubchain; i++)
    {
        next_chain =
            ((i + 1) == m_nsubchain) ? m_subchains(0) : m_subchains(i + 1);
        m_tmp_theta = t_ptr->m_used_theta.col(m_subchains(i));

        /*-----------------Blocking switch----------------*/
        if (para_idx != SIZE_T_MAX)
        {
            m_tmp_theta(para_idx) += Rf_runif(-m_rp, m_rp);
        }
        else
        {
            for (size_t j = 0; j < m_nparameter; ++j)
            {
                m_tmp_theta(j) += Rf_runif(-m_rp, m_rp);
            }
        }

        m_tmp_lp = p_ptr->sumlogprior(m_tmp_theta);
        m_tmp_ll = l_ptr->sumloglike(m_tmp_theta);
        m_tmp_log_pos = m_tmp_lp + m_tmp_ll;
        m_cur_log_pos =
            t_ptr->m_used_lp[next_chain] + t_ptr->m_used_ll[next_chain];
        m_mh_ratio = std::exp(m_tmp_log_pos - m_cur_log_pos);

        update_theta(t_ptr, next_chain, debug);
    }
    if (debug)
    {
        Rcpp::Rcout << "-----End theta migration_Hu\n\n";
    }
}

void de_class::run_chains(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                          bool debug)
{
    std::vector<std::pair<double, std::function<void()>>> migration_actions = {
        {m_sub_migration_prob,
         [&]() { migration(t_ptr, p_ptr, l_ptr, debug); }}};

    for (size_t i = 1; i < t_ptr->m_nsample; ++i)
    {
        double rand_val = Rf_runif(0.0, 1.0);
        bool action_taken = false;

        for (const auto &pair : migration_actions)
        {
            if (rand_val < pair.first)
            {
                pair.second();
                action_taken = true;
                break;
            }
        }
        if (!action_taken)
        {
            crossover(t_ptr, p_ptr, l_ptr, debug);
        }

        t_ptr->store(i);
        t_ptr->print_progress(i);
    }
    Rcpp::Rcout << std::endl;
}

/* -------------------HB Path -------------------------*/
double de_class::sumloghlike(const arma::vec &parameters, size_t chain_idx,
                             std::vector<ThetaPtr> subj_thetas,
                             PriorPtr p_prior)
{

    p_prior->m_p0 = parameters.head(m_half_nparameter);
    p_prior->m_p1 = parameters.tail(m_half_nparameter);
    // Rcpp::Rcout << "parameters " << parameters.t() << std::endl;
    // Rcpp::Rcout << "head " << parameters.head(m_half_nparameter).t()
    //             << std::endl;
    // Rcpp::Rcout << "tail " << parameters.tail(m_half_nparameter).t()
    //             << std::endl;

    // Rcpp::Rcout << "m half nparameter " << m_half_nparameter << "  ";
    // Rcpp::Rcout << "first subj theta: "
    //             << subj_thetas[0]->m_used_theta.col(chain_idx).t() <<
    //             std::endl;

    double out = 0;
    for (size_t j = 0; j < m_nsubject; ++j)
    {
        out +=
            p_prior->sumlogprior(subj_thetas[j]->m_used_theta.col(chain_idx));
    }
    return out;
}

void de_class::run_hchains(ThetaPtr phi_ptr, std::vector<LPtr> l_ptrs,
                           std::vector<ThetaPtr> subj_thetas,
                           PriorPtr hyper_likelihood, PriorPtr h_prior,
                           bool pop_debug, bool sub_debug)
{
    if (pop_debug)
    {
        Rcpp::Rcout << "new run hchain\n";
    }
    for (size_t MC_idx = 1; MC_idx < phi_ptr->m_nsample; ++MC_idx)
    {
        if (m_is_hblocked)
        {
            if (pop_debug)
            {
                Rcpp::Rcout << "Parameter: ";
            }
            for (size_t para_idx = 0; para_idx < m_nparameter; ++para_idx)
            {
                if (pop_debug)
                {
                    std::string pname = phi_ptr->m_pnames[para_idx];
                    Rcpp::Rcout << pname << "\n";
                }
                if (Rf_runif(0.0, 1.0) < m_pop_migration_prob)
                {
                    migration(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                              pop_debug, para_idx);
                }
                else
                {
                    crossover(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                              pop_debug, para_idx);
                }
            }
        }
        else
        {
            if (Rf_runif(0.0, 1.0) < m_pop_migration_prob)
            {
                migration(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                          pop_debug);
            }
            else
            {
                crossover(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                          pop_debug);
            }
        }
        /*---------------------Subject level------------------*/

        for (size_t subj_idx = 0; subj_idx < m_nsubject; ++subj_idx)
        {
            auto l_ptr = l_ptrs[subj_idx];
            auto t_ptr = subj_thetas[subj_idx];

            bool is_subject_debug = sub_debug && subj_idx <= 2;
            if (is_subject_debug)
            {
                Rcpp::Rcout << "Subject " << subj_idx << ", ";
            }
            if (m_is_pblocked)
            {
                if (is_subject_debug)
                {
                    Rcpp::Rcout << "sub-parameter: ";
                }
                for (size_t para_idx = 0; para_idx < m_half_nparameter;
                     ++para_idx)
                {
                    if (is_subject_debug)
                    {
                        std::string pname =
                            subj_thetas[subj_idx]->m_pnames[para_idx];
                        Rcpp::Rcout << pname << "\n";
                    }
                    if (Rf_runif(0.0, 1.0) < m_sub_migration_prob)
                    {
                        migration(phi_ptr, l_ptr, t_ptr, hyper_likelihood,
                                  is_subject_debug, para_idx);
                    }
                    else
                    {
                        crossover(phi_ptr, l_ptr, t_ptr, hyper_likelihood,
                                  is_subject_debug, para_idx);
                    }
                }
            }
            else
            {
                if (Rf_runif(0.0, 1.0) < m_sub_migration_prob)
                {
                    migration(phi_ptr, l_ptr, t_ptr, hyper_likelihood,
                              is_subject_debug);
                }
                else
                {
                    crossover(phi_ptr, l_ptr, t_ptr, hyper_likelihood,
                              is_subject_debug);
                }
            }
            subj_thetas[subj_idx]->store(MC_idx);
        }
        if (pop_debug)
        {
            Rcpp::Rcout << std::endl;
        }
        phi_ptr->store(MC_idx);
        phi_ptr->print_progress(MC_idx);
    }
    Rcpp::Rcout << std::endl;
}
// Hyper level
void de_class::crossover(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                         PriorPtr hyper_likelihood, PriorPtr h_prior,
                         bool debug, size_t para_idx)
{
    if (debug)
    {
        Rcpp::Rcout << "-----Start phi_crossover: "
                    << "Proposed lp, ll, new vs old posterior, & [MH ratio]\n";
    }
    for (size_t i = 0; i < m_nchain; ++i)
    {
        /* Update h like, if subject_theta is changed */
        phi_ptr->m_used_ll(i) = sumloghlike(phi_ptr->m_used_theta.col(i), i,
                                            subj_thetas, hyper_likelihood);
        // Rcpp::Rcout << "Finish first sumloghlike\n";
        m_cur_log_pos = phi_ptr->m_used_ll(i) + phi_ptr->m_used_lp(i);
        m_tmp_theta = phi_ptr->m_used_theta.col(i);
        m_subchains = get_chains(i, 2);

        const arma::vec &theta0 = phi_ptr->m_used_theta.col(m_subchains(0));
        const arma::vec &theta1 = phi_ptr->m_used_theta.col(m_subchains(1));

        /*-----------------Blocking switch----------------*/
        if (para_idx != SIZE_T_MAX)
        {
            m_tmp_theta(para_idx) +=
                Rf_runif(-m_rp, m_rp) +
                m_gamma * (theta0(para_idx) - theta1(para_idx));
        }
        else
        {
            // Rcpp::Rcout << "m_nparameter = " << m_nparameter << "\n";
            // Rcpp::Rcout << "m_tmp_theta" << m_tmp_theta.t() << "\n";
            for (size_t j = 0; j < m_nparameter; ++j)
            {
                m_tmp_theta(j) +=
                    Rf_runif(-m_rp, m_rp) + m_gamma * (theta0(j) - theta1(j));
            }
        }

        /* ---------------------- M-H algorithm----------------------*/
        m_tmp_lp = h_prior->sumlogprior(m_tmp_theta);
        m_tmp_ll = sumloghlike(m_tmp_theta, i, subj_thetas, hyper_likelihood);
        // Rcpp::Rcout << "Finish second sumloghlike\n";

        m_tmp_log_pos = m_tmp_lp + m_tmp_ll;
        m_mh_ratio = std::exp(m_tmp_log_pos - m_cur_log_pos);

        if (std::isnan(m_mh_ratio))
        {
            if (debug)
            {
                Rcpp::Rcout << "  [C" << i << "]: ";
                Rcpp::Rcout << "rejected (mh nan)\n";
            }
        }
        else if (Rf_runif(0, 1) < m_mh_ratio)
        {
            if (debug)
            {
                Rcpp::Rcout << "  [C" << i << "]: ";
                Rcpp::Rcout << m_tmp_lp << ", " << m_tmp_ll << ", ("
                            << m_tmp_log_pos << " vs " << m_cur_log_pos
                            << "), \t[" << m_mh_ratio << "]-> ";
                Rcpp::Rcout << "accepted\n";
            }

            phi_ptr->m_used_theta.col(i) = m_tmp_theta;
            phi_ptr->m_used_lp(i) = m_tmp_lp;
            phi_ptr->m_used_ll(i) = m_tmp_ll;
        }
        else
        {
            if (debug)
            {
                Rcpp::Rcout << "  [C" << i << "]: ";
                Rcpp::Rcout << "rejected\n";
            }
        }
    }
}

void de_class::migration(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                         PriorPtr hyper_likelihood, PriorPtr h_prior,
                         bool debug, size_t para_idx)
{
    m_subchains = get_subchains(); // eg, 0, 1, 3, 4, 8;
    unsigned int next_chain;

    if (debug)
    {
        Rcpp::Rcout << "-----Start phi_migration: ";
        Rcpp::Rcout << "m_nchain: " << m_nchain << std::endl;
        Rcpp::Rcout << "m_subchains:" << m_subchains.t();
        Rcpp::Rcout << "phi_ptr->used_theta cols x rows: "
                    << phi_ptr->m_used_theta.n_cols << " x "
                    << phi_ptr->m_used_theta.n_rows << std::endl;
        if (para_idx != SIZE_T_MAX)
        {
            Rcpp::Rcout << "para_idx: " << para_idx << std::endl;
        }
        Rcpp::Rcout << std::endl;
    }
    for (size_t i = 0; i < m_nsubchain; ++i)
    {
        next_chain =
            ((i + 1) == m_nsubchain) ? m_subchains(0) : m_subchains(i + 1);

        /* Update h like, if subject_theta is changed */
        phi_ptr->m_used_ll(m_subchains(i)) =
            sumloghlike(phi_ptr->m_used_theta.col(m_subchains(i)),
                        m_subchains(i), subj_thetas, hyper_likelihood);

        phi_ptr->m_used_ll(next_chain) =
            sumloghlike(phi_ptr->m_used_theta.col(next_chain), next_chain,
                        subj_thetas, hyper_likelihood);

        m_tmp_theta = phi_ptr->m_used_theta.col(m_subchains(i));

        /*-----------------Blocking switch----------------*/
        if (para_idx != SIZE_T_MAX)
        {
            m_tmp_theta(para_idx) += Rf_runif(-m_rp, m_rp);
        }
        else
        {
            for (size_t j = 0; j < m_nparameter; ++j)
            {
                m_tmp_theta(j) += Rf_runif(-m_rp, m_rp);
            }
        }

        /* ---------------------- M-H algorithm----------------------*/
        m_tmp_lp = h_prior->sumlogprior(m_tmp_theta);
        m_tmp_ll = sumloghlike(m_tmp_theta, m_subchains(i), subj_thetas,
                               hyper_likelihood);
        m_tmp_log_pos = m_tmp_lp + m_tmp_ll;

        m_cur_log_pos =
            phi_ptr->m_used_lp(next_chain) + phi_ptr->m_used_ll(next_chain);
        m_mh_ratio = std::exp(m_tmp_log_pos - m_cur_log_pos);

        if (std::isnan(m_mh_ratio))
        {
            if (debug)
            {
                Rcpp::Rcout << "  [C" << m_subchains(i) << "] vs [C"
                            << next_chain << "]: ";
                Rcpp::Rcout << "rejected (mh nan)\n";
            }
        }
        else if (Rf_runif(0, 1) < m_mh_ratio)
        {
            if (debug)
            {
                Rcpp::Rcout << "  [C" << m_subchains(i) << "] vs [C"
                            << next_chain << "]: ";
                Rcpp::Rcout << m_tmp_log_pos << " vs " << m_cur_log_pos
                            << "), \t[" << m_mh_ratio << "]-> ";
                Rcpp::Rcout << "accepted\n";
            }

            phi_ptr->m_used_theta.col(next_chain) = m_tmp_theta;
            phi_ptr->m_used_lp(next_chain) = m_tmp_lp;
            phi_ptr->m_used_ll(next_chain) = m_tmp_ll;
        }
        else
        {
            if (debug)
            {
                Rcpp::Rcout << "  [C" << m_subchains(i) << "] vs [C"
                            << next_chain << "]: ";
                Rcpp::Rcout << "rejected\n";
            }
        }
    }
    if (debug)
    {
        Rcpp::Rcout << "-----End phi_migration\n\n";
    }
}
// subject level
void de_class::crossover(ThetaPtr phi_ptr, LPtr l_ptr, ThetaPtr t_ptr,
                         PriorPtr p_prior, bool debug, size_t para_idx)
{
    if (debug)
    {
        Rcpp::Rcout << "-----Start theta-phi crossover: \n";
    }

    for (size_t i = 0; i < m_nchain; ++i)
    {
        m_cur_log_pos = t_ptr->m_used_ll(i) + t_ptr->m_used_lp(i);
        m_tmp_theta = t_ptr->m_used_theta.col(i);
        m_subchains = get_chains(i, 2);
        const arma::vec &theta0 = t_ptr->m_used_theta.col(m_subchains(0));
        const arma::vec &theta1 = t_ptr->m_used_theta.col(m_subchains(1));

        /*-----------------Blocking switch----------------*/
        if (para_idx != SIZE_T_MAX)
        {
            m_tmp_theta(para_idx) +=
                Rf_runif(-m_rp, m_rp) +
                m_gamma * (theta0(para_idx) - theta1(para_idx));
        }
        else
        {
            for (size_t j = 0; j < m_half_nparameter; ++j)
            {
                m_tmp_theta(j) +=
                    Rf_runif(-m_rp, m_rp) + m_gamma * (theta0(j) - theta1(j));
            }
        }
        /* ---------------------- Update p_prior----------------------*/
        p_prior->m_p0 = phi_ptr->m_used_theta.col(i).head(m_half_nparameter);
        p_prior->m_p1 = phi_ptr->m_used_theta.col(i).tail(m_half_nparameter);

        /* ---------------------- M-H algorithm----------------------*/
        m_tmp_lp = p_prior->sumlogprior(m_tmp_theta);
        m_tmp_ll = l_ptr->sumloglike(m_tmp_theta);
        m_tmp_log_pos = m_tmp_lp + m_tmp_ll;
        m_mh_ratio = std::exp(m_tmp_log_pos - m_cur_log_pos);
        update_theta(t_ptr, i, debug);
    }
    if (debug)
    {
        Rcpp::Rcout << "-----End theta-phi crossover\n\n";
    }
}

void de_class::migration(ThetaPtr phi_ptr, LPtr l_ptr, ThetaPtr t_ptr,
                         PriorPtr p_prior, bool debug, size_t para_idx)
{
    if (debug)
    {
        Rcpp::Rcout << "-----Start theta-phi migration: \n";
    }

    // auto t_ptr = phi_ptr->m_subj_theta[subj_idx];
    m_subchains = get_subchains(); // eg, 0, 1, 3, 4, 8;
    unsigned int next_chain;

    for (size_t i = 0; i < m_nsubchain; ++i)
    {
        next_chain =
            ((i + 1) == m_nsubchain) ? m_subchains(0) : m_subchains(i + 1);
        m_tmp_theta = t_ptr->m_used_theta.col(m_subchains(i));

        /*-----------------Blocking switch----------------*/
        if (para_idx != SIZE_T_MAX)
        {
            m_tmp_theta(para_idx) += Rf_runif(-m_rp, m_rp);
        }
        else
        {
            for (size_t j = 0; j < m_half_nparameter; ++j)
            {
                m_tmp_theta(j) += Rf_runif(-m_rp, m_rp);
            }
        }
        /* ---------------------- Update p_prior----------------------*/
        p_prior->m_p0 =
            phi_ptr->m_used_theta.col(m_subchains(i)).head(m_half_nparameter);
        p_prior->m_p1 =
            phi_ptr->m_used_theta.col(m_subchains(i)).tail(m_half_nparameter);

        /*------------------First run of the M-H algorithm------------------*/
        m_tmp_lp = p_prior->sumlogprior(m_tmp_theta);
        m_tmp_ll = l_ptr->sumloglike(m_tmp_theta);
        m_tmp_log_pos = m_tmp_lp + m_tmp_ll;

        m_cur_log_pos =
            t_ptr->m_used_lp(next_chain) + t_ptr->m_used_ll(next_chain);
        m_mh_ratio = std::exp(m_tmp_log_pos - m_cur_log_pos);
        update_theta(t_ptr, next_chain, debug);
    }
    if (debug)
    {
        Rcpp::Rcout << "-----End theta-phi migration_Hu\n\n";
    }
}