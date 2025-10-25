// Optimized HB functions with cached hyper-likelihood
// This file contains the fast versions of hierarchical Bayesian sampling functions

#include "de.h"

/* ===========================================================================
   STRATEGY: Cache hyper-likelihood to avoid redundant calculations

   Problem: Original code recalculates sumloghlike() even when subject thetas
            haven't changed. This loops through ALL subjects unnecessarily.

   Solution:
   1. Cache hyper-likelihood after subject-level updates (when it changes)
   2. Reuse cached values at hyper-level (when it's still valid)
   3. Only recalculate when actually needed for M-H algorithm

   Impact: O(n_chains * n_subjects) → O(n_chains) per hyper-level step
   ========================================================================= */

// ============================================================================
// Version 1: WITH timing instrumentation (for profiling)
// ============================================================================

void de_class::crossover_hb_fast(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                                  PriorPtr hyper_likelihood, PriorPtr h_prior,
                                  bool debug, size_t para_idx)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    if (debug)
    {
        Rcpp::Rcout << "-----Start phi_crossover (FAST): "
                    << "Proposed lp, ll, new vs old posterior, & [MH ratio]\n";
    }

    for (size_t i = 0; i < m_nchain; ++i)
    {
        // OPTIMIZATION: Only calculate if cache is invalid (after subject updates)
        // Otherwise, reuse the cached value from last iteration
        if (!m_hlike_cache_valid)
        {
            auto cache_start = std::chrono::high_resolution_clock::now();
            m_cached_hlike(i) = sumloghlike(phi_ptr->m_used_theta.col(i), i,
                                            subj_thetas, hyper_likelihood);
            auto cache_end = std::chrono::high_resolution_clock::now();
            m_timing_stats.likelihood_time +=
                std::chrono::duration<double>(cache_end - cache_start).count();
        }

        phi_ptr->m_used_ll(i) = m_cached_hlike(i);
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
            for (size_t j = 0; j < m_nparameter; ++j)
            {
                m_tmp_theta(j) +=
                    Rf_runif(-m_rp, m_rp) + m_gamma * (theta0(j) - theta1(j));
            }
        }

        /* ---------------------- M-H algorithm----------------------*/
        m_tmp_lp = h_prior->sumlogprior(m_tmp_theta);

        auto ll_start = std::chrono::high_resolution_clock::now();
        m_tmp_ll = sumloghlike(m_tmp_theta, i, subj_thetas, hyper_likelihood);
        auto ll_end = std::chrono::high_resolution_clock::now();
        m_timing_stats.likelihood_time +=
            std::chrono::duration<double>(ll_end - ll_start).count();

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
            m_cached_hlike(i) = m_tmp_ll; // Update cache for accepted proposal
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

    // Cache is now valid for all chains after this crossover
    validate_hlike_cache();

    auto end_time = std::chrono::high_resolution_clock::now();
    m_timing_stats.crossover_time +=
        std::chrono::duration<double>(end_time - start_time).count();
    m_timing_stats.n_crossover++;
}

void de_class::migration_hb_fast(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                                  PriorPtr hyper_likelihood, PriorPtr h_prior,
                                  bool debug, size_t para_idx)
{
    auto start_time = std::chrono::high_resolution_clock::now();

    m_subchains = get_subchains();
    unsigned int next_chain;

    if (debug)
    {
        Rcpp::Rcout << "-----Start phi_migration (FAST): ";
        Rcpp::Rcout << "m_nchain: " << m_nchain << std::endl;
        Rcpp::Rcout << "m_subchains:" << m_subchains.t();
    }

    for (size_t i = 0; i < m_nsubchain; ++i)
    {
        next_chain =
            ((i + 1) == m_nsubchain) ? m_subchains(0) : m_subchains(i + 1);

        // OPTIMIZATION: Only calculate if cache is invalid
        // This eliminates 2 redundant sumloghlike() calls per migration!
        if (!m_hlike_cache_valid)
        {
            auto cache_start = std::chrono::high_resolution_clock::now();
            m_cached_hlike(m_subchains(i)) =
                sumloghlike(phi_ptr->m_used_theta.col(m_subchains(i)),
                           m_subchains(i), subj_thetas, hyper_likelihood);
            m_cached_hlike(next_chain) =
                sumloghlike(phi_ptr->m_used_theta.col(next_chain),
                           next_chain, subj_thetas, hyper_likelihood);
            auto cache_end = std::chrono::high_resolution_clock::now();
            m_timing_stats.likelihood_time +=
                std::chrono::duration<double>(cache_end - cache_start).count();
        }

        phi_ptr->m_used_ll(m_subchains(i)) = m_cached_hlike(m_subchains(i));
        phi_ptr->m_used_ll(next_chain) = m_cached_hlike(next_chain);

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

        auto ll_start = std::chrono::high_resolution_clock::now();
        m_tmp_ll = sumloghlike(m_tmp_theta, m_subchains(i), subj_thetas,
                               hyper_likelihood);
        auto ll_end = std::chrono::high_resolution_clock::now();
        m_timing_stats.likelihood_time +=
            std::chrono::duration<double>(ll_end - ll_start).count();

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
            m_cached_hlike(next_chain) = m_tmp_ll; // Update cache
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

    validate_hlike_cache();

    if (debug)
    {
        Rcpp::Rcout << "-----End phi_migration (FAST)\n\n";
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    m_timing_stats.migration_time +=
        std::chrono::duration<double>(end_time - start_time).count();
    m_timing_stats.n_migration++;
}

void de_class::run_hchains_fast(ThetaPtr phi_ptr, std::vector<LPtr> l_ptrs,
                                std::vector<ThetaPtr> subj_thetas,
                                PriorPtr hyper_likelihood, PriorPtr h_prior,
                                bool pop_debug, bool sub_debug)
{
    // Reset timing statistics
    m_timing_stats = TimingStats();
    auto total_start = std::chrono::high_resolution_clock::now();

    if (pop_debug)
    {
        Rcpp::Rcout << "run_hchains_fast (with caching)\n";
    }

    for (size_t MC_idx = 1; MC_idx < phi_ptr->m_nsample; ++MC_idx)
    {
        /*---------------------Hyper level--------------------*/
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
                    migration_hb_fast(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                                     pop_debug, para_idx);
                }
                else
                {
                    crossover_hb_fast(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                                     pop_debug, para_idx);
                }
            }
        }
        else
        {
            if (Rf_runif(0.0, 1.0) < m_pop_migration_prob)
            {
                migration_hb_fast(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                                 pop_debug);
            }
            else
            {
                crossover_hb_fast(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                                 pop_debug);
            }
        }

        /*---------------------Subject level------------------*/
        // Invalidate cache before subject updates
        invalidate_hlike_cache();

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

        // Cache will be recalculated on next hyper-level update
        if (pop_debug)
        {
            Rcpp::Rcout << std::endl;
        }
        phi_ptr->store(MC_idx);
        phi_ptr->print_progress(MC_idx);
    }
    Rcpp::Rcout << std::endl;

    auto total_end = std::chrono::high_resolution_clock::now();
    m_timing_stats.total_time =
        std::chrono::duration<double>(total_end - total_start).count();
}

// ============================================================================
// Version 2: WITHOUT timing instrumentation (pure performance)
// ============================================================================

void de_class::crossover_hb_fast_notiming(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                                          PriorPtr hyper_likelihood, PriorPtr h_prior,
                                          bool debug, size_t para_idx)
{
    if (debug)
    {
        Rcpp::Rcout << "-----Start phi_crossover (FAST): "
                    << "Proposed lp, ll, new vs old posterior, & [MH ratio]\n";
    }

    for (size_t i = 0; i < m_nchain; ++i)
    {
        // OPTIMIZATION: Reuse cached value if valid
        if (!m_hlike_cache_valid)
        {
            m_cached_hlike(i) = sumloghlike(phi_ptr->m_used_theta.col(i), i,
                                            subj_thetas, hyper_likelihood);
        }

        phi_ptr->m_used_ll(i) = m_cached_hlike(i);
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
            for (size_t j = 0; j < m_nparameter; ++j)
            {
                m_tmp_theta(j) +=
                    Rf_runif(-m_rp, m_rp) + m_gamma * (theta0(j) - theta1(j));
            }
        }

        /* ---------------------- M-H algorithm----------------------*/
        m_tmp_lp = h_prior->sumlogprior(m_tmp_theta);
        m_tmp_ll = sumloghlike(m_tmp_theta, i, subj_thetas, hyper_likelihood);
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
            m_cached_hlike(i) = m_tmp_ll;
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

    validate_hlike_cache();
}

void de_class::migration_hb_fast_notiming(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                                          PriorPtr hyper_likelihood, PriorPtr h_prior,
                                          bool debug, size_t para_idx)
{
    m_subchains = get_subchains();
    unsigned int next_chain;

    if (debug)
    {
        Rcpp::Rcout << "-----Start phi_migration (FAST): ";
        Rcpp::Rcout << "m_nchain: " << m_nchain << std::endl;
        Rcpp::Rcout << "m_subchains:" << m_subchains.t();
    }

    for (size_t i = 0; i < m_nsubchain; ++i)
    {
        next_chain =
            ((i + 1) == m_nsubchain) ? m_subchains(0) : m_subchains(i + 1);

        // OPTIMIZATION: Reuse cached values
        if (!m_hlike_cache_valid)
        {
            m_cached_hlike(m_subchains(i)) =
                sumloghlike(phi_ptr->m_used_theta.col(m_subchains(i)),
                           m_subchains(i), subj_thetas, hyper_likelihood);
            m_cached_hlike(next_chain) =
                sumloghlike(phi_ptr->m_used_theta.col(next_chain),
                           next_chain, subj_thetas, hyper_likelihood);
        }

        phi_ptr->m_used_ll(m_subchains(i)) = m_cached_hlike(m_subchains(i));
        phi_ptr->m_used_ll(next_chain) = m_cached_hlike(next_chain);

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
            m_cached_hlike(next_chain) = m_tmp_ll;
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

    validate_hlike_cache();

    if (debug)
    {
        Rcpp::Rcout << "-----End phi_migration (FAST)\n\n";
    }
}

void de_class::run_hchains_fast_notiming(ThetaPtr phi_ptr, std::vector<LPtr> l_ptrs,
                                        std::vector<ThetaPtr> subj_thetas,
                                        PriorPtr hyper_likelihood, PriorPtr h_prior,
                                        bool pop_debug, bool sub_debug)
{
    if (pop_debug)
    {
        Rcpp::Rcout << "run_hchains_fast_notiming (with caching)\n";
    }

    for (size_t MC_idx = 1; MC_idx < phi_ptr->m_nsample; ++MC_idx)
    {
        /*---------------------Hyper level--------------------*/
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
                    migration_hb_fast_notiming(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                                              pop_debug, para_idx);
                }
                else
                {
                    crossover_hb_fast_notiming(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                                              pop_debug, para_idx);
                }
            }
        }
        else
        {
            if (Rf_runif(0.0, 1.0) < m_pop_migration_prob)
            {
                migration_hb_fast_notiming(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                                          pop_debug);
            }
            else
            {
                crossover_hb_fast_notiming(phi_ptr, subj_thetas, hyper_likelihood, h_prior,
                                          pop_debug);
            }
        }

        /*---------------------Subject level------------------*/
        invalidate_hlike_cache();

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
