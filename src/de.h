#pragma once

#include <ggdmcHeaders/likelihood.h>
#include <ggdmcHeaders/theta.h>
#include <chrono>

constexpr auto SIZE_T_MAX = std::numeric_limits<size_t>::max();

class de_class
{
    /* -------------------------------------------------------------------
    Theta and new HB
    ----------------------------------------------------------------------*/
  private:
    double m_gamma, m_gamma_precursor, m_rp;
    bool m_is_hblocked, m_is_pblocked;
    unsigned int m_nsubject;

    /*-----------Private theta_class functions-------------*/
    void yank_de_input(const DEInput &de_input);
    void update_theta(ThetaPtr t_ptr, size_t chain_idx, bool debug = false);

    arma::mat m_prev_theta; // To track previous parameter values

    /*-----------Optimized: pre-allocated buffers-------------*/
    arma::uvec m_chains_buffer;     // Pre-allocated buffer for chain shuffling
    arma::uvec m_subchains_buffer;  // Pre-allocated buffer for subchain selection
    arma::uvec m_shuffle_indices;   // Pre-allocated shuffle indices

    /*-----------HB: Cache for hyper-likelihood---------------*/
    arma::mat m_cached_hlike;       // Cached hyper-likelihood for each chain
    bool m_hlike_cache_valid;       // Whether cache is valid
    void invalidate_hlike_cache() { m_hlike_cache_valid = false; }
    void validate_hlike_cache() { m_hlike_cache_valid = true; }

  public:
    double m_pop_migration_prob, m_sub_migration_prob;
    unsigned int m_nparameter, m_nchain, m_half_nparameter, m_nsubchain;
    const unsigned int m_min_nsubchain = 2;

    /*------------Computed members ------------*/
    arma::vec m_tmp_theta;
    arma::uvec m_chains, m_subchains;
    double m_tmp_lp, m_tmp_ll, m_tmp_log_pos, m_cur_log_pos, m_mh_ratio;

    /*------------Constructor ------------*/
    de_class(const DEInput &de_input);
    de_class(const DEInput &de_input, unsigned int nsubject);
    ~de_class();

    arma::uvec get_chains(unsigned int k, unsigned int nsubchain);
    arma::uvec get_subchains();

    /*-----------Optimized versions with pre-allocated buffers-------------*/
    void get_chains_fast(unsigned int k, unsigned int nsubchain);
    void get_subchains_fast();

    /* Subject level only */
    void crossover(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                   bool debug = false, size_t para_idx = SIZE_T_MAX);
    void migration(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                   bool debug = false, size_t para_idx = SIZE_T_MAX);

    void run_chains(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                    bool debug = false);

    /*-----------Optimized versions-------------*/
    void crossover_fast(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                        bool debug = false, size_t para_idx = SIZE_T_MAX);
    void migration_fast(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                        bool debug = false, size_t para_idx = SIZE_T_MAX);
    void run_chains_fast(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                         bool debug = false);

    /*-----------Optimized versions without timing overhead-----------*/
    void crossover_fast_notiming(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                                 bool debug = false, size_t para_idx = SIZE_T_MAX);
    void migration_fast_notiming(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                                 bool debug = false, size_t para_idx = SIZE_T_MAX);
    void run_chains_fast_notiming(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                                  bool debug = false);

    /*-----------Timing instrumentation-------------*/
    struct TimingStats
    {
        double total_time = 0.0;
        double crossover_time = 0.0;
        double migration_time = 0.0;
        double likelihood_time = 0.0;
        unsigned int n_crossover = 0;
        unsigned int n_migration = 0;
    };
    TimingStats m_timing_stats;

    /* HB Procedure */
    double sumloghlike(const arma::vec &parameters, size_t chain_idx,
                       std::vector<ThetaPtr> subj_thetas, PriorPtr p_prior);
    void run_hchains(ThetaPtr phi_ptr, std::vector<LPtr> l_ptrs,
                     std::vector<ThetaPtr> subj_thetas,
                     PriorPtr hyper_likelihood, PriorPtr h_prior,
                     bool pop_debug = false, bool sub_debug = false);

    // Hyper level
    void crossover(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                   PriorPtr hyper_likelihood, PriorPtr h_prior,
                   bool debug = false, size_t para_idx = SIZE_T_MAX);
    void migration(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                   PriorPtr hyper_likelihood, PriorPtr h_prior,
                   bool debug = false, size_t para_idx = SIZE_T_MAX);
    // subject level
    void crossover(ThetaPtr phi_ptr, LPtr l_ptr, ThetaPtr t_ptr,
                   PriorPtr p_prior, bool debug = false,
                   size_t para_idx = SIZE_T_MAX);
    void migration(ThetaPtr phi_ptr, LPtr l_ptr, ThetaPtr t_ptr,
                   PriorPtr p_prior, bool debug = false,
                   size_t para_idx = SIZE_T_MAX);

    /* HB Optimized with cached likelihood */
    void run_hchains_fast(ThetaPtr phi_ptr, std::vector<LPtr> l_ptrs,
                          std::vector<ThetaPtr> subj_thetas,
                          PriorPtr hyper_likelihood, PriorPtr h_prior,
                          bool pop_debug = false, bool sub_debug = false);
    void run_hchains_fast_notiming(ThetaPtr phi_ptr, std::vector<LPtr> l_ptrs,
                                   std::vector<ThetaPtr> subj_thetas,
                                   PriorPtr hyper_likelihood, PriorPtr h_prior,
                                   bool pop_debug = false, bool sub_debug = false);

    // Hyper level with cache
    void crossover_hb_fast(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                           PriorPtr hyper_likelihood, PriorPtr h_prior,
                           bool debug = false, size_t para_idx = SIZE_T_MAX);
    void migration_hb_fast(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                           PriorPtr hyper_likelihood, PriorPtr h_prior,
                           bool debug = false, size_t para_idx = SIZE_T_MAX);

    void crossover_hb_fast_notiming(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                                    PriorPtr hyper_likelihood, PriorPtr h_prior,
                                    bool debug = false, size_t para_idx = SIZE_T_MAX);
    void migration_hb_fast_notiming(ThetaPtr phi_ptr, std::vector<ThetaPtr> subj_thetas,
                                    PriorPtr hyper_likelihood, PriorPtr h_prior,
                                    bool debug = false, size_t para_idx = SIZE_T_MAX);

    /* Step 3: Type Conversion Optimization - Subject level */
    void crossover_typeconv(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                           bool debug = false, size_t para_idx = SIZE_T_MAX);
    void migration_typeconv(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                           bool debug = false, size_t para_idx = SIZE_T_MAX);
    void run_chains_typeconv(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                            bool debug = false);
};
