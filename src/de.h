#pragma once

#include <ggdmcHeaders/likelihood.h>
#include <ggdmcHeaders/theta.h>

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

    /* Subject level only */
    void crossover(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                   bool debug = false, size_t para_idx = SIZE_T_MAX);
    void migration(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                   bool debug = false, size_t para_idx = SIZE_T_MAX);

    void run_chains(ThetaPtr t_ptr, PriorPtr p_ptr, LPtr l_ptr,
                    bool debug = false);

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
};
