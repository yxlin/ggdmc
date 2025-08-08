#pragma once

#include <RcppArmadillo.h>
#include <ggdmcHeaders/common_type_casting.h>
#include <ggdmcHeaders/likelihood_type_casting.h>

///////////////////////////////////////////////////
/* ------------- S4 classes  ------------- */
///////////////////////////////////////////////////
Rcpp::S4 new_posterior(const std::shared_ptr<theta_phi> &ptr,
                       const ThetaInput &theta_input)
{
    Rcpp::S4 out("posterior");
    out.slot("theta") = ptr->m_theta;
    out.slot("summed_log_prior") = ptr->m_lp;
    out.slot("log_likelihoods") = ptr->m_ll;

    out.slot("start") = ptr->m_start;
    out.slot("npar") = theta_input.nparameter;
    out.slot("pnames") = theta_input.pnames;
    out.slot("nmc") = theta_input.nmc;
    out.slot("thin") = theta_input.thin;
    out.slot("nchain") = theta_input.nchain;

    return out;
}
