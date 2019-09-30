#ifndef PRIOR_HPP
#define PRIOR_HPP

#include <RcppArmadillo.h>

using namespace Rcpp;

class Prior
{
private:

  enum DistributionType { DEFAULT, TNORM, BETA_LU, GAMMA_L, LNORM_L, UNIF_,
                          CONSTANT, TNORM_TAU };

  DistributionType resolve_option(double type)
  {
    if (type == 1)  return TNORM;
    if (type == 2)  return BETA_LU;
    if (type == 3)  return GAMMA_L;
    if (type == 4)  return LNORM_L;
    if (type == 5)  return UNIF_;
    if (type == 6)  return CONSTANT;
    if (type == 7)  return TNORM_TAU;
    return DEFAULT;  // 0
  }

public:
  unsigned int m_npar;
  arma::vec m_p0, m_p1, m_l, m_u;
  arma::uvec m_d, m_lg;

  Prior (List & pprior);

  ~Prior();

  void dprior(double * pvector, double * out);
  // dprior is important for hierarchical modelling to be accurate

  arma::vec dprior(arma::vec pvector);
  // a wrapper for Armadillo vector

  arma::vec rprior();
  // Used by ininitlise & R's rprior;

  double sumlogprior(arma::vec pvector);

  void print(std::string str) const;
  // debugging function

};

#endif
