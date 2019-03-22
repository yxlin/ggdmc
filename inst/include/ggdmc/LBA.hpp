#ifndef LBA_HPP
#define LBA_HPP

#include <RcppArmadillo.h>

class lba {
public:
  double m_A, m_b, m_mean_v, m_sd_v, m_t0, m_st0; // LBA distribution.
  bool is_posv;

  double *m_meanv_vec, *m_sdv_vec;
  unsigned int m_nmean_v;

  lba (double A, double b, double mean_v, double sd_v, double t0,
       bool posdrift);
  // pdf, cdf

  lba (double A, double b, double * mean_v, double * sd_v, double t0,
       double st0, unsigned int & nmean_v, bool posdrift);
    // rlba_norm. NOTE double * mean_v and double * sd_v

  ~lba();


  void d (std::vector<double> & x, double * output);
  void d (std::vector<double> & x, std::vector<double> & output);
  void d (arma::vec & x, arma::vec & output);

  void p (std::vector<double> & x, double * output);
  void p (std::vector<double> & x, std::vector<double> & output);
  void p (arma::vec & x, arma::vec & output);

  void node1_pdf (arma::vec & x, arma::vec & output);
  void r (unsigned int & n, arma::mat & output);
  void print(const std::string & x) const;

private:
  double dt;

  double remove_t0 (double & x)
  {
    return x - m_t0; // allow negative DT for MCMC converge.
    // if (x < 0) return 0;
    // else       return (x - m_t0);
  }

  double d (double & x)
    // Return probability density function. fpt-pdf
    // NOTE: For unclear reasons, we need to permit unreasonable parameters, like
    // t0 < 0 to proceed to produce very small likelihoods, such as 3.798281e-80.
    // This allows MCMC to converge.
  {
    // if (m_b < m_A) return 1e-10;
    dt = remove_t0(x);

    double xx=m_b/dt, tv=dt*m_mean_v, ts=dt*m_sd_v, denom, term1, term2;

    denom = !is_posv ? 1.0 :
      std::max(R::pnorm(m_mean_v/m_sd_v, 0, 1, true, false), 1e-10);

    term1 = m_mean_v * (R::pnorm((m_b-tv)    / ts, 0, 1, 1, 0) -
                        R::pnorm((m_b-m_A-tv)/ ts, 0, 1, 1, 0));
    term2 = m_sd_v   * (R::dnorm((m_b-m_A-tv)/ ts, 0, 1, 0)  -
                        R::dnorm((m_b-tv)    / ts, 0, 1, 0));

    if (m_A < 1e-10) {
      return std::max(0.0,
                      (m_b/(dt*dt)) * R::dnorm(xx, m_mean_v, m_sd_v, 0) / denom);
    } else {
      return std::max(0.0,
                      (term1 + term2)/ (m_A*denom));
    }
  }

  double p (double & x)
    // Return cumulative distribution function. fpt-cdf
  {
    dt = remove_t0(x);

    double tv=dt*m_mean_v, ts=dt*m_sd_v, denom, term1, term2;

    denom = !is_posv ? 1.0 :
      std::max(R::pnorm(m_mean_v/m_sd_v, 0, 1, true, false), 1e-10);

    term1 = (m_b-m_A-tv) * R::pnorm((m_b - m_A - tv) / ts, 0, 1, true, false) -
            (m_b-tv)     * R::pnorm((m_b       - tv) / ts, 0, 1, true, false);
    term2 = ts * (R::dnorm((m_b - m_A - tv) / ts, 0, 1, false) -
                  R::dnorm((m_b       - tv) / ts, 0, 1, false));
    if (m_A < 1e-10) {
      return std::min(1.0,
             std::max(0.0,
                      R::pnorm(dt, m_mean_v, m_sd_v, false, false) / denom));
    } else {
      return std::min(1.0,
             std::max(0.0, (1.0 + (term1 + term2)/m_A) / denom));
    }
  }

  void rt(arma::vec & output)
    // Return one random set of response times. The size of a set equals to the
    // number of accumulators
  {
    double lower;
    for (size_t i=0; i<m_nmean_v; i++)
    {
      lower =  is_posv ? 0 : R_NegInf;
      tnorm * obj = new tnorm(m_meanv_vec[i], m_sdv_vec[i], lower, R_PosInf);
      output[i] = m_t0 + (m_b - m_A * R::runif(0, 1)) / obj->r();
      delete obj;
    }

    if ( output.has_inf() ) Rcpp::stop("Found infinite in lba class");


  }


};

arma::vec n1PDFfixedt0(arma::vec rt, arma::vec A, arma::vec b,
                       arma::vec mean_v, arma::vec sd_v, arma::vec t0,
                       bool posdrift);
// no st0, so fixed t0

#endif
