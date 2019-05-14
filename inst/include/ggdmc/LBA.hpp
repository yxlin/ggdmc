#ifndef LBA_HPP
#define LBA_HPP

#include <RcppArmadillo.h>

class lba {
public:
  double m_A, m_b, m_mean_v, m_sd_v, m_t0, m_st0; // LBA distribution.
  bool is_posv;

  double *m_meanv_vec, *m_sdv_vec, *m_dt;
  unsigned int m_nmean_v, m_nrt;

  lba (double A, double b, double mean_v, double sd_v, double t0, double st0,
       bool posdrift, arma::vec & rt) :
    m_A(A), m_b(b), m_mean_v(mean_v), m_sd_v(sd_v), m_t0(t0), m_st0(st0),
    is_posv(posdrift)
  {
    m_nrt = rt.size();
    m_dt  = new double[m_nrt];

    for (size_t i = 0; i < m_nrt; i++)
    {
      m_dt[i] = rt[i] - m_t0;
    }

    denom = !is_posv ? 1.0 :
      R::fmax2(R::pnorm(m_mean_v/m_sd_v, 0, 1, 1, 0), 1e-10);

  };
  // pdf, cdf

  lba (double A, double b, double * mean_v, double * sd_v, double t0,
       double st0, unsigned int & nmean_v, bool posdrift) :
    m_A(A), m_b(b), m_t0(t0), m_st0(st0), is_posv(posdrift),
    m_meanv_vec(mean_v), m_sdv_vec(sd_v), m_nmean_v(nmean_v)
  {
    if (m_st0 < 0) Rcpp::stop("st0 must be greater than 0.");
  }
  // rlba_norm. NOTE double * mean_v and double * sd_v

  ~lba() {};

  arma::vec d()
  {
    arma::vec out(m_nrt);

    for (size_t i = 0; i < m_nrt; i++)
    {
      if (m_dt[i] < 0)
      {
         out[i] = 0.0;
      }
      else if (m_A < 1e-10)
      {
        out[i] = R::fmax2(0.0, (m_b / (m_dt[i]*m_dt[i])) * R::dnorm(
          m_b / m_dt[i], m_mean_v, m_sd_v, 0) / denom );
      }
      else
      {
        ts = m_dt[i] * m_sd_v;   // zs
        tv = m_dt[i] * m_mean_v; // zu
        term1 = m_mean_v *(R::pnorm((m_b-tv)/ts, 0, 1, 1, 0) -
          R::pnorm((m_b-m_A-tv)/ts,0,1,1, 0));
        term2 = m_sd_v   *(R::dnorm((m_b-m_A-tv)/ts, 0, 1, 0) -
          R::dnorm((m_b-tv)/ts, 0,1,0));
        out[i] = R::fmax2(0.0, (term1 + term2) / (m_A*denom));
      }

      if (ISNAN(out[i])) out[i] = 0.0;
    }

    delete [] m_dt;
    return out;
  }

  arma::vec p()
  {
    arma::vec out(m_nrt);

    for (size_t i = 0; i < m_nrt; i++)
    {
      if (m_A < 1e-10)
      {
        out[i] = R::fmin2(1.0,
                 R::fmax2(0.0,
                 R::pnorm(m_b / m_dt[i], m_mean_v, m_sd_v, false, false) / denom));
      }
      else
      {
        tv = m_dt[i] * m_mean_v;
        ts = m_dt[i] * m_sd_v;

        term1 = (m_b-m_A-tv) * R::pnorm((m_b - m_A - tv) / ts, 0, 1, true, false) -
                (m_b    -tv) * R::pnorm((m_b       - tv) / ts, 0, 1, true, false);
        term2 = ts * (R::dnorm((m_b - m_A - tv) / ts, 0, 1, false) -
                      R::dnorm((m_b       - tv) / ts, 0, 1, false));

        out[i] = R::fmin2(1.0,
                 R::fmax2(0.0, (1.0 + (term1 + term2)/m_A) / denom));
      }

      if (ISNAN(out[i])) out[i] = 0.0;

    }
    delete [] m_dt;
    return out;
  }


  void r (unsigned int & n, arma::mat & output)
  {
    // output n x 2
    arma::vec tmp(m_nmean_v);

    for (size_t i=0; i<n; i++)
    {
      rt(tmp);
      output(i, 0) = tmp.min();
      output(i, 1) = 1 + tmp.index_min(); // plus 1 to fit R indexing
    }
  }

  void print(const std::string & x) const
  {
    Rcpp::Rcout << x << "[A, b, mean_v, sd_v, t0]: " << m_A << ", " << m_b <<
      ", " << m_mean_v << ", " << m_sd_v << ", " << m_t0 << std::endl;
  }

  bool ValidateParams (bool print)
  {
    using namespace Rcpp;
    bool valid = true;

    if (m_A <= 0)   { valid = false; if (print) Rcout << "invalid parameter A = " << m_A << std::endl; }
    if (m_b < 0)    { valid = false; if (print) Rcout << "invalid parameter b = " << m_b << std::endl; }
    if (m_b < m_A)  { valid = false; if (print) Rcout << "b must greater than A. b = " << m_b << ", A = " << m_A << std::endl; }
    if (m_sd_v < 0) { valid = false; if (print) Rcout << "invalid parameter sd_v = " << m_sd_v << std::endl; }
    if (m_t0 < 0)   { valid = false; if (print) Rcout << "invalid parameter t0 = " << m_t0 << std::endl; }
    if (m_st0 < 0)  { valid = false; if (print) Rcout << "invalid parameter st0 = " << m_st0 << std::endl; }

    return valid;
  }


private:
  double tv, ts, term1, term2, denom;

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

arma::vec n1PDFfixedt0(arma::vec rt, arma::vec A, arma::vec b, arma::vec mean_v,
                       arma::vec sd_v, arma::vec t0, arma::vec st0,
                       bool posdrift);


#endif
