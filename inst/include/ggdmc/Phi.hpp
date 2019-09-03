#ifndef PHI_HPP
#define PHI_HPP

#include <RcppArmadillo.h>

using namespace Rcpp;

class Theta
{
public:
  unsigned int m_start_R, m_store_i, m_nsamp, m_thin, m_nmc;

  arma::cube m_theta;

  arma::mat m_lp, m_ll, m_usetheta;
  arma::vec m_uselp, m_usell;

  Prior      * m_p;
  Likelihood * m_l;

  Theta(unsigned int nmc, unsigned int nchain, unsigned int npar,
        unsigned int thin, Prior * p, Likelihood * l) :
    m_thin(thin), m_nmc(nmc), m_p(p), m_l(l)
  // Start
  {
    using namespace arma;

    mat lp(nchain, nmc);  // sum log prior
    mat ll(nchain, nmc);  // sum log likelihood
    cube theta(npar, nchain, nmc);

    lp.fill(R_NegInf);
    ll.fill(R_NegInf);
    theta.fill(NA_REAL);

    vec tmp;

    for (size_t i = 0; i < nchain; i++)
    {
      size_t j = 1;

      while (lp.row(i).col(0).has_inf() || ll.row(i).col(0).has_inf())
      {
        tmp = p->rprior();

        theta.slice(0).col(i) = tmp;
        lp(i,0) = m_p->sumlogprior(tmp);
        ll(i,0) = m_l->sumloglike(tmp);
        j++;
        if (j > 1e4) { stop("Fail to set up new samples.");}
      }
    }

    /////////////////////////////////////////////////
    m_start_R = 1;
    m_store_i = 0;
    m_nsamp   = 1 + (nmc - 1) * thin;

    m_theta    = theta;
    m_lp       = lp;
    m_ll       = ll;
    m_usetheta = theta.slice(m_store_i);   // npar x nchain;
    m_uselp    = lp.col(m_store_i);        // nchains x 1
    m_usell    = ll.col(m_store_i);        // nchains x 1
  }

  Theta(List & samples, unsigned int nmc, unsigned int thin, Prior * p,
        Likelihood * l, bool add) : m_thin(thin), m_p(p), m_l(l)
  // Restart
  {
    using namespace arma;

    mat lp     = samples["summed_log_prior"];
    mat ll     = samples["log_likelihoods"];
    cube theta = samples["theta"];

    unsigned int npar      = samples["n.pars"];
    unsigned int nchain    = samples["n.chains"];
    unsigned int pnmc      = samples["nmc"];

    if (add) {nmc += pnmc;}

    m_theta = resize(theta, npar, nchain, nmc);
    m_lp    = resize(lp, nchain, nmc);
    m_ll    = resize(ll, nchain, nmc);

    if (add)
    {
      m_theta.slices(pnmc, nmc - 1).fill(NA_REAL);
      m_lp.cols(pnmc, nmc - 1).fill(R_NegInf);
      m_ll.cols(pnmc, nmc - 1).fill(R_NegInf);
      m_start_R = pnmc;
    }
    else
    {
      m_theta.fill(NA_REAL);
      m_lp.fill(R_NegInf);
      m_ll.fill(R_NegInf);

      m_theta.slice(0)   = theta.slice(pnmc - 1);
      m_lp.col(0)        = lp.col(pnmc - 1);
      m_ll.col(0)        = ll.col(pnmc - 1);
      m_start_R = 1;
    }
    ////////////////////////////////////////////////////
    m_nsamp   = 1 + (nmc - m_start_R) * thin;
    m_nmc     = nmc;

    m_store_i  = m_start_R-1;
    m_usetheta = m_theta.slice(m_store_i);   // npar x nchain;
    m_uselp    = m_lp.col(m_store_i);        // nchains x 1
    m_usell    = m_ll.col(m_store_i);        // nchains x 1
  }


  ~Theta()
  {
    delete m_p;
    delete m_l;
    // Rcout << "Theta destructor\n";
  }

  void store(unsigned int i, unsigned int & report, bool print=false)
  {
    if (i % m_thin == 0)
    {
      m_store_i++;
      if ( print && (m_store_i + 1) % report == 0 )
        Rcout << m_store_i + 1 << " ";

      m_lp.col(m_store_i) = m_uselp;   // nchain x nmc
      m_ll.col(m_store_i) = m_usell;
      m_theta.slice(m_store_i)  = m_usetheta; // npar x nchain x nmc
    }
  }

};

class Phi
{
private:
  double sumloghlike(Prior * p, std::vector<Theta *> & ts, unsigned int k)
  {
    double out = 0;
    for (size_t i = 0; i < ts.size(); i++)
    {
      out += p->sumlogprior(ts[i]->m_usetheta.col(k));
    }
    return out;
  }

public:
  unsigned int m_start_R, m_store_i, m_nsamp;

  unsigned int m_nsub, m_npar, m_nchain, m_thin, m_nmc;

  arma::cube m_phi0,    m_phi1;
  arma::mat  m_usephi0, m_usephi1, m_hlp, m_hll;
  arma::vec m_usehlp, m_usehll;
  Prior * m_p;
  Prior * m_lp;
  Prior * m_sp;

  Phi(unsigned int nmc,
      unsigned int nchain,
      unsigned int npar,
      unsigned int nsub,
      unsigned int thin,
      Prior * p, Prior * lp, Prior * sp,
      std::vector<Theta *> & ts) :
    m_nsub(nsub), m_npar(npar), m_nchain(nchain), m_thin(thin), m_nmc(nmc),
    m_p(p), m_lp(lp), m_sp(sp)
  {
    using namespace arma;
    cube phi0 (npar, nchain, nmc);
    cube phi1 (npar, nchain, nmc);
    mat  hlp  (nchain, nmc); // hyper sum log prior
    mat  hll  (nchain, nmc); // hyper log likelihood

    phi0.fill (NA_REAL);
    phi1.fill (NA_REAL);
    hlp.fill  (R_NegInf);
    hll.fill  (R_NegInf);

    // Rcout << "Generating hyper-start points for each chain: ";

    for (size_t i = 0; i < nchain; i++)
    {
      // Rcout << ".";

      phi0.slice(0).col(i) = m_lp->rprior();
      phi1.slice(0).col(i) = m_sp->rprior();
      m_p->m_p0 = phi0.slice(0).col(i);
      m_p->m_p1 = phi1.slice(0).col(i);
      hlp(i,0)  = m_lp->sumlogprior(m_p->m_p0) + m_sp->sumlogprior(m_p->m_p1);
      hll(i,0)  = sumloghlike(m_p, ts, i);
    }
    Rcout << std::endl;

    if (hll.col(0).has_inf())
    {
      Rcout << "\n The hyper parameters generated by hyper prior" <<
        " distributions result in infinite hyper likelihoods. " <<
          "You may try different location and scale prior distributions.\n";
      stop("Model fit stops.");
    }
    /////////////////////////////////////////////////
    m_start_R = 1;
    m_store_i = 0;
    m_nsamp   = 1 + (nmc - 1) * thin;
    m_phi0     = phi0;
    m_phi1     = phi1;
    m_hlp      = hlp;
    m_hll      = hll;

    m_usephi0 = phi0.slice(m_store_i);   // npar x nchain;
    m_usephi1 = phi1.slice(m_store_i);   // npar x nchain;
    m_usehlp  = hlp.col(m_store_i);      // nchains x 1
    m_usehll  = hll.col(m_store_i);      // nchains x 1
  }

  Phi(List & samples,
      unsigned int nmc,
      unsigned int nchain,
      unsigned int npar,
      unsigned int nsub,
      unsigned int thin,
      bool add,
      Prior * p, Prior * lp, Prior * sp) :
    m_nsub(nsub), m_npar(npar), m_nchain(nchain), m_thin(thin), m_p(p),
    m_lp(lp), m_sp(sp)
  // Restart
  {
    using namespace arma;

    List hyper  = samples.attr("hyper");
    List phi    = hyper["phi"];

    mat hlp   = hyper["h_summed_log_prior"];
    mat hll   = hyper["h_log_likelihoods"];
    cube phi0 = phi[0];
    cube phi1 = phi[1];

    unsigned int pnmc = hyper["nmc"];

    if (add) {nmc += pnmc;}

    m_phi0 = resize(phi0, npar, nchain, nmc);
    m_phi1 = resize(phi1, npar, nchain, nmc);
    m_hlp  = resize(hlp, nchain, nmc);
    m_hll  = resize(hll, nchain, nmc);

    if (add)
    {
      m_phi0.slices(pnmc, nmc - 1).fill(NA_REAL);
      m_phi1.slices(pnmc, nmc - 1).fill(NA_REAL);
      m_hlp.cols(pnmc, nmc - 1).fill(R_NegInf);
      m_hll.cols(pnmc, nmc - 1).fill(R_NegInf);
      m_start_R = pnmc;
    }
    else
    {
      m_phi0.fill(NA_REAL);
      m_phi1.fill(NA_REAL);
      m_hlp.fill(R_NegInf);
      m_hll.fill(R_NegInf);

      m_phi0.slice(0) = phi0.slice (pnmc - 1);
      m_phi1.slice(0) = phi1.slice (pnmc - 1);
      m_hlp.col(0)    = hlp.col    (pnmc - 1);
      m_hll.col(0)    = hll.col    (pnmc - 1);
      m_start_R = 1;
    }
    ////////////////////////////////////////////////////
    m_nsamp   = 1 + (nmc - m_start_R) * thin;
    m_nmc     = nmc;

    m_store_i = m_start_R-1;
    m_usephi0 = m_phi0.slice(m_store_i);   // npar x nchain;
    m_usephi1 = m_phi1.slice(m_store_i);   // npar x nchain;
    m_usehlp  = m_hlp.col(m_store_i);      // nchains x 1
    m_usehll  = m_hll.col(m_store_i);      // nchains x 1
  }


  ~Phi()
  {
    delete m_p;
    delete m_lp;
    delete m_sp;
    // Rcout << "Phi destructor\n";
  }

  void store(unsigned int i, unsigned int & report)
  {
    if (i % m_thin == 0)
    {
      m_store_i++;
      if ((m_store_i + 1) % report == 0) Rcout << m_store_i + 1 << " ";

      m_hlp.col(m_store_i)     = m_usehlp;   // nchain x nmc
      m_hll.col(m_store_i)     = m_usehll;
      m_phi0.slice(m_store_i)  = m_usephi0; // npar x nchain x nmc
      m_phi1.slice(m_store_i)  = m_usephi1; // npar x nchain x nmc
    }
  }


};

#endif
