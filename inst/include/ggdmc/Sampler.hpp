#ifndef SAMPLER_HPP
#define SAMPLER_HPP

#include <RcppArmadillo.h>

using namespace Rcpp;

class Sampler
{
private:
  unsigned int m_npar, m_nchain, m_nsubchain;
  double m_gammamult, m_rp, m_ga;

  double cur_logpos, tmp_logpos, tmp_lp, tmp_ll, mh;
  arma::vec m_gamma, tmp0, tmp1;

  arma::uvec PickChains(unsigned int k, unsigned int nchain, arma::uvec chains)
  {
    chains.shed_row(k);
    arma::uvec rchains = arma::shuffle(chains);
    return rchains.rows(0, nchain - 1);
  }

  arma::uvec GetSubchains()
  {
    m_nsubchain = (unsigned int)std::ceil((double)m_nchain * R::runif(0.0, 1.0));
    arma::uvec chains  = arma::shuffle(m_chains);
    return arma::sort(chains.rows(0, m_nsubchain - 1));
  }

  double sumloghlike(arma::vec loc, arma::vec sca,
                     Phi * phi, std::vector<Theta *> & ts,
                     unsigned int k)
  {
    phi->m_p->m_p0 = loc;
    phi->m_p->m_p1 = sca;

    double out = 0;
    for (size_t i = 0; i < ts.size(); i++)
    {
      out += phi->m_p->sumlogprior(ts[i]->m_usetheta.col(k));
    }
    return out;
  }

public:
  arma::uvec m_chains, m_subchains;

  Sampler(unsigned int nchain,
          unsigned int npar,
          double gammamult,
          double rp) :
    m_npar(npar), m_nchain(nchain), m_gammamult(gammamult), m_rp(rp)
  {
    using namespace arma;

    double tmp = m_gammamult / sqrt( (double)(2 * m_npar) );
    vec tmp_gamma(m_npar); tmp_gamma.fill(tmp);
    m_gamma = tmp_gamma;

    m_chains = linspace<uvec>(0, m_nchain - 1, m_nchain);
    m_ga = m_gammamult / std::sqrt(4.0*m_npar);

  }

  ~Sampler()
  {
    // Rcout << "Sampler destructor\n";
  }

  void crossover(unsigned int i, Phi * phi, std::vector<Theta *> & ts)
  // crossover update one parameter at a time
  {

    for (size_t j = 0; j < m_nchain; j++)
    {
      // -------------------------------- update h like ----------------
      phi->m_usehll[ m_chains[j] ] = sumloghlike(phi->m_usephi0.col(m_chains[j]),
                                               phi->m_usephi1.col(m_chains[j]),
                                               phi, ts, m_chains[j]);
      cur_logpos = phi->m_usehlp[ m_chains[j] ] + phi->m_usehll[ m_chains[j]] ;
      // -------------------------------- update h like ----------------

      // -------------------------------- tmp logpos ----------------
      m_subchains = PickChains(m_chains[j], 2, m_chains);
      tmp0 = phi->m_usephi0.col(m_chains[j]);
      tmp1 = phi->m_usephi1.col(m_chains[j]);

      tmp0[i] = R::runif(-m_rp, m_rp) + phi->m_usephi0(i, m_chains[j]) + m_ga *
        ( phi->m_usephi0(i, m_subchains[0]) - phi->m_usephi0(i, m_subchains[1]) );

      tmp1[i] = R::runif(-m_rp, m_rp) + phi->m_usephi1(i, m_chains[j]) + m_ga *
        ( phi->m_usephi1(i, m_subchains[0]) - phi->m_usephi1(i, m_subchains[1]) );

      tmp_lp = phi->m_lp->sumlogprior(tmp0) +
               phi->m_sp->sumlogprior(tmp1);
      tmp_ll = sumloghlike( tmp0, tmp1, phi, ts, m_chains[j] );
      tmp_logpos = tmp_lp + tmp_ll;
    // -------------------------------- tmp logpos ----------------

      mh = std::exp(tmp_logpos - cur_logpos);

      if (!ISNAN(mh) && (R::runif(0, 1) < mh) )
      {
        phi->m_usephi0(i, m_chains[j]) = tmp0[i];
        phi->m_usephi1(i, m_chains[j]) = tmp1[i];
        phi->m_usehlp[m_chains[j]]     = tmp_lp;
        phi->m_usehll[m_chains[j]]     = tmp_ll;
      }

    }    // chain loop

  }

  void crossover(unsigned int i, Theta * t)
  {

    for (size_t j = 0; j < m_nchain; j++)
    {
      cur_logpos = t->m_uselp[m_chains[j]] + t->m_usell[m_chains[j]];

      // -------------------------------- tmp logpos ----------------
      m_subchains = PickChains(m_chains[j], 2, m_chains);
      tmp0 = t->m_usetheta.col(m_chains[j]);

      tmp0[i] = R::runif(-m_rp, m_rp) + t->m_usetheta(i, m_chains[j]) + m_ga *
        ( t->m_usetheta(i, m_subchains[0]) - t->m_usetheta(i, m_subchains[1]) );

      tmp_lp = t->m_p->sumlogprior(tmp0);
      tmp_ll = t->m_l->sumloglike (tmp0);
      tmp_logpos = tmp_ll + tmp_lp;
      // -------------------------------- tmp logpos ----------------

      mh = std::exp(tmp_logpos - cur_logpos);

      if ( !ISNAN(mh) && (R::runif(0, 1) < mh) )
      {
        t->m_usetheta(i, m_chains[j]) = tmp0[i];
        t->m_uselp[m_chains[j]]       = tmp_lp;
        t->m_usell[m_chains[j]]       = tmp_ll;
      }
    }

  }

  void crossover(Theta * t)
  {
    arma::vec noise(m_npar);

    for (size_t i = 0; i < m_nchain; i++)
    {
      cur_logpos = t->m_usell[ m_chains[i] ] + t->m_uselp[ m_chains[i] ];

      m_subchains = PickChains(m_chains[i], 2, m_chains);
      for(size_t j = 0; j < m_npar; j++) noise[j] = R::runif(-m_rp, m_rp);

      tmp0 = noise + t->m_usetheta.col( m_chains[i] ) +
        (m_gamma % ( t->m_usetheta.col( m_subchains[0] ) -
                     t->m_usetheta.col( m_subchains[1] ) ));
      tmp_lp = t->m_p->sumlogprior(tmp0);
      tmp_ll = t->m_l->sumloglike (tmp0);

      tmp_logpos = tmp_lp + tmp_ll;

      // Rcout << "tmp_lp & tmp_ll " << "[" << tmp_lp << " " << tmp_ll << "]" << "\n";
      // if (ISNAN(tmp_logpos)) tmp_logpos = R_NegInf;

      mh = std::exp(tmp_logpos - cur_logpos);
      if ( !ISNAN(mh) && (R::runif(0, 1) < mh) )
      {
        t->m_usetheta.col(m_chains[i]) = tmp0;
        t->m_uselp[ m_chains[i] ]      = tmp_lp;
        t->m_usell[ m_chains[i] ]      = tmp_ll;
      }

    }
  }

  void migrate(Phi * phi, std::vector<Theta *> & ts)
  {
    arma::vec tmp_loc(m_npar), tmp_sca(m_npar);
    m_subchains = GetSubchains();     // eg, 0, 1, 3, 4, 8;
    unsigned int next_chain;
    // , nsubchain = m_subchains.size();

    for (size_t i = 0; i < m_nsubchain; i++)
    {
      next_chain = ((i+1) == m_nsubchain) ? m_subchains(0) : m_subchains(i+1);

      phi->m_usehll[ next_chain ] = sumloghlike(phi->m_usephi0.col(next_chain),
                                                phi->m_usephi1.col(next_chain),
                                                phi, ts, next_chain);
      cur_logpos = phi->m_usehlp[ next_chain ] + phi->m_usehll[ next_chain ];

      ///////////////////////////////////////////////////////////////////////
      for (size_t j = 0; j < m_npar; j++)
      {
        tmp_loc(j) = phi->m_usephi0(j, m_subchains[i]) +
          R::rnorm(phi->m_usephi0(j, m_subchains[i]), m_rp);
        tmp_sca(j) = phi->m_usephi1(j, m_subchains[i]) +
          R::rnorm(phi->m_usephi1(j, m_subchains[i]), m_rp);
      }

      tmp_lp     = phi->m_lp->sumlogprior(tmp_loc) +
                   phi->m_sp->sumlogprior(tmp_sca);
      tmp_ll     = sumloghlike(tmp_loc, tmp_sca, phi, ts, m_subchains[i]);
      tmp_logpos = tmp_lp + tmp_ll;

      mh = std::exp(tmp_logpos - cur_logpos);

      if ( !ISNAN(mh) && (R::runif(0, 1) < mh ) )
      {
        phi->m_usephi0.col(next_chain) = tmp_loc;
        phi->m_usephi0.col(next_chain) = tmp_sca;
        phi->m_usehlp[next_chain]      = tmp_lp;
        phi->m_usehll[next_chain]      = tmp_ll;
      }
    }

  }

  void migrate(Theta * t)
  {
    arma::vec tmp(m_npar);
    m_subchains = GetSubchains();     // eg, 0, 1, 3, 4, 8;
    unsigned int next_chain;
    // , nsubchain = m_subchains.size();

    for(size_t i = 0; i < m_nsubchain; i++)
    {
      next_chain = ((i+1) == m_nsubchain) ? m_subchains[0] : m_subchains[i+1];

      for (size_t j = 0; j < m_npar; j++)
      {
        tmp(j) = t->m_usetheta(j, m_subchains[i]) +
          R::rnorm(t->m_usetheta(j, m_subchains[i]), m_rp);
      }

      tmp_lp     = t->m_p->sumlogprior(tmp);
      tmp_ll     = t->m_l->sumloglike (tmp);
      tmp_logpos = tmp_lp + tmp_ll;

      if (ISNAN(tmp_logpos)) tmp_logpos = R_NegInf;

      cur_logpos = t->m_uselp[next_chain] + t->m_usell[next_chain];

      mh = std::exp(tmp_logpos - cur_logpos);

      if ( !ISNAN(mh) && (R::runif(0, 1) < mh) )
      {
        t->m_usetheta.col(next_chain) = tmp;
        t->m_uselp[next_chain]        = tmp_lp;
        t->m_usell[next_chain]        = tmp_ll;
      }
    }

  }

  void migrate_old(Phi * phi, std::vector<Theta *> & ts)
  {
    using namespace arma;

    double cur_logpos, tmp_logpos, mh;

    m_subchains = GetSubchains();
    // unsigned int nsubchain = subchains.size();

    mat tmp_loc(m_npar, m_nsubchain), tmp_sca(m_npar, m_nsubchain);
    vec cur_lp(m_nsubchain), cur_ll(m_nsubchain),
    tmp_lp(m_nsubchain), tmp_ll(m_nsubchain), noise(m_npar);

    for(size_t i = 0; i < m_nsubchain; i++)
    {
      // -------------------------------- update h like ----------------
      phi->m_usehll[m_subchains[i]] = sumloghlike(
        phi->m_usephi0.col(m_subchains[i]),
        phi->m_usephi1.col(m_subchains[i]),
        phi, ts, m_subchains[i]);
      cur_lp[i]  = phi->m_usehlp[m_subchains[i]];
      cur_ll[i]  = phi->m_usehll[m_subchains[i]];
      // -------------------------------- update h like ----------------

      // -------------------------------- tmp logpos ----------------
      for(size_t j = 0; j < m_npar; j++) noise[j] = R::runif(-m_rp, m_rp);
      tmp_loc.col(i) = phi->m_usephi0.col(m_subchains[i]) + noise; // proposal
      tmp_sca.col(i) = phi->m_usephi1.col(m_subchains[i]) + noise; // proposal

      tmp_lp[i] = phi->m_lp->sumlogprior(tmp_loc.col(i)) +
        phi->m_sp->sumlogprior(tmp_sca.col(i));
      tmp_ll[i] = sumloghlike(tmp_loc.col(i), tmp_sca.col(i), phi, ts, m_subchains[i]);
    }

    tmp_logpos = tmp_ll[m_nsubchain - 1] + tmp_lp[m_nsubchain - 1];
    cur_logpos = cur_ll[0] + cur_lp[0];
    mh = std::exp(tmp_logpos - cur_logpos);

    if ( !ISNAN(mh) && (R::runif(0, 1) < mh ) )
    {
      phi->m_usephi0.col(m_subchains[0]) = tmp_loc.col(m_nsubchain - 1);
      phi->m_usephi1.col(m_subchains(0)) = tmp_sca.col(m_nsubchain - 1);
      phi->m_usehlp[m_subchains(0)]      = tmp_lp[m_nsubchain - 1];
      phi->m_usehll[m_subchains(0)]      = tmp_ll[m_nsubchain - 1];
    }

    if (m_nsubchain != 1)
    {
      for(size_t k = 0; k < (m_nsubchain - 2); k++)
      {
        tmp_logpos = tmp_ll[k] + tmp_lp[k];
        cur_logpos = cur_ll[k + 1] + cur_lp[k + 1];
        mh = std::exp(tmp_logpos - cur_logpos);

        if ( !ISNAN(mh) && (R::runif(0, 1) < mh ) )
        {
          phi->m_usephi0.col(m_subchains(k + 1)) = tmp_loc.col(k);
          phi->m_usephi1.col(m_subchains(k + 1)) = tmp_sca.col(k);
          phi->m_usehlp(m_subchains(k + 1))      = tmp_lp[k];
          phi->m_usehll(m_subchains(k + 1))      = tmp_ll[k];
        }
      }
    }
  }

  void migrate_old(Theta * t)
  {
    using namespace arma;

    double cur_logpos, tmp_logpos, mh;

    m_subchains = GetSubchains();
    // unsigned int nsubchain = subchains.size();

    mat tmp(m_npar, m_nsubchain);
    vec cur_lp(m_nsubchain), cur_ll(m_nsubchain), tmp_lp(m_nsubchain),
        tmp_ll(m_nsubchain), noise(m_npar);

    for(size_t i = 0; i < m_nsubchain; i++)
    {
      for(size_t j = 0; j < m_npar; j++) noise[j] = R::runif(-m_rp, m_rp);

      // npar x nchain;
      tmp.col(i) = t->m_usetheta.col(m_subchains[i]) + noise; // proposal
      cur_lp[i]  = t->m_uselp(m_subchains[i]);
      cur_ll[i]  = t->m_usell(m_subchains[i]);

      tmp_lp[i] = t->m_p->sumlogprior(tmp.col(i));
      tmp_ll[i] = t->m_l->sumloglike (tmp.col(i));
    }

    tmp_logpos = tmp_ll[m_nsubchain - 1] + tmp_lp[m_nsubchain - 1];
    cur_logpos = cur_ll[0] + cur_lp[0];   // migrate to the first subchain

    mh = std::exp(tmp_logpos - cur_logpos);

    if (ISNAN(tmp_logpos)) tmp_logpos = R_NegInf;

    if (!ISNAN(tmp_logpos) && R::runif(0, 1) < mh)
    {
      t->m_usetheta.col(m_subchains[0]) = tmp.col(m_nsubchain - 1);
      t->m_uselp[m_subchains[0]]        = tmp_lp[m_nsubchain - 1];
      t->m_usell[m_subchains[0]]        = tmp_ll[m_nsubchain - 1];
    }

    if (m_nsubchain != 1)
    {
      for(size_t k = 1; k < (m_nsubchain - 1); k++)
      {
        tmp_logpos = tmp_ll(k) + tmp_lp(k);
        cur_logpos = cur_ll(k + 1) + cur_lp(k + 1);

        mh = std::exp(tmp_logpos - cur_logpos);
        if (ISNAN(tmp_logpos)) tmp_logpos = R_NegInf;

        if (!ISNAN(tmp_logpos) && R::runif(0, 1) < mh)
        {
          t->m_usetheta.col(m_subchains[k+1]) = tmp.col(k);
          t->m_uselp[m_subchains[k+1]]        = tmp_lp[k];
          t->m_usell[m_subchains[k+1]]        = tmp_ll[k];
        }
      }
    }
  }


};


#endif
