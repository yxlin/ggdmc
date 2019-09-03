#include <ggdmc.hpp>

using namespace Rcpp;

// --------------------- Initialization -----------------------------
// [[Rcpp::export]]
List init_new(List data,
              List prior,
              unsigned int nchain,
              unsigned int nmc,
              unsigned int thin,
              unsigned int report,
              double rp,
              double gammamult,
              double pm,
              double pm_old,
              bool block)
{

  unsigned int npar = prior.size();

  Design     * d0 = new Design (data);
  Prior      * p0 = new Prior (prior);
  Likelihood * l0 = new Likelihood (data, d0);
  Theta      * t0 = new Theta (nmc, nchain, npar, thin, p0, l0);
  Sampler    * s0 = new Sampler (nchain, npar, gammamult, rp);

  for (size_t i = 1; i < t0->m_nsamp; i++)
  {
    if (R::runif(0.0, 1.0) < pm_old)
    {
      s0->migrate_old(t0);
    }
    else if (R::runif(0.0, 1.0) < pm)
    {
      s0->migrate(t0);
    }
    else
    {
      if (block) { for (size_t j=0; j<npar; j++) s0->crossover(j, t0); }
      else       { s0->crossover(t0); }
    }

    t0->store(i, report, true);
  }
  Rcout << std::endl;

  std::vector<std::string> pnames(npar);
  for(size_t i=0; i < npar; i++) pnames[i] = d0->m_pnames[i];

  List out = List::create(
    Named("theta")            = t0->m_theta,
    Named("summed_log_prior") = t0->m_lp,
    Named("log_likelihoods")  = t0->m_ll,
    Named("data")             = data,
    Named("p.prior")          = prior,
    Named("start")            = t0->m_start_R,
    Named("n.pars")           = npar,
    Named("p.names")          = pnames,
    Named("nmc")              = nmc,
    Named("thin")             = t0->m_thin,
    Named("n.chains")         = nchain);

  delete t0;
  delete s0;
  return out;
}

// [[Rcpp::export]]
List init_old(List samples,
              unsigned int nmc,
              unsigned int thin,
              unsigned int report,
              double rp,
              double gammamult,
              double pm,
              double pm_old,
              bool block,
              bool add)
{
  List samples_in(clone(samples));
  List prior = samples["p.prior"];
  List data  = samples["data"];

  unsigned int nchain = samples["n.chains"];
  unsigned int npar   = samples["n.pars"];

  Design     * d0 = new Design (data);
  Prior      * p0 = new Prior (prior);
  Likelihood * l0 = new Likelihood (data, d0);
  Theta      * t0 = new Theta (samples_in, nmc, thin, p0, l0, add);
  Sampler    * s0 = new Sampler (nchain, npar, gammamult, rp);

  for (size_t i = 1; i < t0->m_nsamp; i++)
  {
    if (R::runif(0.0, 1.0) < pm_old)
    {
      s0->migrate_old(t0);
    }
    else if (R::runif(0.0, 1.0) < pm)
    {
      s0->migrate(t0);
    }
    else
    {
      if (block) { for (size_t j=0; j<npar; j++) s0->crossover(j, t0); }
      else       { s0->crossover(t0); }
    }

    t0->store(i, report, true);
  }
  Rcout << std::endl;


  std::vector<std::string> pnames(npar);
  for(size_t i=0; i < npar; i++) pnames[i] = d0->m_pnames[i];
  List out = List::create(
    Named("theta")            = t0->m_theta,
    Named("summed_log_prior") = t0->m_lp,
    Named("log_likelihoods")  = t0->m_ll,
    Named("data")             = data,
    Named("p.prior")          = prior,
    Named("start")            = t0->m_start_R,
    Named("n.pars")           = npar,
    Named("p.names")          = pnames,
    Named("nmc")              = t0->m_nmc,
    Named("thin")             = t0->m_thin,
    Named("n.chains")         = nchain);

  delete t0;
  delete s0;
  return out;
}

// --------------------- Hierarchical versions -----------------------------
static void update_priors(Theta * t, Phi * phi)
{
  for (size_t i = 0; i < phi->m_nchain; i++)
  {
    phi->m_p->m_p0 = phi->m_usephi0.col(i);
    phi->m_p->m_p1 = phi->m_usephi1.col(i);
    t->m_p->m_p0   = phi->m_usephi0.col(i);
    t->m_p->m_p1   = phi->m_usephi1.col(i);
    t->m_uselp[i]  = phi->m_p->sumlogprior(t->m_usetheta.col(i));
  }
}

// [[Rcpp::export]]
List init_newhier(List prior,
                  List lprior,
                  List sprior,
                  List data,
                  unsigned int nchain,
                  unsigned int nmc,
                  unsigned int thin,
                  unsigned int report,
                  double rp,
                  double gammamult,
                  double pm,
                  double pm_old,
                  bool block)
{
  unsigned int npar = prior.size();
  unsigned int nsub = data.size();

  Prior * p0 = new Prior (prior);
  Prior * lp = new Prior (lprior);
  Prior * sp = new Prior (sprior);

  std::vector<Design     *> ds(nsub);
  std::vector<Prior      *> ps(nsub);
  std::vector<Likelihood *> ls(nsub);
  std::vector<Theta      *> ts(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    List datai = data[i];  // Must cast out first
    ds[i] = new Design (datai);
    ps[i] = new Prior (prior);
    ls[i] = new Likelihood (datai, ds[i]); // diff. RTs, Rs for diff. subjs
    ts[i] = new Theta (nmc, nchain, npar, thin, ps[i], ls[i]);
  }

  Phi     * phi = new Phi (nmc, nchain, npar, nsub, thin, p0, lp, sp, ts);
  Sampler * s0  = new Sampler (nchain, npar, gammamult, rp);
  Rcout << "Start sampling: ";

  for (size_t i = 1; i < phi->m_nsamp; i++)
  {
    if (R::runif(0.0, 1.0) < pm_old)
    {
      s0->migrate_old(phi, ts);
    }
    else if (R::runif(0.0, 1.0) < pm)
    {
      s0->migrate(phi, ts);
    }
    else
    {
      for (size_t j = 0; j < npar; j++) s0->crossover(j, phi, ts);
    }

    /////////////////////////////////////////////////////////////////
    for (size_t k = 0; k < nsub; k++)
    {
      update_priors(ts[k], phi);

      if (R::runif(0.0, 1.0) < pm_old)
      {
        s0->migrate_old(ts[k]);
      }
      else if (R::runif(0.0, 1.0) < pm)
      {
        s0->migrate(ts[k]);
      }
      else if (block)
      {
        for (size_t j=0; j<npar; j++) s0->crossover(j, ts[k]);
      }
      else
      {
        s0->crossover(ts[k]);
      }

      ts[k]->store(i, report);
    }

    phi->store(i, report);
  }

  Rcout << std::endl;

  ////////////////////////////////////////////////////////////////

  std::vector<std::string> pnames = prior.attr("names");
  List out(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    out[i] = List::create(
      Named("theta")            = ts[i]->m_theta,
      Named("summed_log_prior") = ts[i]->m_lp,
      Named("log_likelihoods")  = ts[i]->m_ll,
      Named("data")             = data[i],
      Named("p.prior")          = prior,
      Named("start")            = ts[i]->m_start_R,
      Named("n.pars")           = npar,
      Named("p.names")          = pnames,
      Named("nmc")              = nmc,
      Named("thin")             = thin,
      Named("n.chains")         = nchain);
  }

  List phi_tmp = List::create(
    Named("location") = phi->m_phi0,
    Named("scale")    = phi->m_phi1);

  List ppprior = List::create(
    Named("location") = lprior,
    Named("scale")    = sprior);

  List hyper = List::create(   // 16 elements
    Named("phi")                = phi_tmp,
    Named("h_summed_log_prior") = phi->m_hlp,
    Named("h_log_likelihoods")  = phi->m_hll,
    Named("pp.prior")           = ppprior,
    Named("start")              = phi->m_start_R,
    Named("n.pars")             = npar,
    Named("p.names")            = pnames,
    Named("rp")                 = rp,
    Named("nmc")                = nmc,
    Named("thin")               = thin,
    Named("n.chains")           = nchain);

  out.attr("hyper") = hyper;

  delete s0;
  delete phi;
  for (size_t i = 0; i < nsub; i++) delete ts[i];
  return out;
}

// [[Rcpp::export]]
List init_oldhier(List samples,
                  unsigned int nmc,
                  unsigned int thin,
                  unsigned int report,
                  double rp,
                  double gammamult,
                  double pm,
                  double pm_old,
                  bool block,
                  bool add)
{
  List samples_in(clone(samples));

  List hyper  = samples_in.attr("hyper");
  List pprior = hyper ["pp.prior"];
  List lprior = pprior["location"];
  List sprior = pprior["scale"];

  List subject0 = samples_in[0];
  List prior    = subject0["p.prior"];

  unsigned int npar   = prior.size();
  unsigned int nsub   = samples.size();
  unsigned int nchain = hyper["n.chains"];

  Prior * p0 = new Prior (prior);
  Prior * lp = new Prior (lprior);
  Prior * sp = new Prior (sprior);

  std::vector<Design     *> ds(nsub);
  std::vector<Prior      *> ps(nsub);
  std::vector<Likelihood *> ls(nsub);
  std::vector<Theta      *> ts(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    List subjecti = samples_in[i];
    List datai    = subjecti["data"];  // Must cast out first

    ds[i] = new Design (datai);
    ps[i] = new Prior (prior);
    ls[i] = new Likelihood (datai, ds[i]); // diff. RTs, Rs for diff. subjs
    ts[i] = new Theta (subjecti, nmc, thin, ps[i], ls[i], add);
  }

  Phi     * phi = new Phi (samples_in, nmc, nchain, npar, nsub, thin, add, p0,
                           lp, sp);
  Sampler * s0  = new Sampler (nchain, npar, gammamult, rp);
  Rcout << "Start sampling: ";

  for (size_t i = 1; i < phi->m_nsamp; i++)
  {
    if (R::runif(0.0, 1.0) < pm_old)
    {
      s0->migrate_old(phi, ts);
    }
    else if (R::runif(0.0, 1.0) < pm)
    {
      s0->migrate(phi, ts);
    }
    else
    {
      for (size_t j = 0; j < npar; j++) s0->crossover(j, phi, ts);
    }

    /////////////////////////////////////////////////////////////////
    for (size_t k = 0; k < nsub; k++)
    {
      update_priors(ts[k], phi);

      if (R::runif(0.0, 1.0) < pm_old)
      {
        s0->migrate_old(ts[k]);
      }
      else if (R::runif(0.0, 1.0) < pm)
      {
        s0->migrate(ts[k]);
      }
      else if (block)
      {
        for (size_t j=0; j<npar; j++) s0->crossover(j, ts[k]);
      }
      else
      {
        s0->crossover(ts[k]);
      }

      ts[k]->store(i, report);
    }

    phi->store(i, report);
  }
  Rcout << std::endl;

  ////////////////////////////////////////////////////////////////////////
  std::vector<std::string> pnames = prior.attr("names");
  List out(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    List subjecti = samples_in[i];
    List datai    = subjecti["data"];  // Must cast out first

    out[i] = List::create(
      Named("theta")            = ts[i]->m_theta,
      Named("summed_log_prior") = ts[i]->m_lp,
      Named("log_likelihoods")  = ts[i]->m_ll,
      Named("data")             = datai,
      Named("p.prior")          = prior,
      Named("start")            = ts[i]->m_start_R,
      Named("n.pars")           = npar,
      Named("p.names")          = pnames,
      Named("nmc")              = ts[i]->m_nmc,
      Named("thin")             = ts[i]->m_thin,
      Named("n.chains")         = nchain);
  }

  List phi_tmp = List::create(
    Named("location") = phi->m_phi0,
    Named("scale")    = phi->m_phi1);

  List ppprior = List::create(
    Named("location") = lprior,
    Named("scale")    = sprior);

  hyper = List::create(   // 16 elements
    Named("phi")                = phi_tmp,
    Named("h_summed_log_prior") = phi->m_hlp,
    Named("h_log_likelihoods")  = phi->m_hll,
    Named("pp.prior")           = ppprior,
    Named("start")              = phi->m_start_R,
    Named("n.pars")             = npar,
    Named("p.names")            = pnames,
    Named("rp")                 = rp,
    Named("nmc")                = phi->m_nmc,
    Named("thin")               = phi->m_thin,
    Named("n.chains")           = nchain);

  out.attr("hyper") = hyper;

  delete s0;
  delete phi;
  for (size_t i = 0; i < nsub; i++) delete ts[i];
  return out;

}
