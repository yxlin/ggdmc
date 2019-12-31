//    Copyright (C) <2019>  <Yi-Shin Lin>
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#include <ggdmc.hpp>

using namespace Rcpp;

// --------------------- Initialization -----------------------------
// [[Rcpp::export]]
S4 init_new(S4 dmi,
            S4 prior,
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
  unsigned int npar = prior.slot("npar");

  Design     * d0 = new Design (dmi);
  Prior      * p0 = new Prior (prior);
  Likelihood * l0 = new Likelihood (dmi, d0, 3.0);  // precision = 3.0 for DDM
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

  S4 out ("posterior");
  out.slot("theta")            = t0->m_theta;
  out.slot("summed_log_prior") = t0->m_lp;
  out.slot("log_likelihoods")  = t0->m_ll;
  out.slot("dmi")              = dmi;
  out.slot("prior")            = prior;
  out.slot("start")            = t0->m_start_R;
  out.slot("npar")             = npar;
  out.slot("pnames")           = pnames;
  out.slot("nmc")              = nmc;
  out.slot("thin")             = t0->m_thin;
  out.slot("nchain")           = nchain;


  // List out = List::create(
  //   Named("theta")            = t0->m_theta,
  //   Named("summed_log_prior") = t0->m_lp,
  //   Named("log_likelihoods")  = t0->m_ll,
  //   Named("data")             = data,
  //   Named("p.prior")          = prior,
  //   Named("start")            = t0->m_start_R,
  //   Named("n.pars")           = npar,
  //   Named("p.names")          = pnames,
  //   Named("nmc")              = nmc,
  //   Named("thin")             = t0->m_thin,
  //   Named("n.chains")         = nchain);

  delete t0;
  delete s0;
  return out;
}

// [[Rcpp::export]]
S4 init_old(S4 samples,
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
  S4 samples_in(clone(samples));
  S4 dmi   = samples.slot("dmi");
  S4 prior = samples.slot("prior");

  unsigned int nchain = samples.slot("nchain");
  unsigned int npar   = samples.slot("npar");

  Design     * d0 = new Design (dmi);
  Prior      * p0 = new Prior (prior);
  Likelihood * l0 = new Likelihood (dmi, d0, 3.0);  // precision = 3.0 for DDM
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

  S4 out ("posterior");
  out.slot("theta")            = t0->m_theta;
  out.slot("summed_log_prior") = t0->m_lp;
  out.slot("log_likelihoods")  = t0->m_ll;
  out.slot("dmi")              = dmi;
  out.slot("prior")            = prior;
  out.slot("start")            = t0->m_start_R;
  out.slot("npar")             = npar;
  out.slot("pnames")           = pnames;
  out.slot("nmc")              = t0->m_nmc;
  out.slot("thin")             = t0->m_thin;
  out.slot("nchain")           = nchain;

  // List out = List::create(
  //   Named("theta")            = t0->m_theta,
  //   Named("summed_log_prior") = t0->m_lp,
  //   Named("log_likelihoods")  = t0->m_ll,
  //   Named("data")             = data,
  //   Named("p.prior")          = prior,
  //   Named("start")            = t0->m_start_R,
  //   Named("n.pars")           = npar,
  //   Named("p.names")          = pnames,
  //   Named("nmc")              = t0->m_nmc,
  //   Named("thin")             = t0->m_thin,
  //   Named("n.chains")         = nchain);

  delete t0; // remove prior and likelihood objects. likelihood will remove design object
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
S4 init_newhier(S4 prior,
                S4 lprior,
                S4 sprior,
                List dmi, // a list of multiple dmi
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

  unsigned int npar = prior.slot("npar");
  unsigned int nsub = dmi.size();

  Prior * p0 = new Prior (prior);
  Prior * lp = new Prior (lprior);
  Prior * sp = new Prior (sprior);

  std::vector<Design     *> ds(nsub);
  std::vector<Prior      *> ps(nsub);
  std::vector<Likelihood *> ls(nsub);
  std::vector<Theta      *> ts(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    S4 dmi_i = dmi[i];  // Must cast out first
    ds[i] = new Design (dmi_i);
    ps[i] = new Prior (prior);
    ls[i] = new Likelihood (dmi_i, ds[i], 3.0); // diff. RTs, Rs for diff. subjs
    ts[i] = new Theta (nmc, nchain, npar, thin, ps[i], ls[i]);
  }

  Phi     * phi = new Phi (nmc, nchain, npar, nsub, thin, p0, lp, sp, ts);
  Sampler * s0  = new Sampler (nchain, npar, gammamult, rp);
  // Rcout << "Start sampling: ";

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

  std::vector<std::string> pnames = prior.slot("pnames");
  std::vector<std::string> snames = dmi.attr("names");
  List individuals(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    S4 tmp ("posterior");
    tmp.slot("theta")            = ts[i]->m_theta;
    tmp.slot("summed_log_prior") = ts[i]->m_lp;
    tmp.slot("log_likelihoods")  = ts[i]->m_ll;
    tmp.slot("dmi")              = dmi[i];
    tmp.slot("prior")            = prior;
    tmp.slot("start")            = ts[i]->m_start_R;
    tmp.slot("npar")             = npar;
    tmp.slot("pnames")           = pnames;
    tmp.slot("nmc")              = nmc;
    tmp.slot("thin")             = thin;
    tmp.slot("nchain")           = nchain;

    individuals[i] = tmp;
  }

  // List phi_tmp = List::create(
  //   Named("location") = phi->m_phi0,
  //   Named("scale")    = phi->m_phi1);
  //
  // List ppprior = List::create(
  //   Named("location") = lprior,
  //   Named("scale")    = sprior);

  // List hyper = List::create(   // 16 elements
  //   Named("phi")                = phi_tmp,
  //   Named("h_summed_log_prior") = phi->m_hlp,
  //   Named("h_log_likelihoods")  = phi->m_hll,
  //   Named("pp.prior")           = ppprior,
  //   Named("start")              = phi->m_start_R,
  //   Named("npar")               = npar,
  //   Named("pnames")             = pnames,
  //   Named("rp")                 = rp,
  //   Named("nmc")                = nmc,
  //   Named("thin")               = thin,
  //   Named("nchain")             = nchain);
  // out.attr("hyper") = hyper;
  S4 out ("hyper");
  out.slot("phi_loc")          = phi->m_phi0;
  out.slot("phi_sca")          = phi->m_phi1;
  out.slot("summed_log_prior") = phi->m_hlp;
  out.slot("log_likelihoods")  = phi->m_hll;
  out.slot("prior_loc")        = lprior;
  out.slot("prior_sca")        = sprior;
  out.slot("start")            = phi->m_start_R;
  out.slot("npar")             = npar;
  out.slot("pnames")  = pnames;
  out.slot("nmc")     = nmc;
  out.slot("thin")    = thin;
  out.slot("nchain")  = nchain;
  out.slot("individuals") = individuals;
  out.slot("snames")      = snames;

  delete s0;
  delete phi;
  for (size_t i = 0; i < nsub; i++) delete ts[i];
  return out;
}

// [[Rcpp::export]]
S4 init_oldhier(S4 samples,
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
  S4 samples_in(clone(samples));

  // List hyper  = samples_in.attr("hyper");
  S4 lprior     = samples_in.slot("prior_loc");
  S4 sprior     = samples_in.slot("prior_sca");
  List subjects = samples_in.slot("individuals");
  unsigned int nchain = samples_in.slot("nchain");
  unsigned int npar   = samples_in.slot("npar");

  std::vector<std::string> snames = samples_in.slot("snames");
  unsigned int nsub   = snames.size();

  S4 subject0 = subjects[0]; // posterior class
  S4 prior    = subject0.slot("prior");

  Prior * p0 = new Prior (prior);
  Prior * lp = new Prior (lprior);
  Prior * sp = new Prior (sprior);

  std::vector<Design     *> ds(nsub);
  std::vector<Prior      *> ps(nsub);
  std::vector<Likelihood *> ls(nsub);
  std::vector<Theta      *> ts(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    S4 subjecti = subjects[i];
    S4 dmi_i  = subjecti.slot("dmi");  // This is a dmi and must cast out first
    S4 prior  = subjecti.slot("prior");

    ds[i] = new Design (dmi_i);
    ps[i] = new Prior (prior);
    ls[i] = new Likelihood (dmi_i, ds[i], 3.0); // diff. RTs, Rs for diff. subjs
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
        // arma::uvec pars = arma::randperm(npar);
        // for (size_t j=0; j<npar; j++) s0->crossover(pars[j], ts[k]);
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
  std::vector<std::string> pnames = prior.slot("pnames");
  List individuals(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    S4 subjecti = subjects[i];
    S4 dmi_i    = subjecti.slot("dmi");
    S4 tmp ("posterior");

    tmp.slot("theta")            = ts[i]->m_theta;
    tmp.slot("summed_log_prior") = ts[i]->m_lp;
    tmp.slot("log_likelihoods")  = ts[i]->m_ll;
    tmp.slot("dmi")              = dmi_i;
    tmp.slot("prior")            = prior;
    tmp.slot("start")            = ts[i]->m_start_R;
    tmp.slot("npar")             = npar;
    tmp.slot("pnames")           = pnames;
    tmp.slot("nmc")              = ts[i]->m_nmc;
    tmp.slot("thin")             = ts[i]->m_thin;
    tmp.slot("nchain")           = nchain;

    individuals[i] = tmp;
  }

  // List phi_tmp = List::create(
  //   Named("location") = phi->m_phi0,
  //   Named("scale")    = phi->m_phi1);
  //
  // List ppprior = List::create(
  //   Named("location") = lprior,
  //   Named("scale")    = sprior);

  // hyper = List::create(   // 16 elements
  //   Named("phi")                = phi_tmp,
  //   Named("h_summed_log_prior") = phi->m_hlp,
  //   Named("h_log_likelihoods")  = phi->m_hll,
  //   Named("pp.prior")           = ppprior,
  //   Named("start")              = phi->m_start_R,
  //   Named("npar")             = npar,
  //   Named("pnames")            = pnames,
  //   Named("rp")                 = rp,
  //   Named("nmc")                = phi->m_nmc,
  //   Named("thin")               = phi->m_thin,
  //   Named("nchain")           = nchain);
  //
  // out.attr("hyper") = hyper;

  S4 out ("hyper");
  out.slot("phi_loc")          = phi->m_phi0;
  out.slot("phi_sca")          = phi->m_phi1;
  out.slot("summed_log_prior") = phi->m_hlp;
  out.slot("log_likelihoods")  = phi->m_hll;
  out.slot("prior_loc")        = lprior;
  out.slot("prior_sca")        = sprior;
  out.slot("start")            = phi->m_start_R;
  out.slot("npar")             = npar;
  out.slot("pnames")  = pnames;
  out.slot("nmc")     = phi->m_nmc;
  out.slot("thin")    = phi->m_thin;
  out.slot("nchain")  = nchain;
  out.slot("individuals") = individuals;
  out.slot("snames")      = snames;


  delete s0;
  delete phi;
  for (size_t i = 0; i < nsub; i++) delete ts[i];
  return out;

}

// [[Rcpp::export]]
S4 init_oldhier_from_fixed_model(List samples,
                                  S4 lprior,
                                  S4 sprior,
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
  List subjects(clone(samples)); // from a fixed-effect model fit

  S4 subject0         = subjects[0];
  unsigned int nchain = subject0.slot("nchain");
  unsigned int npar   = subject0.slot("npar");
  std::vector<std::string> snames = subjects.attr("names");
  unsigned int nsub = snames.size();
  S4 prior = subject0.slot("prior");

  Prior * p0 = new Prior (prior);
  Prior * lp = new Prior (lprior);
  Prior * sp = new Prior (sprior);

  std::vector<Design     *> ds(nsub);
  std::vector<Prior      *> ps(nsub);
  std::vector<Likelihood *> ls(nsub);
  std::vector<Theta      *> ts(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    S4 subjecti = subjects[i];
    S4 dmi_i    = subjecti.slot("dmi");  // This is a dmi and must cast out first
    S4 prior    = subjecti.slot("prior");

    ds[i] = new Design (dmi_i);
    ps[i] = new Prior (prior);
    ls[i] = new Likelihood (dmi_i, ds[i], 3.0); // diff. RTs, Rs for diff. subjs
    ts[i] = new Theta (subjecti, nmc, thin, ps[i], ls[i], add);
  }

  // phi is created via the new Phi constructor
  Phi * phi = new Phi (nmc, nchain, npar, nsub, thin, p0, lp, sp, ts);
  Sampler * s0  = new Sampler (nchain, npar, gammamult, rp);
  Rcout << "Start sampling (from a fixed-effect model): ";

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
  Rcout << "\n";

  ////////////////////////////////////////////////////////////////////////
  std::vector<std::string> pnames = prior.slot("pnames");
  List individuals(nsub);

  for (size_t i = 0; i < nsub; i++)
  {
    S4 subjecti = subjects[i];
    S4 dmi_i    = subjecti.slot("dmi");
    S4 tmp ("posterior");

    tmp.slot("theta")            = ts[i]->m_theta;
    tmp.slot("summed_log_prior") = ts[i]->m_lp;
    tmp.slot("log_likelihoods")  = ts[i]->m_ll;
    tmp.slot("dmi")              = dmi_i;
    tmp.slot("prior")            = prior;
    tmp.slot("start")            = ts[i]->m_start_R;
    tmp.slot("npar")             = npar;
    tmp.slot("pnames")           = pnames;
    tmp.slot("nmc")              = ts[i]->m_nmc;
    tmp.slot("thin")             = ts[i]->m_thin;
    tmp.slot("nchain")           = nchain;

    individuals[i] = tmp;
  }
  S4 out ("hyper");
  out.slot("phi_loc")          = phi->m_phi0;
  out.slot("phi_sca")          = phi->m_phi1;
  out.slot("summed_log_prior") = phi->m_hlp;
  out.slot("log_likelihoods")  = phi->m_hll;
  out.slot("prior_loc")        = lprior;
  out.slot("prior_sca")        = sprior;
  out.slot("start")            = phi->m_start_R;
  out.slot("npar")             = npar;
  out.slot("pnames")  = pnames;
  out.slot("nmc")     = phi->m_nmc;
  out.slot("thin")    = phi->m_thin;
  out.slot("nchain")  = nchain;
  out.slot("individuals") = individuals;
  out.slot("snames")      = snames;

  delete s0;
  delete phi;
  for (size_t i = 0; i < nsub; i++) delete ts[i];
  return out;
}

// [[Rcpp::export]]
arma::uvec timesTwo(int x) {
  arma::uvec X = arma::randperm(x);
  return X;
}
