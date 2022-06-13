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

//' Calculate likelihoods
//'
//' These function calculate likelihoods. \code{likelihood_rd} implements
//' the equations in Voss, Rothermund, and Voss (2004). These equations
//' calculate diffusion decision model (Ratcliff & Mckoon, 2008). Specifically,
//' this function implements Voss, Rothermund, and Voss's (2004) equations A1
//' to A4 (page 1217) in C++.
//'
//' @param pvector a parameter vector
//' @param data data model instance
//' @param min_lik minimal likelihood.
//' @param precision a tuning parameter for the precision of DDM likelihood.
//' The larger the value is, the more precise the likelihood is and the slower
//' the computation would be.
//' @return a vector
//' @references Voss, A., Rothermund, K., & Voss, J. (2004).  Interpreting the
//' parameters of the diffusion model: An empirical validation.
//' \emph{Memory & Cognition}, \bold{32(7)}, 1206-1220. \cr\cr
//' Ratcliff, R. (1978). A theory of memory retrival. \emph{Psychological
//' Review}, \bold{85}, 238-255.
//'
//' @examples
//' model <- BuildModel(
//' p.map     = list(A = "1", B = "1", t0 = "1", mean_v = "M", sd_v = "1",
//'             st0 = "1"),
//' match.map = list(M = list(s1 = 1, s2 = 2)),
//' factors   = list(S = c("s1", "s2")),
//' constants = c(st0 = 0, sd_v = 1),
//' responses = c("r1", "r2"),
//' type      = "norm")
//'
//' p.vector <- c(A = .25, B = .35,  t0 = .2, mean_v.true = 1, mean_v.false = .25)
//' dat <- simulate(model, 1e3,  ps = p.vector)
//' dmi <- BuildDMI(dat, model)
//' den <- likelihood(p.vector, dmi)
//'
//' model <- BuildModel(
//' p.map     = list(a = "1", v = "1", z = "1", d = "1", t0 = "1", sv = "1",
//'             sz = "1", st0 = "1"),
//' constants = c(st0 = 0, d = 0),
//' match.map = list(M = list(s1 = "r1", s2 = "r2")),
//' factors   = list(S = c("s1", "s2")),
//' responses = c("r1", "r2"),
//' type      = "rd")
//'
//' p.vector <- c(a = 1, v = 1, z = 0.5, sz = 0.25, sv = 0.2, t0 = .15)
//' dat <- simulate(model, 1e2, ps = p.vector)
//' dmi <- BuildDMI(dat, model)
//' den <- likelihood (p.vector, dmi)
//'
//' @export
// [[Rcpp::export]]
std::vector<double> likelihood (arma::vec pvector, S4 data,
                                double min_lik=1e-10, double precision=3.0)
// used only by R
{
  Design     * obj0 = new Design(data);
  Likelihood * obj1 = new Likelihood(data, obj0, precision);
  arma::vec tmp = obj1->likelihood(pvector);

  std::vector<double> out(obj0->m_nRT);
  for(size_t i=0; i<obj0->m_nRT; i++)
  {
    out[i] = R::fmax2(tmp[i], min_lik);
  }

  delete obj1;
  return out;
}

// [[Rcpp::export]]
arma::mat p_df(arma::vec pvector, std::string cell, std::string mtype,

               std::vector<std::string> pnames,
               std::vector<std::string> parnames,
               std::vector<std::string> dim0,
               std::vector<std::string> dim1,
               std::vector<std::string> dim2,

               std::vector<double> allpar,
               arma::ucube model,

               arma::uvec isr1,
               arma::umat n1idx,
               bool n1order)
// Used only in random.R
{
  Design     * obj0 = new Design(pnames, parnames, dim0, dim1, dim2,
                                 allpar, model);
  Likelihood * obj1 = new Likelihood(mtype, isr1, n1idx, n1order, obj0);

  arma::mat pmat = obj1->get_pmat(pvector, cell); //

  delete obj1; // obj0 is freed in obj1;
  return pmat;

}

// [[Rcpp::export]]
arma::vec ac_(arma::vec x, unsigned int nlag) {

  unsigned int n = x.n_elem;
  unsigned int nm1 = n - 1;
  arma::vec out(nlag);
  arma::mat tmp0 = arma::cor(x, x);
  out(0) = arma::as_scalar(tmp0);
  arma::vec tmp1, tmp2;
  for (size_t i = 1; i < nlag; i++)
  {
    tmp1 = arma::shift(x, (int)i);
    tmp0 = arma::cor(x.rows(i, nm1), tmp1.rows(i, nm1)); // pairwise.complete.obs
    out(i) = arma::as_scalar(tmp0);
  }
  return out;
}

// [[Rcpp::export]]
arma::cube trial_loglik(S4 samples, unsigned int thin_pointwise)
{
  S4 dmi = samples.slot("dmi");
  unsigned int nobs, n, j;
  arma::vec tmp, nmc_thin;

  unsigned int nchain = samples.slot("nchain");
  unsigned int pnmc   = samples.slot("nmc");
  arma::cube theta    = samples.slot("theta"); // npar x nchain x nmc
  Design * d0 = new Design (dmi);
  nobs        = d0->m_nRT;;
  nmc_thin = arma::regspace(thin_pointwise, thin_pointwise, pnmc) - 1;
  n        = nmc_thin.n_elem;

  arma::cube out(nobs, nchain, n); out.fill(R_NegInf);

  Rcout << "Processing chains: ";
  for (size_t k=0; k<nchain; k++)
  {
    for (size_t i=0; i<n; i++)
    { // sample at certain regular iterations to calculate likelihoods
      j    = nmc_thin[i];
      tmp  = likelihood(theta.slice(j).col(k), dmi);
      out.slice(i).col(k) = arma::log(tmp);
    }
    Rcout << ".";
  }

  Rcout << "\n";
  delete d0;
  return out;
}

// List trial_loglik_hier(List samples, unsigned int thin_pointwise)
// {
//   List samples_in(clone(samples));
//
//   List hyper  = samples_in.attr("hyper");
//   List pprior = hyper ["pp.prior"];
//   List lprior = pprior["location"];
//   List sprior = pprior["scale"];
//
//   List subject0 = samples_in[0];
//   List prior    = subject0["p.prior"];
//
//   unsigned int npar   = prior.size();
//   unsigned int nsub   = samples.size();
//   unsigned int nchain = hyper["n.chains"];
//
//   std::vector<Design     *> ds(nsub);
//   std::vector<Likelihood *> ls(nsub);
//
//   for (size_t i = 0; i < nsub; i++)
//   {
//     List subjecti = samples_in[i];
//     List datai    = subjecti["data"];  // Must cast out first
//
//     ds[i] = new Design (datai);
//     ls[i] = new Likelihood (datai, ds[i], 3.0); // diff. RTs, Rs for diff. subjs
//   }
//
//
//   Rcout << "Start sampling: ";
//
//
// }
