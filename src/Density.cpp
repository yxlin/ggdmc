#include <ggdmc.hpp>

using namespace Rcpp;

//' Calculate log likelihoods
//'
//' These function calculate log likelihoods. \code{likelihood_rd} implements
//' the equations in Voss, Rothermund, and Voss (2004). These equations
//' calculate diffusion decision model (Ratcliff & Mckoon, 2008). Specifically,
//' this function implements Voss, Rothermund, and Voss's (2004) equations A1
//' to A4 (page 1217) in C++.
//'
//' @param pvector a parameter vector
//' @param data data model instance
//' @param min_lik minimal likelihood.
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
std::vector<double> likelihood (arma::vec pvector, List data,
                                double min_lik=1e-10)
// used only by R
{
  Design     * obj0 = new Design(data);
  Likelihood * obj1 = new Likelihood(data, obj0);
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

