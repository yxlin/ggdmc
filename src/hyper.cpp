#include <ggdmc.hpp>
using namespace Rcpp;



//' Extract Start Posterior Sample
//'
//' Extract the theta's of the first MCMC iteration across chains and
//' participants. Note that the ps array in DMC is a nchain x nsubject x
//' nparameter array. Armadillo operates on slice (the third dimension), so
//' chain dimension has to be on slice.
//'
//' @param samples posterior samples
//' @return a nsubject x npar x nchain array
//' @export
// [[Rcpp::export]]
arma::cube GetTheta0(List samples) {
  List samples0 = samples[0];
  unsigned int nsub   = samples.size();
  unsigned int npar   = samples0["n.pars"];
  unsigned int nchain = samples0["n.chains"];
  unsigned int start  = samples0["start"];  // Note start is an R index
  arma::cube out(nsub, npar, nchain);

  for (size_t i = 0; i < nsub; i ++) {
    List subject     = samples[i];
    arma::cube theta = subject["theta"];    // nchain x npar x nmc
    for (size_t j = 0; j < nchain; j++) {
      out.slice(j).row(i) = theta.slice(start - 1).row(j);
    }
  }
  return out; // nsub x npar x nchain
}

arma::cube GetTheta0_new(List samples) {
  List samples0 = samples[0];
  unsigned int nsub   = samples.size();
  unsigned int npar   = samples0["n.pars"];
  unsigned int nchain = samples0["n.chains"];
  unsigned int start  = samples0["start"];  // Note start is an R index
  arma::cube out(nchain, npar, nsub);

  for (size_t i = 0; i < nsub; i ++) {
    List subject     = samples[i];
    arma::cube theta = subject["theta"];    // nchain x npar x nmc
    for (size_t j = 0; j < nchain; j++) {
      out.slice(j).row(i) = theta.slice(start - 1).row(j);
    }
  }
  return out; // nsub x npar x nchain
}

arma::cube GetUsethetas(arma::field<arma::mat> usethetas) {
  arma::mat usetheta0 = usethetas(0);   // nchain x npar
  unsigned int nsub = usethetas.n_elem;
  unsigned int npar = usetheta0.n_cols;
  unsigned int nchain = usetheta0.n_rows;
  arma::cube out(nsub, npar, nchain);

  for (size_t i = 0; i < nsub; i++) {
    for (size_t j = 0; j < nchain; j++) {
      out.slice(j).row(i) = usethetas(i).row(j);
    }
  }
  return out;
}

// [[Rcpp::export]]
void StartIteration(List samples) {
  List subjecti;
  unsigned int start_C;
  unsigned int nsub = samples.size();
  arma::mat theta;
  arma::vec logprior, loglike;

  for (size_t i = 0; i < nsub; i++) {
    subjecti = samples[i];
    arma::cube thetai   = subjecti["theta"] ; // nchain x npar x nmc
    arma::mat logpriori = subjecti["summed_log_prior"] ; // nmc x nchain
    arma::mat loglikei  = subjecti["log_likelihoods"] ;  // nmc x nchain
    unsigned int start_R  = subjecti["start"] ; // R index
    start_C  = start_R - 1 ;        // C index
    theta    = thetai.slice(start_C) ; // nchain x npar
    logprior = vectorise(logpriori.row(start_C)) ;   // nchain x 1
    loglike  = vectorise(loglikei.row(start_C)) ;   // nchain x 1
  }
}

void Get_n(List samples, unsigned int& nallpar, unsigned int& ncondition,
            unsigned int& nresponse) {
  List subject0 = samples[0];
  List pprior   = subject0["p.prior"];
  List data     = subject0["data"];
  NumericVector modelAttr = data.attr("model");

  arma::vec allpar = modelAttr.attr("all.par");
  List modelDim    = modelAttr.attr("dimnames");
  std::vector<std::string> dim1 = modelDim[0] ; // row; s1.r1, etc
  std::vector<std::string> dim3 = modelDim[2] ; // row; s1.r1, etc
  ncondition = dim1.size();
  nallpar    = allpar.n_elem;
  nresponse  = dim3.size();

}



