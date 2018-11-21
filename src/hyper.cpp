#include <ggdmc.hpp>
using namespace Rcpp;

arma::vec UpdatePriors(arma::mat theta, std::vector<std::string> dists,
  arma::mat p1, arma::mat p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {
  // theta = nchain x npar

  unsigned int nchain = theta.n_rows;
  arma::vec out(nchain);
  for (size_t i = 0; i < nchain; i++) {
    out(i) = sumlogprior(arma::trans(theta.row(i)), dists,
      arma::trans(p1.row(i)), arma::trans(p2.row(i)), lower, upper, islog);
  }

  return out ;
}


arma::mat UpdatePriors2(arma::cube theta, std::vector<std::string> dists,
  arma::mat p1, arma::mat p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {

  unsigned int nchain = theta.n_slices;
  unsigned int nsub   = theta.n_rows;
  arma::vec pvec, location, scale;
  // arma::mat out(nsub, nchain);
  arma::mat out(nchain, nsub);
  // p1 and p2 are usephi = nchain x npar
  // arma::uvec rchains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));
  // unsigned int k0;

  for (size_t i = 0; i < nchain; i++) {
    // k0 = rchains(i);
    location = arma::trans(p1.row(i));
    scale    = arma::trans(p2.row(i));

    for (size_t j = 0; j < nsub; j++) {
      pvec = arma::trans(theta.slice(i).row(j));
      out(i, j) = sumlogprior(pvec, dists, location, scale, lower, upper, islog);
    }
  }

  return out;
}

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
double sumloghprior(arma::vec location, arma::vec scale,
  std::vector<std::string> ldists, std::vector<std::string> sdists,
  arma::vec lp1, arma::vec sp1, arma::vec lp2, arma::vec sp2, arma::vec llower,
  arma::vec slower, arma::vec lupper, arma::vec supper, arma::uvec llog,
  arma::uvec slog) {
  return sumlogprior(location, ldists, lp1, lp2, llower, lupper, llog) +
  sumlogprior(scale, sdists, sp1, sp2, slower, supper, slog);
}

// [[Rcpp::export]]
double sumloghlike(arma::mat thetak, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {
  double out = 0; // thetak: nsub x npar
  for(size_t i = 0; i < thetak.n_rows; i++) {
    out += sumlogprior(arma::trans(thetak.row(i)), dists, p1, p2, lower, upper,
      islog);
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



