#include <ggdmc.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec dprior_(arma::vec pvec, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {
  unsigned int npar = dists.size();
  std::string dist1 ("tnorm");  // Available pdf' in DMC
  std::string dist2 ("beta_lu");
  std::string dist3 ("gamma_l");
  std::string dist4 ("lnorm_l");
  std::string dist5 ("constant");
  arma::vec out(npar); out.fill(NA_REAL);
  double x, l, u, tmp;

  for (size_t i = 0; i < npar; i++)
  {
    if ( dists[i].compare(dist1) == 0 ) {         // tnorm
      l = std::isnan(lower[i]) ? -INFINITY : lower[i];
      u = std::isnan(upper[i]) ? INFINITY : upper[i];
      x = pvec[i];
      tmp = dtn_scalar(x, p1[i], p2[i], l, u, islog[i]);
      out[i] = std::isnan(tmp) ? 1e-10 : tmp;
    } else if ( dists[i].compare(dist2) == 0) {  // beta_lu
      // Rcout << "beta_lu " << std::endl;
      l = std::isnan(lower[i]) ? 0 : lower[i];
      u = std::isnan(upper[i]) ? 1 : upper[i];
      x = (pvec[i] - l) / (u - l);
      tmp = islog[i] ? R::dbeta(x, p1[i], p2[i], islog[i]) - std::log(u - l) :
        R::dbeta(x, p1[i], p2[i], islog[i]) / (u - l);
      out[i] = std::isnan(tmp) ? 1e-10 : tmp;

    } else if ( dists[i].compare(dist3) == 0) {  // gamma_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      x = (std::isinf(l) || std::isinf(u)) ? pvec[i] : pvec[i] - l;
      out[i] = R::dgamma(x, p1[i], p2[i], islog[i]);
    } else if ( dists[i].compare(dist4) == 0) {  // lnorm_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      x = (std::isinf(l) || std::isinf(u)) ? pvec[i] : pvec[i] - l;
      out[i] = R::dlnorm(x, p1[i], p2[i], islog[i]);
    } else if (dists[i].compare(dist5) == 0) {  // constant
      out[i] = islog[i] ? 0 : 1;
    } else {
      Rcout << "Distribution type not yet supported" << "\n";
      out[i] = 1e-10;
    }
  }
  return out;
}

// [[Rcpp::export]]
double sumlogprior(arma::vec pvec, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {

  arma::vec den = dprior_(pvec, dists, p1, p2, lower, upper, islog);
  den.replace(arma::datum::inf, 1e-10);
  // den.replace(-arma::datum::inf, 1e-10);  // replace each -INFINITY with 1e-10
  // den.replace(arma::datum::nan, 1e-10);  // replace each nan with 1e-10
  double out = arma::accu(den);
  // if (std::isnan(out)) { out = -INFINITY; }
  return out;
}

//' Probability densities of prior distributions
//'
//' \code{sumlogpriorNV} calculate sum log-likelihood. \code{rprior} generates
//' random observation from prior distributions
//'
//' @param pvec parameter vector
//' @param prior prior distributions
//'
//' @return a vector
//' @export
// [[Rcpp::export]]
NumericVector dprior(NumericVector pvec, List prior) {
  List l1;   // a container to loop through inner list
  std::vector<std::string> pnames = pvec.names() ;
  NumericVector out = NumericVector(pvec.size());

  std::string distType1 ("tnorm");  // Available pdf' in DMC
  std::string distType2 ("beta_lu");
  std::string distType3 ("gamma_l");
  std::string distType4 ("lnorm_l");
  std::string distType5 ("constant");

  bool islog;
  double x, p1, p2, lower, upper, den;
  unsigned int npar = pvec.size();

  for (size_t i = 0; i < npar; i++) {
    l1 = prior[pnames[i]];
    std::string distName = l1.attr("dist");
    p1 = as<double>(l1[0]);  // mean; shape1; shape; meanlog
    p2 = as<double>(l1[1]);  // sd;   shape2; scale; sdlog

    // Do do.call
    if (distName.compare(distType1) == 0) {         // tnorm
      lower = as<double>(l1[2]);
      upper = as<double>(l1[3]);
      islog = as<bool>(l1[4]);
      x = pvec[i];
      den = dtn_scalar(x, p1, p2, lower, upper, islog);
    } else if (distName.compare(distType2) == 0) {  // beta_ul, shape1, shape2
      lower = as<double>(l1[2]);
      upper = as<double>(l1[3]);
      islog = as<bool>(l1[4]);
      x = (pvec[i] - lower) / (upper - lower);
      den = !islog ? R::dbeta(x, p1, p2, 0) / (upper - lower) :
        R::dbeta(x, p1, p2, 1) - std::log(upper - lower);
    } else if (distName.compare(distType3) == 0) {  // gamma_l, shape, scale
      lower = as<double>(l1[2]);
      islog = as<bool>(l1[3]);
      x = pvec[i] - lower;
      den = R::dgamma(x, p1, p2, islog);
    } else if (distName.compare(distType4) == 0) {  // lnorm_l, meanlow, sdlow
      lower = as<double>(l1[2]);
      islog = as<bool>(l1[3]);
      x = pvec[i] - lower;
      den = R::dlnorm(x, p1, p2, islog);
    } else if (distName.compare(distType5) == 0) {  // constant
      islog = Rcpp::as<bool>(l1[1]);
      den = islog ? 0 : 1;
    } else {
      Rcout << "Distribution type not yet supported" << "\n";
      den = 1e-10;
    }
    out[i] = den;
  }

  out.attr("names") = pnames;
  return out;
}

//' @rdname rprior
//' @export
// [[Rcpp::export]]
NumericVector rprior_scalar(List prior) {
  List l1;   // a list to hold each parameter's setting inside priorList
  unsigned int npar = prior.size();
  NumericVector out = Rcpp::NumericVector(npar);
  std::vector<std::string> pnames = prior.attr("names");

  std::string distType1 ("tnorm");
  std::string distType2 ("beta_lu");
  std::string distType3 ("gamma_l");
  std::string distType4 ("lnorm_l");
  std::string distType5 ("constant");

  double p1, p2, lower, upper;

  for (size_t i = 0; i < npar; i++) {
    l1 = prior[pnames[i]];
    std::string distName = l1.attr("dist");

    p1 = as<double>(l1[0]);  // parameter1: mean; shape1; shape; meanlog
    p2 = as<double>(l1[1]);  // parameter2: sd;   shape2; scale; sdlog
    lower = as<double>(l1[2]);
    upper = as<double>(l1[3]);

    if (distName.compare(distType1) == 0) {         // tnorm
      out[i] = rtn_scalar(p1, p2, lower, upper);
    } else if (distName.compare(distType2) == 0) {  // beta_ul
      out[i] = lower + R::rbeta(p1, p2) * (upper - lower);
    } else if (distName.compare(distType3) == 0) {  // gamma_l
      out[i] = lower + R::rgamma(p1, p2);
    } else if (distName.compare(distType4) == 0) {  // lnorm_l
      out[i] = lower + R::rlnorm(p1, p2);
    } else if (distName.compare(distType5) == 0) {  // constant
      out[i] = p1;
    } else {
      Rcout << "Distribution type not yet supported\n";
      out[i] = 1e-10;
    }
  }

  out.attr("names") = pnames;
  return out;
}

//' @rdname rprior
//' @export
// [[Rcpp::export]]
NumericMatrix rprior_mat(List prior, unsigned int n) {
  unsigned int nrow = n;
  unsigned int npar = prior.size() ;
  NumericMatrix m = na_matrix(nrow, npar) ;
  NumericVector a;

  for (size_t i = 0; i < nrow; i++) {
    a = rprior_scalar(prior);
    for (size_t j = 0; j < npar; j++) m(i, j) = a[j] ;
  }

  CharacterVector pnames = prior.names();
  colnames(m) = pnames;
  return m ;
}

// [[Rcpp::export]]
arma::vec rprior_vec(std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper) {
  unsigned int npar = dists.size();
  std::string dist1 ("tnorm"); // Available pdf's Andrew has implemented
  std::string dist2 ("beta_lu");
  std::string dist3 ("gamma_l");
  std::string dist4 ("lnorm_l");
  std::string dist5 ("constant");
  arma::vec out = arma::vec(npar).fill(NA_REAL);
  double l, u;

  // [p1 p2]: [mean sd]; [shape1 shape2]; [shape scale]; [meanlog sdlog]
  for (size_t i = 0; i < npar;  i++) {
    if ( dists[i].compare(dist1) == 0 ) {         // tnorm
      l = std::isnan(lower[i]) ? -INFINITY : lower[i];
      u = std::isnan(upper[i]) ?  INFINITY : upper[i];
      out[i] = rtn_scalar(p1[i], p2[i], l, u);
    } else if ( dists[i].compare(dist2) == 0) {  // beta_ul
      l = std::isnan(lower[i]) ? 0 : lower[i];
      u = std::isnan(upper[i]) ? 1 : upper[i];
      out[i] = l + R::rbeta(p1[i], p2[i]) * (u - l);
    } else if ( dists[i].compare(dist3) == 0) {  // gamma_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      out[i] = R::rgamma(p1[i], p2[i]) + l;
    } else if ( dists[i].compare(dist4) == 0 ) {  // lnorm_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      out[i] = R::rlnorm(p1[i], p2[i]) + l;
    } else {  // constant
      out[i] = p1[i];
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat rprior(unsigned int n, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper) {
  unsigned int npar = dists.size();
  arma::mat m = arma::mat(npar, n).zeros();
  for (size_t i = 0; i < n; i++) m.col(i) = rprior_vec(dists, p1, p2, lower, upper);
  return m;
}

//' @rdname dprior
//' @export
// [[Rcpp::export]]
double sumlogpriorNV(arma::vec pvec, List prior) {
  NumericVector pvecNV  = as<NumericVector>(wrap(pvec)) ;
  std::vector<std::string> pnames = prior.names() ;
  pvecNV.names() = pnames;
  return sum(dprior(pvecNV, prior)); // sum is in Rcpp
}
