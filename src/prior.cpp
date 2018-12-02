#include <ggdmc.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec dprior_(arma::vec pvec, arma::vec dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog) {
  unsigned int npar = dists.size();
  // std::string dist1 ("tnorm");  
  // std::string dist2 ("beta_lu");
  // std::string dist3 ("gamma_l");
  // std::string dist4 ("lnorm_l");
  // std::string dist5 ("unif_");
  // std::string dist6 ("constant");
  // std::string dist7 ("tnorm2");
  arma::vec out(npar); out.fill(NA_REAL);
  double x, l, u, tmp;

  for (size_t i = 0; i < npar; i++)
  {
    if ( std::isnan(p1(i)) ) {
      // NA go here; NA will be converted to 0 (unsigned int type)
      out[i] = 0;
    } else if ( dists(i) == 1 ) {         // tnorm
      l = std::isnan(lower[i]) ? R_NegInf : lower[i];
      u = std::isnan(upper[i]) ? R_PosInf : upper[i];
      x = pvec[i];
      tmp = dtn_scalar(x, p1[i], p2[i], l, u, islog[i]);
      out[i] = std::isnan(tmp) ? -23.02585 : tmp;
    } else if ( dists(i) == 2 ) {  // beta_lu
      // Rcout << "beta_lu " << std::endl;
      l = std::isnan(lower[i]) ? 0 : lower[i];
      u = std::isnan(upper[i]) ? 1 : upper[i];
      x = (pvec[i] - l) / (u - l);
      tmp = islog[i] ? R::dbeta(x, p1[i], p2[i], islog[i]) - std::log(u - l) :
        R::dbeta(x, p1[i], p2[i], islog[i]) / (u - l);
      out[i] = std::isnan(tmp) ? -23.02585 : tmp;

    } else if ( dists(i) == 3 ) {  // gamma_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      x = (std::isinf(l) || std::isinf(u)) ? pvec[i] : pvec[i] - l;
      out[i] = R::dgamma(x, p1[i], p2[i], islog[i]);
    } else if (  dists(i) == 4 ) {  // lnorm_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      x = (std::isinf(l) || std::isinf(u)) ? pvec[i] : pvec[i] - l;
      out[i] = R::dlnorm(x, p1[i], p2[i], islog[i]);
    } else if ( dists(i) == 5 ) {  // unif_

      out[i] = R::dunif(pvec[i], p1[i], p2[i], islog[i]);
      
    } else if ( dists(i) == 6 ) { // constant
      // x = pvec[i];
      // tmp = dtn_scalar(x, p1[i], 1e-10, l, u, islog[i]);
      // 
      // out[i] = std::isnan(tmp) ? -23.02585 : tmp;
      // tmp = (pvec[i] == p1[i]) ? 1 : 1e-10;
      // out[i] = islog[i] ? 0 : 1;
      out[i] = 0;
    } else if ( dists(i) == 7) { // tnorm mean precision
      
      l = std::isnan(lower[i]) ? R_NegInf : lower[i];
      u = std::isnan(upper[i]) ? R_PosInf : upper[i];
      x = pvec[i];
      tmp = dtn_scalar2(x, p1[i], p2[i], l, u, islog[i]);
      out[i] = std::isnan(tmp) ? -23.02585 : tmp;

    } else {
      Rcout << "Distribution type not yet supported" << "\n";
      out[i] = R_NegInf;
    }
  }
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
  List L0;   // a container to loop through inner list
  std::vector<std::string> pnames = pvec.names() ;
  unsigned int npar = pvec.size();
  NumericVector out(npar);
  
  // std::string distType0 (NA);  
  // std::string distType1 ("tnorm");  
  // std::string distType2 ("beta_lu");
  // std::string distType3 ("gamma_l");
  // std::string distType4 ("lnorm_l");
  // std::string distType5 ("unif_");
  // std::string distType6 ("constant");
  
  bool lg;
  double x, p1, p2, lower, upper, den;
  
  for (size_t i = 0; i < npar; i++) {
    L0 = prior[pnames[i]];
    unsigned int distName = L0.attr("dist");
    p1 = L0[0];  // mean; shape1; shape; meanlog
    p2 = L0[1];  // sd;   shape2; scale; sdlog
    lower = L0[2];
    upper = L0[3];
    lg    = L0[4];
    
    if (distName == 0) {
      // Rcout << "Found NA\n";
      den = lg ? -arma::datum::inf : 1e-10;
    } else if (distName == 1) {         // tnorm
      x = pvec[i];
      den = dtn_scalar(x, p1, p2, lower, upper, lg);
    } else if (distName == 2) {  // beta_ul, shape1, shape2
      x = (pvec[i] - lower) / (upper - lower);
      den = !lg ? R::dbeta(x, p1, p2, 0) / (upper - lower) :
        R::dbeta(x, p1, p2, 1) - std::log(upper - lower);
    } else if (distName == 3) {  // gamma_l, shape, scale
      x = pvec[i] - lower;
      den = R::dgamma(x, p1, p2, lg);
    } else if (distName == 4) {  // lnorm_l, meanlow, sdlow
      
      x = pvec[i] - lower;
      den = R::dlnorm(x, p1, p2, lg);
    } else if (distName == 5) {  // unif_
      den = R::dunif(pvec[i], p1, p2, lg);
    } else if (distName == 6) { // constant
      
      double tmp = (pvec[i] == p1) ? 1 : 1e-10;
      den = lg ? std::log(tmp) : tmp;
    } else if (distName == 7) {
      x = pvec[i];
      den = dtn_scalar2(x, p1, p2, lower, upper, lg);
    } else {
      Rcout << "Distribution type not yet supported\n";
      den = R_NegInf;
    }
    out[i] = den;
  }
  
  out.attr("names") = pnames;
  return out;
}

//' @rdname rprior
//' @export
// [[Rcpp::export]]
NumericVector rprior_vec(List prior) {
  unsigned int npar = prior.size();
  NumericVector out = NumericVector(npar);
  std::vector<std::string> pnames = prior.attr("names");

  double dist;
  double p1, p2, lower, upper;
  List L0;   
  
  for (size_t i = 0; i < npar; i++) {
    L0 = prior[pnames[i]];
    dist = L0.attr("dist");
    p1 = L0[0];  // parameter1: mean; shape1; shape; meanlog; min
    p2 = L0[1];  // parameter2: sd;   shape2; scale; sdlog;   max
    lower = L0[2];
    upper = L0[3];

    if ( std::isnan(dist) ) {
      // Rcout << "NA found\n";
      out[i] = NA_REAL;
    } else if (dist == 1.0) {         // tnorm
      out[i] = rtn_scalar(p1, p2, lower, upper);
      
    } else if (dist == 2.0) {  // beta_ul
      out[i] = lower + R::rbeta(p1, p2) * (upper - lower);
      
    } else if (dist == 3.0) {  // gamma_l; 
      out[i] = lower + R::rgamma(p1, p2);
      
    } else if (dist == 4.0) {  // lnorm_l
      out[i] = lower + R::rlnorm(p1, p2);
      
    } else if (dist == 5.0) {  // unif_
      out[i] = R::runif(p1, p2);
      
    } else if (dist == 6.0) { // constant
      out[i] = p1;
      
    } else if (dist ==7.0) { // tnorm2
      out[i] = rtn_scalar2(p1, p2, lower, upper);
      
    } else {
      Rcout << "Distribution type not supported\n";
      out[i] = NA_REAL;
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
    a = rprior_vec(prior);
    for (size_t j = 0; j < npar; j++) m(i, j) = a[j] ;
  }
  
  CharacterVector pnames = prior.names();
  colnames(m) = pnames; // colnames is Rcpp function?
  return m ;
}

// [[Rcpp::export]]
arma::vec rprior_vec_(arma::vec dists, arma::vec p1, arma::vec p2, 
                      arma::vec lower, arma::vec upper) {
  // replace DMC modified r-function; used in initialise.cpp internally
  unsigned int npar = dists.size();
  arma::vec out = arma::vec(npar).fill(NA_REAL);
  double l, u;

  // [p1 p2]: [mean sd]; [shape1 shape2]; [shape scale]; [meanlog sdlog]
  for (size_t i = 0; i < npar;  i++) {
    if ( std::isnan(dists(i)) ) {
      out[i] = NA_REAL;
    } else if ( dists(i) == 1 ) {         // tnorm
      l = std::isnan(lower[i]) ? R_NegInf : lower[i];
      u = std::isnan(upper[i]) ? R_PosInf : upper[i];
      out[i] = rtn_scalar(p1[i], p2[i], l, u);
    } else if ( dists(i) == 2 ) {  // beta_ul
      l = std::isnan(lower[i]) ? 0 : lower[i];
      u = std::isnan(upper[i]) ? 1 : upper[i];
      out[i] = l + R::rbeta(p1[i], p2[i]) * (u - l);
    } else if ( dists(i) == 3 ) {  // gamma_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      out[i] = R::rgamma(p1[i], p2[i]) + l;
    } else if ( dists(i) == 4 ) {  // lnorm_l
      l = std::isnan(lower[i]) ? 0 : lower[i];
      out[i] = R::rlnorm(p1[i], p2[i]) + l;
    } else if ( dists(i) == 5 ) {
      out[i] = R::runif(p1[i], p2[i]);
    } else if ( dists(i) == 6 ){  // constant
      out[i] = p1[i];
    } else if ( dists(i) == 7 ) { // tnorm2
      l = std::isnan(lower[i]) ? R_NegInf : lower[i];
      u = std::isnan(upper[i]) ? R_PosInf : upper[i];
      out[i] = rtn_scalar2(p1[i], p2[i], l, u);
      
    } else {
      Rcout << "Distribution type not supported\n";
      out[i] = NA_REAL;
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::mat rprior_mat_(unsigned int n, arma::vec dists, arma::vec p1, 
                      arma::vec p2, arma::vec lower, arma::vec upper) {
  unsigned int npar = dists.size();
  arma::mat m = arma::mat(npar, n).fill(NA_REAL);
  for (size_t i = 0; i < n; i++) {
    m.col(i) = rprior_vec_(dists, p1, p2, lower, upper); 
  }
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

// [[Rcpp::export]]
double sumlogprior(arma::vec pvec, arma::vec dists,
                   arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
                   arma::uvec lg) {
  
  arma::vec den = dprior_(pvec, dists, p1, p2, lower, upper, lg);
  den.replace(R_PosInf, 1e-10);
  // den.replace(-arma::datum::inf, 1e-10);  // replace each -INFINITY with 1e-10
  // den.replace(arma::datum::nan, 1e-10);  // replace each nan with 1e-10
  double out = arma::accu(den);
  if (std::isnan(out)) { out = R_NegInf; }
  return out;
}

// [[Rcpp::export]]
List RestorePrior(List pprior, arma::vec p1, arma::vec p2, arma::vec lower, 
                  arma::vec upper, arma::uvec lg, arma::vec dists, 
                  std::vector<std::string> pnames) {
  
  unsigned int npar = p1.n_elem;
  std::vector< Rcpp::List > v(npar);

  List L0, L1;   
  for(size_t i = 0; i < npar; i++) {
    L0 = pprior[pnames[i]];
    if (std::isnan(p1(i))) {
      L1 = List::create(
        _["p1"] = L0[0],
        _["p2"] = L0[1],
        _["lower"] = L0[2],
        _["upper"] = L0[3],
        _["lg"] = L0[4]
      );
    } else {
      L1 = List::create(
        _["p1"] = p1(i),
        _["p2"] = p2(i),
        _["lower"] = lower(i),
        _["upper"] = upper(i),
        _["lg"] = lg(i)
      );
    }
    L1.attr("dist") = dists[i];
    L1.attr("untrans") = "identity";
    v[i] = L1;
  }
  List out = wrap(v);
  out.attr("names") = pnames;
  out.attr("class") = "prior";
  return out;
}

void GetPrior(List prior, arma::vec& dist, arma::vec& p1, arma::vec& p2, 
              arma::vec& lower, arma::vec& upper, arma::uvec& lg) {
  
  std::vector<std::string> pnames = prior.attr("names");
  unsigned int npar = pnames.size();
  List L1;
  dist.set_size(npar);
  p1.set_size(npar);
  p2.set_size(npar);
  lower.set_size(npar);
  upper.set_size(npar);
  lg.set_size(npar);
  double disti;
  
  for (size_t i = 0; i < npar; i++) {
    L1 = prior[pnames[i]];
    disti   = L1.attr("dist");
    dist(i) = disti;
    p1(i)   = L1[0];
    p2(i)   = L1[1];
    lower(i) = L1[2];
    upper(i) = L1[3];
    lg(i)    = L1[4];
  }
}

arma::vec UpdatePriors(arma::mat theta, arma::vec dists,
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


// [[Rcpp::export]]
double sumloghprior(arma::vec location, arma::vec scale, arma::vec ldists, 
                    arma::vec sdists, arma::vec lp1, arma::vec sp1, 
                    arma::vec lp2, arma::vec sp2, arma::vec llower, 
                    arma::vec slower, arma::vec lupper, arma::vec supper, 
                    arma::uvec llog, arma::uvec slog) {
  
  return sumlogprior(location, ldists, lp1, lp2, llower, lupper, llog) +
    sumlogprior(scale, sdists, sp1, sp2, slower, supper, slog);
}

// [[Rcpp::export]]
double sumloghlike(arma::mat thetak, arma::vec dists,
                   arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
                   arma::uvec islog) {
  double out = 0; // thetak: nsub x npar
  for(size_t i = 0; i < thetak.n_rows; i++) {
    out += sumlogprior(arma::trans(thetak.row(i)), dists, p1, p2, lower, upper,
                       islog);
  }
  
  // if ( std::isinf(out) ) out = 1e-10;
  return out;
}
