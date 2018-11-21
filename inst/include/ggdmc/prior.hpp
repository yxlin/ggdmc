#include <RcppArmadillo.h>

arma::vec dprior(arma::vec pvec,
  std::vector<std::string> dists, arma::vec p1, arma::vec p2,
  arma::vec lower, arma::vec upper, arma::uvec islog);

Rcpp::NumericVector dpriorNV(Rcpp::NumericVector pvec, Rcpp::List prior);

double summedlogpriorNV(arma::vec pvec, Rcpp::List prior);

double sumlogprior(arma::vec pvec,
  std::vector<std::string> dists, arma::vec p1, arma::vec p2,
  arma::vec lower, arma::vec upper, arma::uvec islog);

Rcpp::NumericVector rprior_scalar(Rcpp::List pPrior);

Rcpp::NumericMatrix rprior_mat(Rcpp::List pPrior, unsigned int n) ;

arma::vec rprior_vec(std::vector<std::string> dists, arma::vec p1, arma::vec p2,
                     arma::vec lower, arma::vec upper);

arma::mat rprior(unsigned int n, std::vector<std::string> dists, arma::vec p1,
                 arma::vec p2, arma::vec lower, arma::vec upper);


