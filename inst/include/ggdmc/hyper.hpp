#include <RcppArmadillo.h>

arma::cube GetTheta0(Rcpp::List samples);
arma::cube GetUsethetas(arma::field<arma::mat> subusethetas);


void Get_n(Rcpp::List samples, unsigned int &ncondition, unsigned int& nallpar,
           unsigned int& nresponse);

Rcpp::List DrawHPrior(Rcpp::List subject0, arma::mat usephi0, arma::mat usephi1,
                      bool debug);

