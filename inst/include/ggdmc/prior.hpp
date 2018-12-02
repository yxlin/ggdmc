#include <RcppArmadillo.h>

arma::vec dprior_(arma::vec pvec, arma::vec dists,
                  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
                  arma::uvec islog);
Rcpp::NumericVector dprior(Rcpp::NumericVector pvec, Rcpp::List prior);  


Rcpp::NumericVector rprior_vec(Rcpp::List prior);
Rcpp::NumericMatrix rprior_mat(Rcpp::List prior, unsigned int n) ;
arma::vec rprior_vec_(arma::vec dists, arma::vec p1, arma::vec p2,
                     arma::vec lower, arma::vec upper);
arma::mat rprior_mat_(unsigned int n, arma::vec dists, arma::vec p1,
                 arma::vec p2, arma::vec lower, arma::vec upper);

double sumlogpriorNV(arma::vec pvec, Rcpp::List prior);
double sumlogprior(arma::vec pvec, arma::vec dists, arma::vec p1, arma::vec p2,
                   arma::vec lower, arma::vec upper, arma::uvec islog);

Rcpp::List RestorePrior(Rcpp::List pprior, arma::vec p1, arma::vec p2, 
                        arma::vec lower, arma::vec upper, arma::uvec islog,
                        arma::vec dists, std::vector<std::string> pnames);

void GetPrior(Rcpp::List pprior, arma::vec& dist, arma::vec& p1,
              arma::vec& p2, arma::vec& lower, arma::vec& upper, arma::uvec& lg);

arma::vec UpdatePriors(arma::mat theta, arma::vec dists, arma::mat p1, 
                       arma::mat p2, arma::vec lower, arma::vec upper, 
                       arma::uvec islog);

double sumloghprior(arma::vec location, arma::vec scale, arma::vec ldists,
                    arma::vec sdists, arma::vec lp1, arma::vec sp1, 
                    arma::vec lp2, arma::vec sp2, arma::vec llower, 
                    arma::vec slower, arma::vec lupper, arma::vec supper, 
                    arma::uvec llog, arma::uvec slog);

double sumloghlike(arma::mat thetak, arma::vec dists,
                   arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
                   arma::uvec islog);
