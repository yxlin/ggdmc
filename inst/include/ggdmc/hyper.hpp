#include <RcppArmadillo.h>

arma::cube GetTheta0(Rcpp::List samples);
arma::cube GetUsethetas(arma::field<arma::mat> subusethetas);
arma::vec UpdatePriors(arma::mat theta, std::vector<std::string> dists,
  arma::mat p1, arma::mat p2, arma::vec lower, arma::vec upper,
  arma::uvec islog);

arma::mat UpdatePriors2(arma::cube theta, std::vector<std::string> dists,
  arma::mat p1, arma::mat p2, arma::vec lower, arma::vec upper,
  arma::uvec islog);

double sumloghprior(arma::vec location, arma::vec scale,
  std::vector<std::string> ldists, std::vector<std::string> sdists,
  arma::vec lp1, arma::vec sp1, arma::vec lp2, arma::vec sp2, arma::vec llower,
  arma::vec slower, arma::vec lupper, arma::vec supper, arma::uvec llog,
  arma::uvec slog);

double sumloghlike(arma::mat thetak, std::vector<std::string> dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog);

void Get_n(Rcpp::List samples, unsigned int &ncondition,
                 unsigned int& nallpar, unsigned int& nresponse);

double sumloghlike(arma::mat data, std::vector<std::string> dist, arma::vec p1,
                   arma::vec p2, arma::vec lower, arma::vec upper,
                   arma::uvec islog);

Rcpp::List DrawHPrior(Rcpp::List subject0, arma::mat usephi0, arma::mat usephi1,
                      bool debug);

