#include <RcppArmadillo.h>

bool checkDDM(std::vector<double> pVec) ;

Rcpp::NumericMatrix getAccumulatorMatrix(Rcpp::NumericVector pVec,
  std::string cell, Rcpp::NumericVector model, bool n1order);

arma::vec density_rd(arma::vec pVec,
  std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, double precision);

arma::vec density_norm(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift);

arma::vec density_norm_pda(arma::vec pVec, std::vector<std::string> pnames,
                           arma::vec allpar, std::vector<std::string> parnames,
                           arma::ucube model,
                           std::string type,
                           std::vector<std::string> dim1,
                           std::vector<std::string> dim2,
                           std::vector<std::string> dim3,
                           arma::umat n1idx, arma::uvec ise, arma::umat cellidx,
                           arma::vec RT, arma::uvec matchcell, arma::uvec isr1,
                           unsigned int npda, double bw, bool debug);

arma::vec density_plba1(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int nsim, double bw,
  unsigned int ncore, bool debug);


arma::vec density_plba1_gpu(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int nsim,
  double bw, unsigned int ncore, unsigned int gpuid, unsigned int nthread,
  bool debug);

arma::vec density_cnorm_pda(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  unsigned int nsim, double bw, bool debug);

double sumloglike(arma::vec pvec, std::vector<std::string> pnames,
                  arma::vec allpar, std::vector<std::string> parnames,
                  arma::ucube model, std::string type,
                  std::vector<std::string> dim1,
                  std::vector<std::string> dim2,
                  std::vector<std::string> dim3,
                  arma::umat n1idx, arma::uvec ise, arma::umat cellidx,
                  arma::vec RT, arma::uvec matchcell, arma::uvec isr1, bool posdrift,
                  unsigned int nsim, double bw, unsigned int ncore,
                  unsigned int gpuid, bool debug);

