#include <RcppArmadillo.h>

void InitializeOneSubject(Rcpp::List samples, arma::umat& rj);
void InitializeSubjectRJ(Rcpp::List samples, arma::field<arma::umat>& rj);

void TransformSubjects(Rcpp::List samples, arma::field<arma::cube>& thetas,
  arma::field<arma::mat>& usethetas, arma::field<arma::mat>& logpriors,
  arma::field<arma::vec>& uselogpriors, arma::field<arma::mat>& loglikes,
  arma::field<arma::vec>& useloglikes, arma::uvec& store_i,
  std::vector<std::string>& types, arma::field<arma::vec>& allpars,
  arma::field<arma::umat>& n1idxes, arma::field<arma::uvec>& matchcells,
  arma::field<arma::uvec>& emptycells, arma::field<arma::umat>& cellidxes,
  arma::field<std::vector<std::string>>& parnames,
  arma::field<std::vector<std::string>>& dim1s,
  arma::field<std::vector<std::string>>& dim2s,
  arma::field<std::vector<std::string>>& dim3s,
  arma::field<arma::uvec>& isr1s, arma::uvec& posdrift,
  arma::field<arma::ucube>& models,
  arma::uvec& npdas, arma::vec& bws, arma::uvec& gpuids,
  arma::field<arma::vec>& RTs);

void GetPrior(Rcpp::List pprior, std::vector<std::string>& dists, arma::vec& p1,
  arma::vec& p2, arma::vec& lower, arma::vec& upper, arma::uvec& islog);

arma::uvec GetIsR1(Rcpp::NumericVector modelAttr, std::string type);

void CheckPnames(Rcpp::List samples);
void CheckHyperPnames(Rcpp::List samples);

Rcpp::NumericMatrix na_matrix(unsigned int nr, unsigned int nc) ;

arma::uvec getArmaIdx(Rcpp::NumericVector x, unsigned int trueOrFalse) ;


arma::umat cellIdx2Mat(Rcpp::List data);

