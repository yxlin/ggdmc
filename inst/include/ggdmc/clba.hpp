#include <RcppArmadillo.h>

arma::vec remove_t0(arma::vec rt, double t0);

arma::vec fptcdf(arma::vec rt, double A, double b, double mean_v, double sd_v,
  double t0, bool posdrift);

arma::vec fptpdf(arma::vec rt, double A, double b, double mean_v, double sd_v,
  double t0, bool posdrift);

arma::vec n1PDFfixedt0(arma::vec rt, arma::vec A, arma::vec b, arma::vec mean_v,
                       arma::vec sd_v, arma::vec t0, bool posdrift);

arma::vec n1PDFfixedt0_pda(arma::vec rt, double A, double b,
                           arma::mat mean_v, arma::vec sd_v, double t0,
                           unsigned int n, double h, bool debug);

arma::mat rlba_norm(unsigned int n, arma::vec A, arma::vec b, arma::mat mean_v,
                    arma::vec sd_v, arma::vec t0, arma::vec st0, bool posdrift,
                    bool return_ttf, bool debug);

arma::vec n1PDF_cnorm(arma::vec rt, arma::vec A, arma::vec b, arma::vec t0,
  arma::vec mean_v, arma::vec sd_v, arma::vec st0, double corr_v,
  unsigned int n, double h, bool debug);

Rcpp::NumericVector n1PDF_gpu(arma::vec x, double A, double b, arma::vec mean_v,
  arma::vec sd_v, double t0, int n, unsigned int nthread, unsigned int gpuid,
  double bw, bool debug);

Rcpp::NumericVector n1PDF_plba0_gpu(arma::vec x, double A, double b,
  arma::vec mean_v, arma::vec sd_v, double t0, arma::vec mean_w, double rD,
  double swt, unsigned int n, unsigned int nthread, unsigned int gpuid,
  double bw, bool debug);

Rcpp::NumericVector n1PDF_plba1_gpu(arma::vec x, double A, double b,
  arma::vec mean_v, arma::vec sd_v, double t0, arma::vec mean_w, double rD,
  double swt, unsigned int n, unsigned int nthread, unsigned int gpuid,
  double bw, bool debug);
