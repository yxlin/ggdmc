#include <RcppArmadillo.h>

void set_seed(unsigned int seed);

double rtn_scalar(double mean,  double sd, double l, double u);

double dtn_scalar(double x, double mean, double sd, double lower,
  double upper, bool lp);

double ptn_scalar(double q, double mean, double sd, double lower, double upper,
                  bool lt, bool lp);

arma::vec dtnorm(arma::vec x, double mean, double sd, double lower,
  double upper, bool lp);

arma::vec rtnorm(unsigned int n, double mean, double sd, double lower,
                 double upper);

arma::vec ptnorm(arma::vec q, double mean, double sd, double lower,
  double upper, bool lt, bool lp);
