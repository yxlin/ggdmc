#ifndef VONMISES_HPP
#define VONMISES_HPP
#include <RcppArmadillo.h>

arma::vec rvonmises(unsigned int n, double mu, double kappa);

#endif