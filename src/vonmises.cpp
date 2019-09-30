//    Copyright (C) <2019>  <Yi-Shin Lin>
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#include <ggdmc.hpp>

using namespace Rcpp;

inline double pvm_(double q, double kappa, double tol) {
  bool flag = true;
  double p = 1;
  double sum = 0;
  while (flag) {
    double term = (R::bessel_i(kappa, p, 1.) * std::sin(p * q)) / p;
    sum += term;
    p++;
    if (std::abs(term) < tol) flag = false;
  }

  double out = q/M_2PI + sum / (M_PI * R::bessel_i(kappa, 0, 1.));
  return out;
}

inline arma::vec fmod(arma::vec dividend, double divisor)
{
  // a float point modulus operator taking armadillo vector and a double-type
  // divisor
  return dividend - arma::floor(dividend / divisor)*divisor;
}

//' Generate random deviates from a von Mises distribution
//'
//' This function generates random numbers in radian unit from a von Mises
//' distribution using the location (ie mean) parameter, mu and the
//' concentration (ie precision) parameter kappa.
//'
//' A random number for a circular normal distribution has the form:\cr
//' \deqn{f(theta; mu, kappa) = 1 / (2*pi*I0(kappa)) * exp(kappa*cos(theta-mu))}
//' theta is between 0 and 2*pi.
//'
//' \code{I0(kappa)} in the normalizing constant is the modified Bessel
//' function of the first kind and order zero.
//'
//' @param n number of observations. Must be a scalar.
//' @param mu mean direction of the distribution. Must be a scalar.
//' @param kappa concentration parameter. A positive value
//' for the concentration parameter of the distribution. Must be a scalar.
//'
//' @return a column vector
//' @references
//' \enumerate{
//' Ulric Lund, Claudio Agostinelli, et al's  (2017). R package 'circular':
//' Circular Statistics (version 0.4-91).
//' \url{https://r-forge.r-project.org/projects/circular/}
//' }
//' @examples
//' n  <- 1e2
//' mu <- 0
//' k  <- 10
//'
//' \dontrun{
//' vm1 <- circular:::RvonmisesRad(n, mu, k)
//' vm2 <- rvm(n, mu, k)
//' vm3 <- circular:::conversion.circular(circular:::circular(vm1))
//' vm4 <- circular:::conversion.circular(circular:::circular(vm2))
//' plot(vm3)
//' plot(vm4)
//' }
//' @export
// [[Rcpp::export]]
arma::vec rvonmises(unsigned int n, double mu, double kappa) {

  double U, a, b, r, z, f, c, tmp;
  arma::vec out(n);
  arma::vec::iterator i = out.begin() ;

  // If kappa is small, sample angles from a uniform distribution [0 2*pi]
  if (kappa < 1e-10 ) {
    do { *i =  R::runif(0.0, M_2PI); i++; } while (i < out.end());
  } else {
    a = 1.0 + std::sqrt(1.0 + 4.0 * kappa * kappa);
    b = (a - std::sqrt(2.0 * a)) / (2.0 * kappa);
    r = (1.0 + b * b) / (2.0 * b);

    do {
      z  = std::cos(M_PI * R::runif(0.0, 1.0));
      f  = (1.0 + r*z) / (r + z);
      c  = kappa*(r - f);
      U = R::runif(0.0, 1.0);

      if (c * (2.0 - c) > U) {
        tmp = (R::runif(0.0, 1.0) > 0.50) ? std::acos(f) + mu : -std::acos(f) + mu;
        *i = tmp - std::floor(tmp/M_2PI)*M_2PI; // store in out
        i++;
      } else {
        if (std::log(c/U) + 1.0 >= c) {
          tmp = (R::runif(0.0, 1.0) > 0.50) ? std::acos(f) + mu : -std::acos(f) + mu;
          *i = tmp - std::floor(tmp/M_2PI)*M_2PI;
          i++;
        }
      }
    } while(i < out.end());
  }

  return out;
}


//' @rdname rvonmises
//' @export
// [[Rcpp::export]]
arma::vec dvonmises(arma::vec x, double mu, double kappa)
{
  arma::vec out(x.n_elem);

  if (kappa == 0) {
    out.fill(1/M_2PI);
  } else if (kappa < 1e5) {
    out = arma::exp(kappa * arma::cos(x - mu)) /
      ( M_2PI * R::bessel_i(kappa, 0,  1.) );
  } else {
    for (size_t i = 0; i < x.n_elem; i++) {
      double num = x(i) - mu;
      double tmp = std::fmod(num, M_2PI);
      out(i) = tmp == 0 ?  INFINITY : 0;
    }
  }

  return out;
}

//' @rdname rvonmises
//' @export
// [[Rcpp::export]]
arma::vec pvonmises(arma::vec q, double mu, double kappa, double tol = 1e-20)
{
  arma::vec qmod = fmod(q, M_2PI);
  unsigned int n = q.n_elem;
  // double mu_mod = std::fmod(mu, M_2PI);
  arma::vec out(n); out.fill(NA_REAL);
  if (mu == 0) {
    for (size_t i = 0; i < qmod.n_elem; i++)
      out(i) = pvm_(qmod(i), kappa, tol);
  } else {
    double upper, lower;
    for (size_t i = 0; i < qmod.n_elem; i++) {
      if (qmod(i) <= mu) {
        upper = std::fmod(qmod(i) - mu, M_2PI);
        if (upper == 0) upper = M_2PI;
        lower = std::fmod(-mu, M_2PI);
        out(i) = pvm_(upper, kappa, tol) - pvm_(lower, kappa, tol);
      } else {
        upper = qmod(i) - mu;
        lower = std::fmod(mu, M_2PI);
        out(i) = pvm_(upper, kappa, tol) + pvm_(lower, kappa, tol);
      }
    }
  }

  return out;
}

