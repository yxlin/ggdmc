#ifndef TNORM_HPP
#define TNORM_HPP

#include <RcppArmadillo.h> // for std, cmath and many other supports via Rcpp

#define SQRT_2PI   2.5066282746310007e+0 /* sqrt(2 x pi) */

class tnorm {
public:
  // truncated normal distribution.
  double m, s, l, u;  // mean, sd, lower, upper, and precision
  bool lp, lt;        // log probability, lower tail

  tnorm (double mu, double sig, double lower, double upper, bool lg);

  tnorm (double mu, double sig, double lower, double upper, bool lg,
         bool lower_tail);
  tnorm (double mu, double sig, double lower, double upper);

  double d (double x);
  // Return probability density function.
  double p (double x);
  // Return cumulative distribution function.
  double r ();
  // Return random deviates

  void d (std::vector<double> & x, std::vector<double> & output);
  void p (std::vector<double> & x, std::vector<double> & output);

  double d2 (double x);
  // Return probability density function, using precision.



private:
  double rtnorm0(const double &l, const double &u) {
    // Accept-Reject Algorithm 0; Naive method A-R method
    bool invalid = true;
    double z;
    while (invalid) {
      z = R::rnorm(0.0, 1.0);
      if (z <= u && z >= l) break;
    }
    return z;
  }

  double rtnorm1(const double &l, const double &u) {
    // Algorithm 1; 'expl'; use when lower > mean; upper = INFINITY; p 122, right
    bool invalid = true;
    double z, r, num; // a stands for alphaStar in Robert (1995)
    double a = 0.5 * (std::sqrt(l*l + 4.0) + l);

    while (invalid) {
      z   = (-1.0/a)*std::log(R::runif(0.0, 1.0)) + l; // control lower boundary
      num = R::runif(0.0, 1.0);
      r   = std::exp(-0.5 * (z - a)*(z - a));
      if (num <= r && z <= u) break;
    }
    return z ;
  }

  double rtnorm2(const double &l, const double &u) {
    // Algorithm 2; 'expu'; use when upper < mean; lower = -INFINITY.
    bool invalid = true;
    double z, r, num;
    double a = 0.5 * (std::sqrt(u*u + 4.0) - u);

    while (invalid) {
      z   = (-1.0/a)*std::log(R::runif(0.0, 1.0)) - u; // control lower boundary
      num = R::runif(0.0, 1.0);
      r   = std::exp(-0.5 * (z - a)*(z - a));
      if (num <= r && z <= -l) break;
    }
    return -z;  // note the negative
  }

  double rtnorm3(const double &l, const double &u) {
    // Algorithm 3; page 123. 2.2. Two-sided truncated normal dist.
    bool invalid = true;
    double z, r, num; // a stands for alphaStar in Robert (1995)
    double l2 = l*l;
    double u2 = u*u;

    while (invalid) {
      z = R::runif(l, u) ;
      if (l > 0) {
        r = std::exp(0.5 * (l2 - z*z));
      } else if (u < 0) {
        r = std::exp(0.5 * (u2 - z*z));
      } else  {
        r = std::exp( -0.5 * z * z ) ;
      }
      num = R::runif(0.0, 1.0);
      if (num <= r) break;
    }
    return z ;
  }

};

#endif
