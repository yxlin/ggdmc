#include <ggdmc.hpp>

using namespace Rcpp;

#define SQRT_2PI   2.5066282746310007e+0 /* sqrt(2 x pi) */

void set_seed(unsigned int seed) {
  Environment base_env("package:base");
  Function set_seed_r = base_env["set.seed"];
  set_seed_r(seed);
}

inline double rtnorm0(const double l, const double u) {
  // Accept-Reject Algorithm 0; Naive method A-R method
  bool invalid = true;
  double z;
  while (invalid) {
    z = R::rnorm(0.0, 1.0);
    if (z <= u && z >= l) break;
  }
  return z;
}

inline double rtnorm1(const double l, const double u) {
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

inline double rtnorm2(const double l, const double u) {
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

inline double rtnorm3(const double l, const double u) {
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

// [[Rcpp::export]]
double rtn_scalar(double mean,  double sd, double l, double u) {
  double z, stdl, stdl2, stdu, stdu2, eq_a1, eq_a2; // Standardised lower and upper
  bool a0, a1, a2;
  stdl  = (l - mean) / sd; // l == stdlower, u == stdupper
  stdu  = (u - mean) / sd;
  stdl2 = stdl*stdl;
  stdu2 = stdu*stdu;

  // Accept-Reject Algorithm 0;
  // Algorithm (1): Use Proposition 2.3 with only lower truncation. upper==INFINITY
  // rejection sampling with exponential proposal. Use if lower > mean
  // Algorithm (2): Use -x ~ N_+ (-mu, -mu^+, sigma^2) on page 123. lower==-INFINITY
  // rejection sampling with exponential proposal. Use if upper < mean.
  // Algorithm (3, else): rejection sampling with uniform proposal.
  // Use if bounds are narrow and central.
  a0 = (stdl < 0 && u == INFINITY) || (stdl == -INFINITY && stdu > 0) ||
    (std::isfinite(stdl) && std::isfinite(u) && stdl < 0 && stdu > 0 && (stdu - stdl) > SQRT_2PI);
  eq_a1 = stdl + (2.0 * std::sqrt(M_E) / (stdl + std::sqrt(stdl2 + 4.0))) *
    (std::exp( 0.25 * (2.0*stdl - stdl*std::sqrt(stdl2 + 4.0))));
  a1 = (stdl >= 0) && (stdu > eq_a1);
  eq_a2 = -stdu + (2.0 * std::sqrt(M_E) / (-stdu + std::sqrt(stdu2 + 4.0))) *
    (std::exp(0.25 * (2.0*stdu + stdu*std::sqrt(stdu2 + 4.0))));
  a2 = (stdu <= 0) && (-stdl > eq_a2);

  if (a0) {
    z = rtnorm0(stdl, stdu);
  } else if (a1) {
    z = rtnorm1(stdl, stdu);
  } else if (a2) {
    z = rtnorm2(stdl, stdu);
  } else {
    z = rtnorm3(stdl, stdu);
  }
  return(z*sd + mean);
}

double dtn_scalar(double x, double mean, double sd, double lower,
  double upper, bool lp) {

  double out, numer, denom;
  if ((x >= lower) && (x <= upper)) {
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    denom = R::pnorm(upper, mean, sd, true, false) -
      R::pnorm(lower, mean, sd, true, false);
    numer = R::dnorm(x, mean, sd, lp);
    out = lp ? (numer - std::log(denom)) : (numer/denom);
  } else {
    out = lp ? -INFINITY : 0;
  }
  return(out);
}

double ptn_scalar(double q, double mean, double sd, double lower, double upper,
  bool lt, bool lp) {

  double out, denom, qtmp;
  if (lt) {
    out = (q < lower) ? 0 : 1;
  } else {
    out = (q < lower) ? 1 : 0;
  }
  if ((q >= lower) && (q <= upper)) {
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    denom = R::pnorm(upper, mean, sd, true, false) -
      R::pnorm(lower, mean, sd, true, false);
    qtmp  = lt ? (R::pnorm(q, mean, sd, true, false) - R::pnorm(lower, mean, sd, true, false)) :
      (R::pnorm(upper, mean, sd, true, false) - R::pnorm(q, mean, sd, true, false));
    out  = lp ? (std::log(qtmp)-std::log(denom)) : (qtmp/denom);
  }
  return(out);
}


//' Truncated Normal Distribution
//'
//' Random number generation, probability density and cumulative density
//' functions for truncated normal distribution.
//'
//' @param x,q vector of quantiles;
//' @param n number of observations. n must be a scalar.
//' @param mean mean (must be scalar).
//' @param sd standard deviation (must be scalar).
//' @param lower lower truncation value (must be scalar).
//' @param upper upper truncation value (must be scalar).
//' @param lt lower tail. If TRUE (default) probabilities are \code{P[X <= x]},
//' otherwise, \code{P[X > x]}.
//' @param log log probability. If TRUE (default is FALSE) probabilities p are
//' given as \code{log(p)}.
//' @return a column vector.
//' @examples
//' ## rtn example
//' dat1 <- rtnorm(1e5, 0, 1, 0, Inf)
//' hist(dat1, breaks = "fd", freq = FALSE, xlab = "",
//'      main = "Truncated normal distributions")
//'
//' ## dtn example
//' x <- seq(-5, 5, length.out = 1e3)
//' dat1 <- dtnorm(x, 0, 1, -2, 2, 0)
//' plot(x, dat1, type = "l", lwd = 2, xlab = "", ylab= "Density",
//'      main = "Truncated normal distributions")
//'
//' ## ptn example
//' x <- seq(-50, 10, length.out = 1e3)
//' mean <- 0
//' sd <- 1
//' lower <- 0
//' upper <- 5
//' dat1 <- ptnorm(x, 0, 1, 0, 5, log = TRUE)
//' @export
// [[Rcpp::export]]
arma::vec dtnorm(arma::vec x, double mean, double sd, double lower,
  double upper, bool log = false) {
  if (upper < lower) stop("upper must be greater than lower.");
  if (sd < 0) stop("sd must be greater than 0.\n");
  if (sd==R_NegInf || sd==R_PosInf) stop("sd must have a finite value.\n");
  if (mean==R_NegInf || mean==R_PosInf) stop("mean must have a finite value.\n");

  arma::vec out(x.n_elem);
  for (size_t i = 0; i < x.n_elem; i++)
    out(i) = dtn_scalar(x(i), mean, sd, lower, upper, log);
  return out;
}

//' @rdname dtnorm
//' @export
// [[Rcpp::export]]
arma::vec rtnorm(unsigned int n, double mean, double sd, double lower,
                 double upper) {
  arma::vec out(n);
  for (size_t i = 0; i < n; i++)
    out(i) = rtn_scalar(mean, sd, lower, upper);
  return out;
}

//' @rdname dtnorm
//' @export
// [[Rcpp::export]]
arma::vec ptnorm(arma::vec q, double mean, double sd, double lower,
  double upper, bool lt = true, bool log = false) {
  if (upper < lower) {Rcpp::stop("'upper' must be greater than 'lower'.");}
  if (sd < 0)        {Rcpp::stop("'sd' must be greater than 0.\n");}
  if (sd==R_NegInf   || sd==R_PosInf)   {Rcpp::stop("'sd' must have a finite value.\n");}
  if (mean==R_NegInf || mean==R_PosInf) {Rcpp::stop("'mean' must have a finite value.\n");}

  arma::vec out(q.n_elem);
  for(size_t i = 0; i < q.n_elem; i++)
    out[i] = ptn_scalar(q[i], mean, sd, lower, upper, lt, log);
  return out;
}
