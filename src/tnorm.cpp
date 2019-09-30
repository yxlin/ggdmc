#include <ggdmc.hpp>

tnorm::tnorm (double mu, double sig, double lower, double upper, bool lg) :
  m(mu), s(sig), l(lower), u(upper), lp(lg)
// Constructor. Initialize with mu and sigma. The default with no arguments
// is tnorm(0, 1, 0, Inf).
{
  // Bayesian optimization must allow bad values.
  // if (sig < 0.)
  // {
  //   Rcpp::Rcout << "Invalid sigma = " <<  sig << std::endl;
  //   Rcpp::stop("Bad sigma in dtnorm");
  // }
}

tnorm::tnorm (double mu, double sig, double lower, double upper, bool lg,
       bool lower_tail) :
  m(mu), s(sig), l(lower), u(upper), lp(lg), lt(lower_tail)
{
  if (sig < 0.)
  {
    Rcpp::Rcout << "Invalid sigma = " <<  sig << std::endl;
    Rcpp::stop("Bad sigma in ptnorm");
  }
}

tnorm::tnorm (double mu, double sig, double lower, double upper) :
  m(mu), s(sig), l(lower), u(upper)
{
  if (sig < 0.)
  {
    Rcpp::Rcout << "Invalid sigma = " <<  sig << std::endl;
    Rcpp::stop("Bad sigma in rtnorm");
  }
}


double tnorm::d (double x)
{
  double out, numer, denom;

  // if (s < 0)
  // {
  //   out = lp ? R_NegInf : 0;
  // }
  if  ((x >= l) && (x <= u))
  {
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    denom = R::pnorm(u, m, s, true, false) - R::pnorm(l, m, s, true, false);
    numer = R::dnorm(x, m, s, lp);
    out = lp ? (numer - std::log(denom)) : (numer/denom);
  }
  else
  {
    out = lp ? R_NegInf : 0;
  }

  return out;
}

void tnorm::d (std::vector<double> & x, std::vector<double> & output)
{
  for(size_t i=0; i<x.size(); i++)
  {
    output[i] = d(x[i]);
  }
}

double tnorm::d2 (double x)
{
  double out, numer, denom, sd;
  // s is precision

  if ((x >= l) && (x <= u)) {
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    sd = 1/std::sqrt(s);
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    denom = R::pnorm(u, m, sd, true, false) - R::pnorm(l, m, sd, true, false);
    numer = R::dnorm(x, m, sd, lp);
    out = lp ? (numer - std::log(denom)) : (numer/denom);
  } else {
    out = lp ? R_NegInf : 0;
  }
  return(out);

}

double tnorm::p (double x)
{
  double out, denom, qtmp;
  if (lt) {
    out = (x < l) ? 0 : 1;
  } else {
    out = (x < l) ? 1 : 0;
  }

  if ( (x >= l) && (x <= u) )
  {
    // 4th arg: lower.tail (lt)=1; 5th arg: log.p (lg)=0
    denom = R::pnorm(u, m, s, true, false) - R::pnorm(l, m, s, true, false);

    qtmp  = lt ?
    (R::pnorm(x, m, s, true, false) - R::pnorm(l, m, s, true, false)) :
      (R::pnorm(u, m, s, true, false) - R::pnorm(x, m, s, true, false)) ;

    out  = lp ? (std::log(qtmp)-std::log(denom)) : (qtmp/denom);
  }
  return out;
}

void tnorm::p (std::vector<double> & x, std::vector<double> & output)
{
  for(size_t i=0; i<x.size(); i++)
  {
    output[i] = p(x[i]);
  }
}
double tnorm::r ()
{
  // Standardised lower and upper
  double z, stdl, stdl2, stdu, stdu2, eq_a1, eq_a2;
  bool a0, a1, a2;
  stdl  = (l - m) / s; // l == stdlower, u == stdupper
  stdu  = (u - m) / s;
  stdl2 = stdl*stdl;
  stdu2 = stdu*stdu;

  // Accept-Reject Algorithm 0;
  // Algorithm (1): Use Proposition 2.3 with only lower truncation.
  // upper==INFINITY
  // rejection sampling with exponential proposal. Use if lower > mean
  // Algorithm (2): Use -x ~ N_+ (-mu, -mu^+, sigma^2) on page 123.
  // lower==-INFINITY
  // rejection sampling with exponential proposal. Use if upper < mean.
  // Algorithm (3, else): rejection sampling with uniform proposal.
  // Use if bounds are narrow and central.
  a0 = (stdl < 0 && u == R_PosInf) || (stdl == R_NegInf && stdu > 0) ||
    (std::isfinite(stdl) && std::isfinite(u) && stdl < 0 && stdu > 0 &&
    (stdu - stdl) > SQRT_2PI);
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
  return z*s + m;
}

//' Truncated Normal Distribution
//'
//' Random number generation, probability density and cumulative density
//' functions for truncated normal distribution.
//'
//' @param x,q vector of quantiles;
//' @param n number of observations. n must be a scalar.
//' @param p1 mean (must be scalar).
//' @param p2 standard deviation (must be scalar).
//' @param lower lower truncation value (must be scalar).
//' @param upper upper truncation value (must be scalar).
//' @param lt lower tail. If TRUE (default) probabilities are \code{P[X <= x]},
//' otherwise, \code{P[X > x]}.
//' @param lg log probability. If TRUE (default is FALSE) probabilities p are
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
//' x <- seq(-10, 10, length.out = 1e2)
//' mean <- 0
//' sd <- 1
//' lower <- 0
//' upper <- 5
//' dat1 <- ptnorm(x, 0, 1, 0, 5, lg = TRUE)
//' @export
// [[Rcpp::export]]
std::vector<double> dtnorm(std::vector<double> x, double p1, double p2,
                           double lower, double upper, bool lg = false) {
  if (upper < lower) Rcpp::stop("upper must be greater than lower.");
  if (p2 < 0) Rcpp::stop("sd must be greater than 0.\n");
  if (p2 == R_NegInf || p2 == R_PosInf) Rcpp::stop("sd must have a finite value.\n");
  if (p1 == R_NegInf || p1 == R_PosInf) Rcpp::stop("mean must have a finite value.\n");

  std::vector<double> out(x.size());

  tnorm * obj = new tnorm(p1, p2, lower, upper, lg);
  obj->d(x, out);

  delete obj;
  return out;
}


//' @rdname dtnorm
//' @export
// [[Rcpp::export]]
std::vector<double> rtnorm(unsigned int n, double p1, double p2, double lower,
                 double upper) {

  std::vector<double> out(n);
  tnorm * obj = new tnorm(p1, p2, lower, upper);
  for(size_t i = 0; i <n; i++) out[i] = obj->r();
  delete obj;
  return out;
}

//' @rdname dtnorm
//' @export
// [[Rcpp::export]]
std::vector<double> ptnorm(std::vector<double> q, double p1, double p2,
                           double lower, double upper, bool lt = true,
                           bool lg = false) {
  if (upper < lower) {Rcpp::stop("'upper' must be greater than 'lower'.");}
  if (p2 < 0)        {Rcpp::stop("'sd' must be greater than 0.\n");}
  if (p2 == R_NegInf   || p2 == R_PosInf)   {Rcpp::stop("'sd' must have a finite value.\n");}
  if (p1 == R_NegInf || p1 == R_PosInf) {Rcpp::stop("'mean' must have a finite value.\n");}

  std::vector<double>  out(q.size());
  tnorm * obj = new tnorm(p1, p2, lower, upper, lg, lt);
  obj->d(q, out);
  delete obj;
  return out;
}



