/* The original C version of the constant drift-diffusion density function is
* from fast-dm 30.2 density.c - compute the densities g- and g+ of the first
* exit time. Copyright (C) 2012  Andreas Voss, Jochen Voss.
* The vectorised and parallel integration functions are created by
* Yi-Shin Lin 2016 */
#include <ggdmc.hpp>
#include <vector>

// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP; e.g., OS X  clang
#ifdef _OPENMP
#include <omp.h>
#endif

#define EPSILON 1e-6
double TUNE_PDE_DT_MIN = 1e-6;
double TUNE_PDE_DT_MAX = 1e-6;
double TUNE_PDE_DT_SCALE = 0.0;

double  TUNE_DZ;
double  TUNE_DV;
double  TUNE_DT0;

double  TUNE_INT_T0;
double  TUNE_INT_Z;

int  precision_set = 0;

void set_precision (double p)
{
  /* Try to achieve an accuracy of approximately 10^{-p} for the CDF. */
  TUNE_PDE_DT_MIN = std::pow(10.0, -0.400825*p-1.422813);
  TUNE_PDE_DT_MAX = std::pow(10.0, -0.627224*p+0.492689);
  TUNE_PDE_DT_SCALE = std::pow(10.0, -1.012677*p+2.261668);
  TUNE_DZ = std::pow(10.0, -0.5*p-0.033403);
  TUNE_DV = std::pow(10.0, -1.0*p+1.4);
  TUNE_DT0 = std::pow(10.0, -0.5*p-0.323859);

  TUNE_INT_T0 = 0.089045 * std::exp(-1.037580*p);
  TUNE_INT_Z = 0.508061 * std::exp(-1.022373*p);

  precision_set = 1;
}

double integrate(double (*F)(double, std::vector<double>&),
                 std::vector<double>& pVec, double a, double b,
                 double step_width)
{
  double width = b-a ;  // integration width
  int N = std::max(4, (int)(width / step_width)) ; // N at lease == 4
  double step = width / N;
  double x ;
  double out = 0 ;
  for(x = a + 0.5*step; x < b; x += step) { out += step * F(x, pVec) ; }
  return out ;
}

double g_minus_small_time(double DT, double zr, int N)
{
  int i;
  double sum = 0 ;
  for(i = -N/2; i <= N/2; i++)
  {
    double d = 2*i + zr ; \
    sum += exp(-d*d / (2*DT)) * d ;
  }
  return sum / sqrt(2*M_PI*DT*DT*DT) ;
}

double g_minus_large_time(double DT, double zr, int N)
{
  int i ;
  double sum = 0 ;
  for(i = 1; i <= N; i++)
  {
    double d = i * M_PI ;
    sum += std::exp(-0.5 * d*d * DT) * sin(d*zr) * i ;
  }
  return sum * M_PI ;
}

double g_minus_no_var(double DT, double a, double zr, double v)
{
  int N_small, N_large;
  double simple, factor, eps, out;
  double ta = DT/(a*a);

  factor = std::exp(-a*zr*v - 0.5*v*v*DT) / (a*a); // Front term in A3
  eps = EPSILON / factor;

  N_large = (int)ceil(1/(M_PI*sqrt(DT)));
  if (M_PI*ta*eps < 1)
  {
    N_large = std::max(N_large,
                       (int)ceil(std::sqrt(-2.0*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
  }

  if (2*sqrt(2*M_PI*ta)*eps < 1)
  {
    N_small = (int)ceil(std::max(sqrt(ta) + 1,
                        2 + std::sqrt(-2.0*ta*log(2.0*eps*std::sqrt(2.0*M_PI*ta)))));
  } else
  {
    N_small = 2;
  }

  if (N_small < N_large)
  {
    simple = g_minus_small_time(DT/(a*a), zr, N_small);
  } else
  {
    simple = g_minus_large_time(DT/(a*a), zr, N_large);
  }

  out = std::isinf(factor) ? 0 : (factor * simple) ;
  return out;
}

double integral_v_g_minus(double zr, std::vector<double>& pVec)
{
  /* 0  1  2  3  4   5   6  7    8   9
  * a  v  zr d  sz  sv  t0 st0  DT  precision
  * position 2 has been converted to zr in g_plus and g_minus
  * position 8 has been converted to DT in _z_g_minus
  */
  double DT = pVec[8] ; // Let's make it clear. It's DT
  double a  = pVec[0] ;
  double v  = pVec[1] ;
  double sv = pVec[5] ; // this should be sv
  int N_small, N_large ;
  double simple, factor, eps, out ;
  double ta = DT/(a*a) ;

  factor = 1 / (a*a * std::sqrt(DT * sv*sv + 1)) *
    std::exp(-0.5 * (v*v*DT + 2.0*v*a*zr - a*zr*a*zr*sv*sv) / (DT*sv*sv + 1.0)) ;
  eps = EPSILON / factor;

  N_large = (int)ceil(1.0 / (M_PI*std::sqrt(DT))) ;
  if (M_PI*ta*eps < 1)
  {
    N_large = std::max(N_large,
                       (int)std::ceil(std::sqrt(-2.0*std::log(M_PI*ta*eps) / (M_PI*M_PI*ta)))) ;
  }

  if (2*std::sqrt(2.0*M_PI*ta)*eps < 1)
  {
    N_small = (int)ceil(std::max(std::sqrt(ta)+1,
                        2+std::sqrt(-2.0*ta*log(2.0*eps*std::sqrt(2.0*M_PI*ta))))) ;
  } else
  {
    N_small = 2 ;
  }

  if (std::isinf(factor))
  {
    out = 0 ;
  } else if (sv == 0)
  {
    out = g_minus_no_var(DT, a, zr, v) ;
  } else if (N_small < N_large)
  {
    simple = g_minus_small_time(DT/(a*a), zr, N_small) ;
    out = factor * simple ;
  } else
  {
    simple = g_minus_large_time(DT/(a*a), zr, N_large) ;
    out = factor * simple ;
  }
  return out ;
}

double integral_z_g_minus(double DT, std::vector<double>& pVec)
{
  /* 0  1  2  3  4   5   6  7    8   9
  * a  v  zr d  sz  sv  t0 st0  RT  precision
  * position 2 has been converted to zr in g_plus and g_minus
  */
  double out ;
  pVec[8] = DT ; // Convert RT to DT, so integral_v_g_minus's pVec[8] is DT

  if (DT <= 0) // if DT <= 0
  {
    out = 0 ;
  } else if (pVec[4] == 0) // this should be sz
  {
    out = integral_v_g_minus(pVec[2], pVec);
  } else
  {
    double a = pVec[2] - 0.5*pVec[4] ;   // zr - 0.5*szr; uniform variability
    double b = pVec[2] + 0.5*pVec[4] ;   // zr + 0.5*szr
    double step_width = TUNE_INT_Z ;     // TUNE_INT_Z
    /* pVec goes into _v_g_minus becomes
    * 0  1  2  3   4   5   6  7    8   9
    * a  v  zr d  sz  sv  t0 st0  DT  precision
    */
    out = integrate(integral_v_g_minus, pVec, a, b, step_width) / pVec[4] ;
  }
  return out ;
}

double integral_t0_g_minus(double DT, std::vector<double>& pVec)
{
  /* 0  1  2  3   4   5   6  7    8   9
  * a  v  zr d  sz  sv  t0 st0  RT  precision
  * position 2 has been converted to zr in g_plus and g_minus
  */
  double out;
  if (pVec[7] == 0)
  {
    out = integral_z_g_minus(DT, pVec); // should send t as RT
  } else
  {
    double a = DT - 0.5*pVec[7] ;   // DT - 0.5*st0
    double b = DT + 0.5*pVec[7] ;   // DT + 0.5*st0
    double step_width = TUNE_INT_T0 ;
    out = integrate(integral_z_g_minus, pVec, a, b, step_width) / pVec[7] ;
  }
  return out;
}

//' Calculate Drift-diffusion Probability Density
//'
//' \code{g_minus} and \code{g_plus} implement A1 to A4 equations in Voss,
//' Rothermund, and Voss (2004). These equations calculate Ratcliff's
//' drift-diffusion model (1978). This source codes are derived from
//' Voss & Voss's fast-dm 30.2 in density.c.
//'
//' Two parallel functions \code{g_minus_parallel} and \code{g_plus_parallel},
//' using OpenMP libraries to do numerical integration. They resolve the
//' problem when high precision (> 10) is required.
//'
//' @param pVec a 9-element parameter (double) vector. The user has to follow
//' the sequence strictly. a, v, zr, d, sz, sv, t0, st0, RT, precision.
//' @references Voss, A., Rothermund, K., & Voss, J. (2004). Interpreting the
//' parameters of the diffusion model: A empirical validation
//' \emph{Memory and Cognition}, \bold{32(7)}, 1206--1220. \cr\cr
//' Ratcliff, R (1978). A theory of memory retrieval. \emph{Psychology Review},
//' \bold{85(2)}, 59--108.
//' @export
//' @examples
//' pvec1 <- c(a=2, v=2.5, zr=0.5, d=0, sz=0.3, sv=1, t0=0.3, st0=0,
//'            RT=.550, precision=2.5)
//' g_minus(pvec1) ## 0.04965882
//' g_plus(pvec1) ## 2.146265
//'
//' pvec2 <- c(a=2, v=2.5, zr=0.5, d=.2, sz=0.3, sv=1, t0=0.3, st0=.1,
//'            RT=.550, precision=2.5)
//'
//' g_minus(pvec2) ## 0.04194364
//' g_plus(pvec2)  ## 1.94957
//'
// [[Rcpp::export]]
double g_minus(std::vector<double> pVec)
{
  /* Let just use one vector, passing also RT and precision
  * 0  1   2   3   4   5   6  7    8   9
  * a  v  zr  d  sz  sv  t0 st0  RT  precision (dmc uses names)
  */
  set_precision(pVec[9]) ;
  double DT = pVec[8] - pVec[6] - 0.5*pVec[3] ;
  return integral_t0_g_minus(DT, pVec);
}

//' @rdname g_minus
//' @export
// [[Rcpp::export]]
double g_plus(std::vector<double> pVec)
{
  /* Let just use one vector, passing also RT and precision
  ----------Name sequence correction----------
  * 0  1  2   3  4   5   6  7    8   9
  * a  v  zr  d  sz  sv  t0 st0  RT  precision (dmc uses names)
  * Note rtdists take z as zr.
  * Should not pass pVec's mem location, because here pVec 1(v) and 2(zr)
  * is flipped and assign back to its internal mem location
  */
  set_precision(pVec[9]) ;
  double DT = pVec[8] - pVec[6] - 0.5*pVec[3] ;
  // pVec[2] = pVec[2] / pVec[0] ; // make z relative to a (becoming zr)
  pVec[2] = 1 - pVec[2] ;  // positive bound reverses zr to 1-zr
  pVec[1] = -pVec[1] ;     // positive bound reverses v to -v
  return integral_t0_g_minus(DT, pVec);
}


