/* (1) pde.c - numerically solve the Fokker-Planck equation
 * (2) phi.c - the CDF and inverse CDF of the standard normal distribution
 * (3) cdf.c, FCalculator.h, FControler.h & CDF_***.h - compute the CDF for
 * the diffusion model
 * (4) construct-samples.c, Sampling.hpp - Contains main call for random
 * sampling
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301 USA.
 */

// See also Feller (1971, p358 & p359)

#ifndef FCALCULATOR_H
#define FCALCULATOR_H

#include <RcppArmadillo.h> // for std, cmath and many other supports via Rcpp
#include <gsl/gsl_integration.h>


#define xrenew(T,OLD,N) ( (T *)xrealloc(OLD, (N)*sizeof(T)) )

class F_calculator
// A parent class for F_plain, F_sz, F_sv & F_st0. All four derived classes
// should all return F_calculator type.
{
public:
  int N, plus;
  void *data; // Cast to a derived class in their member functions

  void (*start)   (F_calculator *, int plus);
  void (*free)    (F_calculator *);
  const double *(*get_F)  (F_calculator *, double t);
  double (*get_z) (const F_calculator *, int i);
};

class F_plain_data //  plain: no variability
{
public:
  double  a, v, t0, d;	/* parameters (except z) */
  double  dz;				    /* z step-size */
  double  t_offset;		  /* time adjustment, resulting from t0 and d */
  double  t;				    /* adjusted time, corresponding to the vector F */
  double *F;  // state at time t + t_offset; ie CDFs

  double TUNE_PDE_DT_MIN;
  double TUNE_PDE_DT_MAX;
  double TUNE_PDE_DT_SCALE;

  double F_limit(double z);
};

class F_sz_data // sz
{
public:
  F_calculator *base_fc;    // gives the values we average over
  double *avg;              // the computed averages
  int  k;                   // the average involves 2*k+1 cells
  double  q;                // unused part of the outermost cells
  double  f;                // scale factor for the integration
};

class F_sv_data // sv
{
public:
  int  nv;                // number of points in integration
  std::vector<F_calculator*> base_fc; // F_calculators for different v
  double *avg;
};

class F_st0_data // st0
{
public:
  F_calculator *base_fc;
  double  st0;		// variability of t0
  int     M;			// number of stored grid lines
  double  start;	// t-value of first stored grid line
  double  dt;		  // t-spacing of stored grid lines
  double *values;	// array: stored grid lines (length M*(N+1))
  char   *valid;	// which lines in 'values' are valid
  int     base;		// first grid line starts at pos. base*(N+1)
  double *avg;		// the computed average (size N+1)
};

class CDF // For FCalculator_new
{
private:
  double logMill(double x) // log of Mill's ratio
  {
    double out = x >= 1e4 ? -std::log(x) :
      R::pnorm(x, 0, 1, 0, 1) - R::dnorm(x, 0, 1, 1);
    return out;
  }

public:
  double a;     // 0 Boundary separation
  double v;     // 1 Mean of the drift
  double mean_zr;    // 2 Mean of diffusion starting point relative to boundaries
  double d;     // 3 Difference between boundaries of non-decision time
  double szr;   // 4 width of zr distribution
  double sv;    // 5 standard deviation of v distribution
  double mean_t0;    // 6 mean non-decision time
  double st0;   // 7 width of t0 distribution
  double s;     // 8 standard deviation; sqrt(diffusion constant)
  double sv2, v2, a2, zr2, eps;

  double DT, lower, upper, zr;    // align with gsl integration routines

  // Specialized constructor
  CDF(std::vector<double> params)
  {
    a   = params[0];
    v   = params[1];
    sv  = params[5];
    szr = params[4];

    mean_zr = params[2];
    if (szr == 0) zr = mean_zr;

    sv2 = sv*sv;
    v2  = v*v;
    a2  = a*a;
    eps = std::sqrt(DBL_EPSILON);
  };
  CDF(Parameters * params)
  {
    a   = params->a;
    v   = params->v;
    sv  = params->sv;
    szr = params->szr;

    mean_zr = params->zr;
    if (szr == 0) zr = mean_zr; // if szr != 0, zr will be unset.

    sv2 = sv*sv;
    v2  = v*v;
    a2  = a*a;
    eps = std::sqrt(DBL_EPSILON);

  };

  double G0_over_zr (double t, double zr)
  {
    // Distribution function at lower barrier (Eq. 3 of the article)
    double logphi, rj, x0, x1, gj;
    double sqt  = std::sqrt(t);
    double sqet = sqt * std::sqrt(1. + sv2*t);

    double G = 0;
    unsigned int j = 0;
    bool flag = true;

    while (flag)
    {
      // Odd
      rj = (double)j*a + a*zr;
      logphi = R::dnorm(rj/sqt, 0, 1, 1);
      x0 = (rj - v*t + sv2*(rj + a*zr) * t) / sqet;
      x1 = (rj + v*t + sv2*(rj - a*zr) * t) / sqet;
      gj = std::exp(logphi + logMill(x0)) + std::exp(logphi + logMill(x1));
      G += gj;

      flag = gj >= eps;
      if (!flag) break;
      j++;

      // Even
      rj = (double)j*a + a*(1.-zr);
      logphi = R::dnorm(rj/sqt, 0, 1, 1);
      x0 = (rj - v*t + sv2*(rj + a*zr) * t) / sqet;
      x1 = (rj + v*t + sv2*(rj - a*zr) * t) / sqet;
      gj = std::exp(logphi + logMill(x0)) + std::exp(logphi + logMill(x1));
      G -= gj;

      j++;
    }
    return std::exp(.5*(-v2*t - 2.*v*a*zr + sv2*a2*zr*zr) / (1. + sv2*t) ) * G;
  }

  double G0 (double t)
  {
    // Distribution function at lower barrier (Eq. 3 of the article)
    // when mean_zr == zr
    double logphi, rj, x0, x1, gj;
    double sqt  = std::sqrt(t);
    double sqet = sqt * std::sqrt(1. + sv2*t);

    double G = 0;
    unsigned int j = 0;
    bool flag = true;

    while (flag)
    {
      // Odd
      rj = (double)j*a + a*zr;
      logphi = R::dnorm(rj/sqt, 0, 1, 1);
      x0 = (rj - v*t + sv2*(rj + a*zr) * t) / sqet;
      x1 = (rj + v*t + sv2*(rj - a*zr) * t) / sqet;
      gj = std::exp(logphi + logMill(x0)) + std::exp(logphi + logMill(x1));
      G += gj;

      flag = gj >= eps;
      if (!flag) break;
      j++;

      // Even
      rj = (double)j*a + a*(1.-zr);
      logphi = R::dnorm(rj/sqt, 0, 1, 1);
      x0 = (rj - v*t + sv2*(rj + a*zr) * t) / sqet;
      x1 = (rj + v*t + sv2*(rj - a*zr) * t) / sqet;
      gj = std::exp(logphi + logMill(x0)) + std::exp(logphi + logMill(x1));
      G -= gj;

      j++;
    }
    return std::exp(.5*(-v2*t - 2.*v*a*zr + sv2*a2*zr*zr) / (1. + sv2*t) ) * G;
  }


};

// A specific two-step constructor, two member functions and a destructor
// of the F_calculator class
F_calculator * F_new (Parameters *params);
void           F_start (F_calculator *fc, int boundary);
const double * F_get_F (F_calculator *fc, double t);
double         F_get_z (const F_calculator *fc, int i);
void           F_delete (F_calculator *fc);
double         F_get_val (F_calculator *fc, double t, double z);

/*---------------------------------------------------------------------------
 Original found in construct-samples.c
 ---------------------------------------------------------------------------*/
Rcpp::List sampling(int s_size, Parameters * params, bool random_flag);

#endif
