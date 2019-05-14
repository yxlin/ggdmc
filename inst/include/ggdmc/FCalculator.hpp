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

// A specific two-step constructor,  two member functions and a destructor
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
