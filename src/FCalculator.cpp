/* - Rcpp::sampling and Rcpp::r_fastdm are adapted from src folder in
 * rtdists 0.9-0 by Singmann, Brown, Gretton, Heathcote, Voss, Voss & Terry.
 *
 * - cdf.c - compute the CDF for the diffusion model (F* functions)
 * - pde.c - numerically solve the Fokker-Planck equation (solve_tridiag,
 * make_step & advance_to)
 * - xmalloc.c - memory allocation with error checking (*xrealloc)
 * - construct-samples.c - construct sample data for the article
 * (compare_doubles, find_slot)
 *
 *
 * ------------------------------------------------------------------------
 * A verbatim copy of Jochen Voss & Andreas Voss's copyright declaration.
 * ------------------------------------------------------------------------
 * Copyright (C) 2006  Jochen Voss, Andreas Voss.
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

#include <ggdmc.hpp>

/*---------------------------------------------------------------------------
 Original found in pde.c
 ---------------------------------------------------------------------------*/

// xmalloc.c to fixed the memory location, but dynamically allocated memory
// space in solve_tridiag.
static void *xrealloc (void *ptr, size_t newsize)
{
  if (newsize == 0)
  {
    if (ptr)  free (ptr);
    return NULL;
  }

  ptr = ptr ? realloc (ptr, newsize) : malloc(newsize);
  if (ptr == NULL)  Rcpp::stop("memory exhausted");
  return  ptr;
}

static void solve_tridiag(int n, const double *rhs, double *res, double left,
                double mid, double right)
/* Solve an n by n tridiagonal system of linear equations.
 *
 * The matrix has 'mid' on the diagonal, 'left' on the subdiagonal and
 * 'right' on the superdiagonal.
 */
{
    // These two statics fix identical memory locations for the variables,
    // (1) pointer-to-double *tmp and (2) int tmp_len, when this function is
    // terminated and the program is still running. That is, when calculating
    // a CDF, the procedure of solving its PDE will be called several times.

    // The function requests n-1 double's, so fix the location of a pointer,
    // instead of n-1 double.
    static double *tmp = NULL;
    static int tmp_len = 0;
    double p, old_res, old_tmp;
    int  i;

    // Show to yourself that they are indeed at the same location while been
    // called repeatedly.
    // std::cout << tmp << " = &tmp, " << &tmp_len << " = &tmp_len\n";

    if (n-1 > tmp_len) {
      /* Reallocating/freeing the scratch buffer for every call to
      * 'solve_tridiag' caused about 10% of the total CPU load during
      * some fast-dm runs.  To avoid this problem, re-use the
      * same buffer between runs if possible. */
      tmp = xrenew(double, tmp, n-1); //
      // still we want to scratch the content of
      tmp_len = n-1;
      // the two variables as explained by V&V's note above.
    }

    /* step 1: solving forward */
    tmp[0] = old_tmp = right / mid;
    res[0] = old_res = rhs[0] / mid;
    for (i=1; i<n-1; ++i) {
      p = 1.0/(mid-left*old_tmp);
      res[i] = old_res = (rhs[i] - left*old_res)*p;
      tmp[i] = old_tmp = right*p;
    }
    p = 1.0/(mid-left*old_tmp);
    res[n-1] = (rhs[n-1] - left*old_res)*p;

    /* step 2: solving backward */
    for (i=n-1; i>0; --i)  res[i-1] -= tmp[i-1]*res[i];

    // Does this work? No it does not
    // delete [] tmp;

}

static void make_step (int N, double *vector, double dt, double dz, double v)
/* Advance the numerical solution of the PDE by one step in time,
* using the Crank-Nicolson scheme.  The time step size is 'dt', the
* space grid size is 'dz'.  */
{

    double * tmp_vector = new double[N+1];
    double  left, mid, right;
    int	i;

    left  =  (1-dz*v) / (2*dz*dz);
    mid   =  -1 / (dz*dz);
    right =  (1+dz*v) / (2*dz*dz);

    tmp_vector[1] = (dt*left * vector[0] +
      (1+0.5*dt*mid) * vector[1] +
      0.5*dt*right * vector[2]);

    for (i=2; i<N-1; i++) {
      tmp_vector[i] = (0.5*dt*left * vector[i-1] +
        (1+0.5*dt*mid) * vector[i] +
        0.5*dt*right * vector[i+1]);
    }

    tmp_vector[N-1] = (0.5*dt*left * vector[N-2] +
      (1+0.5*dt*mid) * vector[N-1] +
      dt*right * vector[N]);

    solve_tridiag(N-1, tmp_vector+1, vector+1,
                  -0.5*dt*left, 1-0.5*dt*mid, -0.5*dt*right);

    delete [] tmp_vector;
}

static void advance_to (int N, double *vector, double t0, double t1, double dz,
              double v, double TUNE_PDE_DT_MIN, double TUNE_PDE_DT_MAX,
              double TUNE_PDE_DT_SCALE)
/* Advance the state 'vector' of the PDE from time 't0' to time 't1' */
{
    bool done = false;
    do {
      double  dt = TUNE_PDE_DT_MIN + TUNE_PDE_DT_SCALE*t0;
      if (dt > TUNE_PDE_DT_MAX)  dt = TUNE_PDE_DT_MAX;
      if (t0 + dt >= t1) {
        dt = t1 - t0;
        done = 1;
      }
      make_step (N, vector, dt, dz, v);
      t0 += dt;
    } while (! done);
}

/*---------------------------------------------------------------------------
 Original found in cdf.c
 ---------------------------------------------------------------------------*/
// Four different destructors
static void F_plain_delete (F_calculator *fc);
static void F_sz_delete    (F_calculator *fc);
static void F_sv_delete    (F_calculator *fc);
static void F_st0_delete   (F_calculator *fc);

static const double * F_plain_get_F (F_calculator *fc, double t);
static const double * F_sz_get_F    (F_calculator *fc, double t);
static const double * F_sv_get_F    (F_calculator *fc, double t);
static const double * F_st0_get_F   (F_calculator *fc, double t);

static double F_plain_get_z (const F_calculator *fc, int i);
static double F_sz_get_z    (const F_calculator *fc, int i);
static double F_sv_get_z    (const F_calculator *fc, int i);
static double F_st0_get_z   (const F_calculator *fc, int i);

// Use by sampling
static void F_plain_start(F_calculator *fc, int plus);
static void F_sz_start   (F_calculator *fc, int plus);
static void F_sv_start   (F_calculator *fc, int plus);
static void F_st0_start  (F_calculator *fc, int plus);

/*** plain: no variability **************************************************/
double F_plain_data::F_limit(double z)
{
  if (std::fabs(v) < 1e-8) { return 1 - z/a; }
  else {return ( std::exp(-2*v*z)-std::exp(-2*v*a) ) / (1-std::exp(-2*v*a) );}
}

static F_calculator *F_plain_new (Parameters * params)
{
  F_calculator * fc   = new F_calculator;
  F_plain_data * data = new F_plain_data;

  // V&V's note "N must be even, otherwise the case szr == 1 fails"
  int N = 2*(int)(params->a*0.5/params->TUNE_DZ+0.5);
  if (N<4) N = 4;

  fc->N = N;
  fc->plus = -1;

  data->F  = new double[N+1];

  data->a  = params->a;
  data->v  = params->v;
  data->t0 = params->t0;
  data->d  = params->d;

  data->TUNE_PDE_DT_MIN   = params->TUNE_PDE_DT_MIN;
  data->TUNE_PDE_DT_MAX   = params->TUNE_PDE_DT_MAX;
  data->TUNE_PDE_DT_SCALE = params->TUNE_PDE_DT_SCALE;

  data->dz = params->a/N;

  fc->data  = data;
  fc->start = F_plain_start;
  fc->free  = F_plain_delete;
  fc->get_F = F_plain_get_F;
  fc->get_z = F_plain_get_z;

  return  fc;
}

static const double *F_plain_get_F (F_calculator *fc, double t)
{
  F_plain_data *data = (F_plain_data*)fc->data;
  t -= data->t_offset;

  if (t > data->t)
  {
    advance_to (fc->N, data->F, data->t, t, data->dz, data->v,
                data->TUNE_PDE_DT_MIN, data->TUNE_PDE_DT_MAX,
                data->TUNE_PDE_DT_SCALE);
    data->t = t;
  }

  return data->F;
}

static double F_plain_get_z (const F_calculator *fc, int i)
  // Used by F_sz because data must be casted to appropriate type.
  // Syntax will not look nice if remove this function
{
  F_plain_data *data = (F_plain_data*)fc->data;
  return i * data->dz;
}

static void F_plain_delete (F_calculator *fc)
{
  F_plain_data *data = (F_plain_data*)fc->data;
  delete [] data->F;
  delete data;
  delete fc;
}

static void F_plain_start (F_calculator *fc, int plus)
{
  F_plain_data *data = (F_plain_data*)fc->data;
  int  N = fc->N;

  fc->plus = plus;
  data->t_offset = data->t0 - data->d * (plus == 1? 0.5 : -0.5);
  data->t = 0;

  data->F[0] = (plus == 1) ? 1 : 0;
  double z;
  for (int i=1; i<N; i++)
  {
    z = i*data->dz;   // original F_get_z (fc, i);
    data->F[i] = data->F_limit(z);
  }
  data->F[N] = (plus == 1) ? 1 : 0;
}


/*  sz ********************************************** */
static F_calculator *F_sz_new (Parameters *params)
{ // 3 news; 1 array;
  F_calculator * base_fc = F_plain_new (params);
  F_calculator * fc = new F_calculator;
  F_sz_data * data = new F_sz_data;

  double sz, tmp, dz;
  int N, k;

  sz = params->szr*params->a;
  if (sz < params->TUNE_SZ_EPSILON) return base_fc;

  N = base_fc->N;
  dz = F_plain_get_z(base_fc, 1) - F_plain_get_z(base_fc, 0);
  tmp = sz/(2*dz);
  k = (int)(std::ceil(tmp) + 0.5); // step in z space
  if (2*k > N) Rcpp::stop ("2*k > N"); // assert is slient in R

  fc->N = N-2*k;
  fc->plus = -1;
  data->avg = new double[fc->N+1];

  data->base_fc = base_fc;
  data->k = k;
  data->q = k - tmp;
  data->f = dz/sz;

  fc->data = data;
  fc->start = F_sz_start;
  fc->free = F_sz_delete;
  fc->get_F = F_sz_get_F;
  fc->get_z = F_sz_get_z;

  return  fc;
}

static void F_sz_delete (F_calculator *fc)
{
  F_sz_data *data = (F_sz_data *)fc->data;
  F_delete(data->base_fc);
  delete [] data->avg;
  delete data;
  delete fc;
}

static const double *F_sz_get_F (F_calculator *fc, double t)
{
  F_sz_data *data = (F_sz_data *)fc->data;

  const double * F;
  double  tmp, q, f;
  int  i, j, m;

  F = F_get_F(data->base_fc, t);

  m = 2*data->k;
  q = data->q;
  f = data->f;

  if (m >= 3)
  {
    for (i=0; i<=fc->N; ++i)
    {
      tmp  = F[i] * 0.5*(1-q)*(1-q);
      tmp += F[i+1] * (1-0.5*q*q);
      for (j=i+2; j<i+m-1; ++j)  tmp += F[j];
      tmp += F[i+m-1] * (1-0.5*q*q);
      tmp += F[i+m] * 0.5*(1-q)*(1-q);
      data->avg[i] = tmp * f;
    }
  } else {
    /* m == 2 */
    for (i=0; i<=fc->N; ++i)
    {
      tmp = F[i] * 0.5*(1-q)*(1-q);
      tmp += F[i+1] * (1-q*q);
      tmp += F[i+2] * 0.5*(1-q)*(1-q);
      data->avg[i] = tmp * f;
    }
  }
  /* m == 1 is impossible here */

  return data->avg;
}

static double F_sz_get_z (const F_calculator *fc, int i)
{
  F_sz_data *data = (F_sz_data *)fc->data;
  return  F_get_z (data->base_fc, i+data->k);
}

static void F_sz_start (F_calculator *fc, int plus)
{
  F_sz_data *data = (F_sz_data *)fc->data;
  fc->plus = plus;
  F_start (data->base_fc, plus);
}

/*** sv **************************************************/
static F_calculator * F_sv_new (Parameters *params)
{ // 4 new's; two array
  if (params->sv < params->TUNE_SV_EPSILON) return F_sz_new (params);

  F_calculator * fc = new F_calculator;
  F_sv_data * data  = new F_sv_data;

  int nv,j;
  double x;

  nv = (int)(params->sv/params->TUNE_DV + 0.5);
  if (nv < 3) nv = 3;

  // Create a temp copy of the parameters
  std::vector<F_calculator*> base_fc(nv);
  // F_calculator * (*base_fc) = new F_calculator*[nv];

  Parameters temp_params = *params;
  temp_params.sv = 0;   // Integrate across svs

  for (j=0; j<nv; ++j)
  {
    x = R::qnorm ((0.5+j)/nv, 0., 1., 1, 0);
    temp_params.v = params->sv*x + params->v;
    base_fc[j] = F_sz_new(&temp_params);
  }

  fc->N = base_fc[0]->N;
  fc->plus = -1;

  data->avg = new double[fc->N+1];
  data->nv = nv;
  data->base_fc = base_fc;

  fc->data  = data;
  fc->start = F_sv_start;
  fc->free  = F_sv_delete;
  fc->get_F = F_sv_get_F;
  fc->get_z = F_sv_get_z;

  return  fc;
}

static const double * F_sv_get_F (F_calculator *fc, double t)
{
  F_sv_data * data = (F_sv_data *)fc->data;
  const double *F;
  double *avg = data->avg;
  int  i, j;

  // Calculate avg in the F_sz_data
  F = F_get_F(data->base_fc[0], t);

  // retrieve avg in the F_sz_data
  for (i=0; i<=fc->N; ++i) avg[i] = F[i];

  for (j=1; j<data->nv; ++j)
  {
    F = F_get_F(data->base_fc[j], t);
    for (i=0; i<=fc->N; ++i)  avg[i] += F[i];
  }
  for (i=0; i<=fc->N; ++i)  avg[i] /= data->nv;


  return avg;
}

static void F_sv_delete (F_calculator *fc)
{
  F_sv_data *data = (F_sv_data *)fc->data;
  for (int j=0; j<data->nv; ++j) F_delete(data->base_fc[j]);
  delete [] data->avg;
  delete data;
  delete fc;
}

static double F_sv_get_z (const F_calculator *fc, int i)
{
  F_sv_data *data = (F_sv_data *)fc->data;
  return  F_get_z (data->base_fc[0], i);
}

static void F_sv_start (F_calculator *fc, int plus)
{
  F_sv_data *data = (F_sv_data *)fc->data;
  int  j;
  fc->plus = plus;
  // To set plus in each of the F_sz base_fc
  for (j=0; j<data->nv; ++j) F_start (data->base_fc[j], plus);
}

/*** st0 **************************************************/
static const double * F_st0_get_row(const F_calculator * fc, int j)
{
  const F_st0_data *data = (F_st0_data *)fc->data;
  int  M, N, idx;
  double *row;

  M = data->M;
  N = fc->N;
  if (0 > j || j > M) Rcpp::stop("j not in 0 ~ M (inclusive)");

  idx = (data->base + j)%M;
  row = data->values + idx*(N+1);

  if (! data->valid[idx])
  {
    double t;
    const double * F;
    t = data->start + j*data->dt;
    F = F_get_F(data->base_fc, t);

    memcpy(row, F, (N+1)*sizeof(double));
    data->valid[idx] = 1;
  }

  return row;
}

static void add_vec(long n, double a, const double * x, double * y)
{
#ifdef HAVE_LIBBLAS
  // Rcpp::Rcout <<"HAVE_LIBBLAS " << std::endl;
  extern void daxpy_(long *Np, double *DAp, const double *X, long *INCXp,
                     double *Y, long *INCYp);
  long inc = 1;
  daxpy_(&n, &a, x, &inc, y, &inc);
#else /* ! HAVE_LIBBLAS */
  // Rcpp::Rcout <<"No LIBBLAS " << std::endl;
  int i;
  if (a == 1) { for (i=0; i<n; ++i) y[i] += x[i]; }
  else { for (i=0; i<n; ++i) y[i] += a*x[i]; }
#endif /* ! HAVE_LIBBLAS */
}

static F_calculator * F_st0_new (Parameters * params)
{ // 5 new's here; two are individual pointers and 3 are arrays
  // 1 base_fc is from a previous layer
  F_calculator * base_fc = F_sv_new (params);
  if (params->st0 <= params->TUNE_DT0*1e-6) return base_fc;

  F_calculator * fc = new F_calculator;
  F_st0_data * data = new F_st0_data;

  int M = (int)(params->st0/params->TUNE_DT0 + 1.5);
  if (M < 3)  M = 3;
  int N = base_fc->N;

  data->st0     = params->st0;
  data->dt      = params->st0/(M-2);
  data->base_fc = base_fc;
  data->M       = M;
  data->base   = 0;
  data->values = new double[M*(N+1)];
  data->valid  = new char[M];
  data->avg    = new double[N+1];

  fc->N    = N;
  fc->plus = -1;
  fc->data = data;

  fc->start = F_st0_start;
  fc->free  = F_st0_delete;
  fc->get_F = F_st0_get_F;
  fc->get_z = F_st0_get_z;

  return  fc;
}

static const double * F_st0_get_F (F_calculator *fc, double t)
{
  F_st0_data *data = (F_st0_data*)fc->data;
  double  a, b;
  const double *row;
  double  q, r, *avg;
  int  N, shift;
  int  i, j, m;

  a = t - 0.5*data->st0;
  b = t + 0.5*data->st0;
  N = fc->N;

  // how many of the precalculated rows can we keep?
  if (a - data->start >= data->M*data->dt)
  {
    shift = data->M; // beware of integer overflows for small dt
  } else {
    shift = (int)((a - data->start)/data->dt);
    if (shift < 0) Rcpp::stop("shift < 0; F_st0_get_F problem");
  }

  for (j=0; j<shift; ++j) data->valid[(data->base+j)%data->M] = 0;

  if (shift < data->M)
  {
    data->start += shift*data->dt;
    data->base = (data->base+shift)%data->M;
  } else {
    data->start = a;
  }

  /* compute the average over the rows from a to b */
  avg = data->avg;
  for (i=0; i<=N; ++i)  avg[i] = 0;

  // A independent block to confine scope
  { // tmp allocated and scope begins
    double tmp = (b - data->start)/data->dt;
    m = (int)(std::ceil(tmp) + 0.5);      // V&V's note:
    if (m >= data->M)  m = data->M-1; // "protect against rounding errors"
    q = (a - data->start)/data->dt;
    r = m - tmp;
  } // tmp expires


  if (m >= 3)
  {
    row = F_st0_get_row(fc, 0);
    add_vec(N+1, 0.5*(1-q)*(1-q), row, avg);

    row = F_st0_get_row(fc, 1);
    add_vec(N+1, 1-0.5*q*q, row, avg);

    for (j=2; j<m-1; ++j)
    {
      row = F_st0_get_row(fc, j);
      add_vec(N+1, 1, row, avg);
    }

    row = F_st0_get_row(fc, m-1);
    add_vec(N+1, 1-0.5*r*r, row, avg);

    row = F_st0_get_row(fc, m);
    add_vec(N+1, 0.5*(1-r)*(1-r), row, avg);
  }
  else if (m == 2)
  {
    row = F_st0_get_row(fc, 0);
    add_vec(N+1, 0.5*(1-q)*(1-q), row, avg);

    row = F_st0_get_row(fc, 1);
    add_vec(N+1, 1-0.5*(q*q+r*r), row, avg);

    row = F_st0_get_row(fc, 2);
    add_vec(N+1, 0.5*(1-r)*(1-r), row, avg);
  }
  else if (m == 1)
  {
    row = F_st0_get_row(fc, 0);
    add_vec(N+1, 0.5*((1-q)*(1-q)-r*r), row, avg);

    row = F_st0_get_row(fc, 1);
    add_vec(N+1, 0.5*((1-r)*(1-r)-q*q), row, avg);
  }

  for (i=0; i<=N; ++i) avg[i] *= data->dt/(b-a);
  return avg;
}

static void F_st0_delete (F_calculator *fc)
{
  F_st0_data *data = (F_st0_data *)fc->data;
  F_delete (data->base_fc); // base_fc is a fc created by F_sv_new
  delete [] data->valid;
  delete [] data->values;
  delete [] data->avg;
  delete data;
  delete fc;
}

static double F_st0_get_z (const F_calculator *fc, int i)
{
  F_st0_data *data = (F_st0_data*)fc->data;
  return  F_get_z (data->base_fc, i);
}

static void F_st0_start (F_calculator *fc, int plus)
{
  F_st0_data *data = (F_st0_data *)fc->data;
  int j;

  fc->plus = plus;
  F_start (data->base_fc, plus);
  data->start = -DBL_MAX;

  // initially mark all of the cache as invalid
  for (j = 0; j < data->M; ++j) data->valid[j] = 0;
}

// The entry point  ----------------------------------------------------
// Change F_st0_new to eg F_sv_new, F_sz_new, or F_plain_new to test code from
// a different orion layer.
F_calculator * F_new (Parameters *params)
{
  return F_st0_new (params);
}

void F_start (F_calculator *fc, int boundary)
  // Set the initial condition for the PDE.
  // If upper boundary prepare to calculate the CDF for hitting a before 0,
  // otherwise prepare to calculate the CDF for hitting 0 before a.
  //
  // CONVERSION NOTE: Changed from enum boundary b to int boundary, where 0 =
  // lower and 1 = upper
{
  fc->start (fc, boundary);
}

const double *F_get_F (F_calculator *fc, double t)
{
  return fc->get_F (fc, t);
}

void F_delete (F_calculator *fc) { fc->free (fc); }

double F_get_z (const F_calculator *fc, int i)
  // Get the z-value corresponding to index i.
{
  return  fc->get_z (fc, i);
}

double F_get_val (F_calculator *fc, double t, double z)
  // Get the value of the CDF for the parameters given when creating fc.
  // The function uses linear interpolation for z-values between the grid points.
  // Don't use this function for parameter fitting, since it is not very fast
  // (use 'F_get_F' instead).
  //
  // ? CONVERSION NOTE: Speed difference seems to be within +/- 10% for
  // rt/parameter combinations so far tested
  //
{
  const double *F;
  double   z0, z1;
  double     p, x;
  int  N = fc->N;
  int  i;

  F = F_get_F (fc, t);

  if (N == 0)
  {
    x = F[0];
  }
  else
  {
    z0 = F_get_z (fc, 0);
    z1 = F_get_z (fc, N);
    i = (int)(N*(z-z0)/(z1-z0));

    if (i < N)
    {
      z0 = F_get_z (fc, i);
      z1 = F_get_z (fc, i+1);
      p = (z1-z) / (z1-z0);
      x = p*F[i] + (1-p)*F[i+1];
    }
    else
    {
      x = F[N];
    }
  }
  return  x;
}

/*---------------------------------------------------------------------------
 Original found in construct-samples.c
 ---------------------------------------------------------------------------*/
static int compare_doubles (const void *a, const void *b)
{
  const double *da = (double *)a;
  const double *db = (double *)b;

  if (*da < *db)  return -1;
  if (*da > *db)  return 1;
  return  0;
}

static int find_slot(double target, const double *value, int l, int r)
{
  int m = (l+r)/2;
  if (m == l) { return l; }
  else if ( value[m] > target) { return find_slot(target, value, l, m); }
  else    { return find_slot(target, value, m, r); }
}

/*---------------------------------------------------------------------------
 Original found in rtdists 0.9-0
 ---------------------------------------------------------------------------*/

Rcpp::List sampling(int s_size, Parameters * params, bool random_flag)
{
  double Fs_min=1, Fs_max=0, t_min=-.5, t_max=.5, dt, t;
  int  i, N;

  /* get the F-values for the samples */
  double * Fs = new double[s_size];

  if (random_flag)
  {
    for (i=0; i<s_size; ++i)
    {
      Fs[i] = Rf_runif(0, 1);
      if (Fs[i] > Fs_max) Fs_max = Fs[i];
      if (Fs[i] < Fs_min) Fs_min = Fs[i];
    }
  }
  else
  {
    // Generate equally-spaced samples
    for (i=0; i<s_size; ++i) { Fs[i] = (i+0.5)/s_size; }
    Fs_min = Fs[0];
    Fs_max = Fs[s_size-1];
  }

  /////////////////////////////////////////////////////////////////////////////
  // Create a new F calculator
  F_calculator *fc = F_new(params);
  double z = params->zr*params->a;

  // get the required t-range; BOUNDARY_UPPER & BOUNDARY_LOWER in Parameters.h
  F_start (fc, BOUNDARY_UPPER);
  while ( F_get_val(fc, t_max, z) < Fs_max )	t_max += 0.1;

  F_start (fc, BOUNDARY_LOWER);
  while ( F_get_val (fc, -t_min, z) > Fs_min)	t_min -= 0.1;

  N = (int)((t_max-t_min)/0.001 + 0.5); // get a table of F-values
  dt = (t_max-t_min)/N;
  double * F = new double[N+1];

  F_start(fc, BOUNDARY_UPPER);
  for (i=0; i<=N; ++i)
  {
    t = t_min+i*dt;
    if (t < 0) continue;
    F[i] = F_get_val (fc, t, z);
  }

  F_start (fc, BOUNDARY_LOWER);
  for (i=N; i>=0; --i)
  {
    t = -(t_min+i*dt);
    if (t < 0) continue;
    F[i] = F_get_val (fc, t, z);
  }
  F_delete (fc);

  // protect against rounding errors: make F increasing and restrict to valid range
  for (i=0; i<=N; ++i)
  {
    if (F[i] < 0)	F[i] = 0;
    if (F[i] > 1)	F[i] = 1;
  }

  std::qsort(F, N+1, sizeof(double), compare_doubles); // use quick sort in cstdlib
  if (F[0] > Fs_min)		F[0] = Fs_min;
  if (F[N] < Fs_max)		F[N] = Fs_max;

  std::vector<double> out_RTs(s_size);
  std::vector<unsigned int> out_bounds(s_size);

  for (i=0; i<s_size; ++i)
  {
    double y = Fs[i];
    int k = find_slot(y, F, 0, N);
    t = t_min + (k + (y-F[k]) / (F[k+1]-F[k]))*dt;

    if (F[k] > y || y > F[k+1]) Rcpp::stop("y not in the range");

    out_bounds[i] = t >= 0;
    out_RTs[i] = std::fabs(t);
  }

  delete [] F;
  delete [] Fs;

  return Rcpp::List::create(Rcpp::Named("rt") = out_RTs,
                            Rcpp::Named("boundary") = out_bounds);
}


// [[Rcpp::export]]
Rcpp::List r_fastdm (unsigned int num_values, std::vector<double> params,
                     double precision=3, bool stop_on_error=true)
// R-style sampling from the DM - returns a List consisting of RTs and
// boundaries

{
  if ((num_values < 1) || (num_values > MAX_INPUT_VALUES))
  {
    Rcpp::stop("Number of samples requested exceeds maximum of %d.\n", MAX_INPUT_VALUES);
  }

  Parameters * g_Params = new Parameters (params, precision);

  if (!g_Params->ValidateParams(stop_on_error))
  {
    if (stop_on_error)
    {
      Rcpp::stop("Error validating parameters.\n");
    }
    else
    {
      Rcpp::NumericVector out_RTs(num_values, 0.0);
      Rcpp::NumericVector out_bounds(num_values, 0.0);
      return Rcpp::List::create(Rcpp::Named("rt") = out_RTs,
                                Rcpp::Named("boundary") = out_bounds);
    }
  }

  // Pass through to Sampling.hpp
  Rcpp::List out = sampling (num_values, g_Params, true);

  delete g_Params;

  return out;
}
