#include <ggdmc.hpp>

/*---------------------------------------------------------------------------
 Original public functions in Parameters.h. The version with std::vector is
 an overloaded constructor. Show function is for debugging purpose.
 ---------------------------------------------------------------------------*/
Parameters::Parameters(std::vector<double> params, double precision)
{
  a   = params[PARAM_a];
  v   = params[PARAM_v];
  zr  = params[PARAM_zr];
  d   = params[PARAM_d];
  szr = params[PARAM_szr];
  sv  = params[PARAM_sv];
  t0  = params[PARAM_t0];
  st0 = params[PARAM_st0];

  SetPrecision (precision);
}

Parameters::Parameters(double * params, double precision)
{
  a   = params[PARAM_a];
  v   = params[PARAM_v];
  zr  = params[PARAM_zr];
  d   = params[PARAM_d];
  szr = params[PARAM_szr];
  sv  = params[PARAM_sv];
  t0  = params[PARAM_t0];
  st0 = params[PARAM_st0];

  SetPrecision (precision);
}

bool Parameters::ValidateParams (bool print)
{
  using namespace Rcpp;
  bool valid = true;

  if (a <= 0)                         { valid = false; if (print) Rcout << "error: invalid parameter a = " << a << std::endl;  }
  if (szr < 0 || szr > 1)             { valid = false; if (print) Rcout << "error: invalid parameter szr = " << szr << std::endl; }
  if (st0 < 0)                        { valid = false; if (print) Rcout << "error: invalid parameter st0 = " << st0 << std::endl; }
  if (sv < 0)                         { valid = false; if (print) Rcout << "error: invalid parameter sv = " << sv << std::endl; }
  if (t0 - fabs(0.5*d) - 0.5*st0 < 0) { valid = false; if (print) Rcout << "error: invalid parameter combination t0 = " << t0 << ", d = " << d << ", st0 =" << st0 << std::endl; }
  if (zr - 0.5*szr <= 0)              { valid = false; if (print) Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (zr + 0.5*szr >= 1)              { valid = false; if (print) Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}

  return valid;
}

void Parameters::Show(std::string str) const
{
  Rcout << str << ":\n";
  Rcout << "[a\tv\tt0\td]  = " << "[" << a << "\t" << v << "\t" << t0
            << "\t" << d << "]" << std::endl;
  Rcout << "[szr\tsv\tst0\tzr] = " << "[" << szr << "\t" << sv << "\t"
            << st0 << "\t" << zr << "]" << std::endl;
}

/*---------------------------------------------------------------------------
 static functions originally found in Density.h. They are confined within
 this file.
---------------------------------------------------------------------------*/
static double g_minus_small_time (double t, double zr, int N)
  // A3 Formula on p1217; Appendix Mathematical Details V & V (2004)
  // See also Feller (1971, p359 & p370 Problem 22); Note zr = z/a. v = 0 & a = 1
{
  if (t <= 0) Rcpp::stop("t must be greater than 0.");

  int i; // i must be int not unsigned int
  double d, s = 0;
  for (i = -N/2; i <= N/2; i++) // i must be int; it cannot be unsigned
  {
    d = 2*i + zr;
    s += exp(-d*d / (2*t)) * d;
  }

  return s / sqrt(6.283185307179586231996*t*t*t);
}

static double g_minus_large_time (double t, double zr, int N)
  // A4, p1217
{
  int i;
  double d, s = 0;
  for (i = 1; i <= N; i++)
  {
    d = i * M_PI;
    s += exp(-0.5*d*d*t) * sin(d*zr) * i;
  }
  return s * M_PI;
}

static double g_minus_no_var (double t, double a, double zr, double v)
  // Depending on g_minus_small_time and g_minus_large_time
  // Determine the refinement for the infinite serie
{
  int N_small, N_large;
  double simple, factor, eps, ta = t/(a*a);

  factor = exp(-a*zr*v - 0.5*(v*v)*t) / (a*a); // Front term in A3
  if (std::isinf(factor)) { return 0; }
  eps = EPSILON / factor;

  N_large = (int)ceil(1./(M_PI*sqrt(t)));

  if (M_PI*ta*eps < 1.) { // std::max as imax
    N_large = std::max(N_large,
                       (int)ceil(sqrt(-2.0*log(M_PI*ta*eps) / ((M_PI*M_PI)*ta))));
  }

  if (2.*sqrt(2.*M_PI*ta)*eps < 1.) { // std::max as fmax
    N_small = (int)ceil(std::max(sqrt(ta) + 1,
                        2. + sqrt(-2.0*ta*log(2.0*eps*sqrt(2.0*M_PI*ta)))));
  } else {
    N_small = 2;
  }

  if (N_small < N_large) {
    simple = g_minus_small_time(ta, zr, N_small); // Note it's ta not t
  } else {
    simple = g_minus_large_time(ta, zr, N_large);
  }
  return factor * simple;
}

static double integral_v_g_minus (double t, double zr, Parameters *params)
  // integrate over a range of v based on a formula; ie sv != 0
{
  int N_small, N_large;
  double a=params->a, v=params->v, sv=params->sv;
  double simple, factor, eps, ta = t/(a*a);

  // The factor is where difference is
  // Here uses a*zr to get z (zr = z/a). Must not change the multiplication
  // sequence of a*zr*a*zr*sv*sv; otherwise it will produce rounding errors
  factor = 1 / (a*a * sqrt(t * sv*sv + 1)) *
    exp(-0.5 * (v*v*t + 2*v*a*zr - a*zr*a*zr*sv*sv) / (t*sv*sv+1));

  // Early exit 1
  if (std::isinf(factor)) {return 0;}
  eps = EPSILON / factor;

  // Early exit 2
  if (sv == 0) {return g_minus_no_var(t, a, zr, v);}

  // Below is identical as in g_minus_no_var
  N_large = (int)ceil(1./(M_PI*sqrt(t)));
  if (M_PI*ta*eps < 1.) {
    N_large = std::max(N_large,
                       (int)ceil(sqrt(-2.*log(M_PI*ta*eps) / (M_PI*M_PI*ta))));
  }

  if (2.*sqrt(2.*M_PI*ta)*eps < 1.) {
    N_small = (int)ceil(std::max(sqrt(ta) + 1.,
                        2. + sqrt(-2.*ta*log(2.*eps*sqrt(2.*M_PI*ta)))));
  } else {
    N_small = 2;
  }

  if (N_small < N_large) {
    simple = g_minus_small_time(ta, zr, N_small);
  } else {
    simple = g_minus_large_time(ta, zr, N_large);
  }
  return factor * simple;
}

static double integrate_v_over_zr (Parameters *params, double a, double b,
                                   double t, double step_width)
  // used by integral_z_g_minus
{
  double x, s=0, width=b-a;
  int N = std::max(4, (int) (width / step_width)); // usually less than 10
  double step=width/N;
  for(x = a+0.5*step; x < b; x += step)
  {
    s += step * integral_v_g_minus (t, x, params); // width * height = area
  }
  return s;
}

static double integral_z_g_minus (double t, Parameters *params)
  // integral over a uniform, fixed range of zr, depending on
  // integrate_v_over_zr and integral_v_g_minus.
{
  double out;
  if (t <= 0) { return 0; };

  if (params->szr < params->TUNE_SZ_EPSILON) {
    out = integral_v_g_minus(t, params->zr, params);
  } else {
    out = integrate_v_over_zr(
      params,
      params->zr - .5*params->szr,
      params->zr + .5*params->szr,
      t, params->TUNE_INT_Z) / params->szr;
  }
  return out;
}

static double integrate_z_over_t (Parameters *params, double a, double b,
                                  double step_width)
{
  double x, s=0, width=b-a;
  int N = std::max(4, (int) (width / step_width));
  double step=width/N;
  for(x = a+0.5*step; x < b; x += step)
  {
    s += step * integral_z_g_minus(x, params);
  }
  return s;
}

static double integral_t0_g_minus (double t, Parameters *params)
  // integral over a uniform, fixed range of t, depending on integrate_z_over_t,
  // and integral_z_g_minus, which depends on integrate_v_over_zr and
  // integral_v_g_minus.
{
  double out;
  if (params->st0 < params->TUNE_ST0_EPSILON) {
    out = integral_z_g_minus(t, params);
  } else {
    out = integrate_z_over_t(
      params,
      t - .5*params->st0,
      t + .5*params->st0,
      params->TUNE_INT_T0) / params->st0;
  }
  return out;
}

/*---------------------------------------------------------------------------
 Non-static functions decleared in Parameters.hpp. Other functions can call
 them.
 ---------------------------------------------------------------------------*/
double g_minus(double t, Parameters *params)
{
  return integral_t0_g_minus (t - params->t0 - 0.5*params->d, params);
}

double g_plus(double t, Parameters *params)
{
  Parameters params_(*params);
  params_.zr = 1 - params->zr;
  params_.v = - params->v;
  return integral_t0_g_minus (t - params_.t0 + 0.5*params_.d, &params_);
}
