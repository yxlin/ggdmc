/* The is modified from the C version of the constant drift-diffusion density
 * function in fast-dm 30.2 density.c - compute the densities g- and g+ of the
 * first exit time. Copyright (C) 2012  Andreas Voss, Jochen Voss. */
#define _USE_MATH_DEFINES // For Windows users
#include <RcppArmadillo.h>
// #include <cmath>
// #include <vector>

void set_precision (double p) ;

// Change Voss's C-style struct function to C++ STL vector
double integrate(double (*F)(double, std::vector<double>&),
  std::vector<double>& pVec,
  double a, double b, double step_width) ;

double integrate_parallel(double (*F)(double, std::vector<double>&),
  std::vector<double>& pVec, double a, double b,
  double step_width) ;

/* Formula A3; Note t here is DT; They pass t/a^2 to get zr. That is, z/a not
A3's z; i is A3's n */
double g_minus_small_time(double DT, double zr, int N) ;

// Formula A4; ; Note t here is DT
double g_minus_large_time(double DT, double zr, int N) ;

double g_minus_no_var(double DT, double a, double zr, double v) ;

// When sv is not 0, change the "factor" to incorporate it.
// integral_v_g_minus will choose either g_minus_no_var or itself
double integral_v_g_minus(double zr, std::vector<double>& pVec) ;

double integral_z_g_minus(double DT, std::vector<double>& pVec) ;
double integral_z_g_minus_parallel(double DT, std::vector<double>& pVec) ;

double integral_t0_g_minus(double DT, std::vector<double>& pVec) ;
double integral_t0_g_minus_parallel(double DT, std::vector<double>& pVec) ;

double g_minus(std::vector<double> pVec) ;
double g_plus(std::vector<double> pVec) ;

double g_minus_parallel(std::vector<double> pVec) ;
double g_plus_parallel(std::vector<double> pVec) ;

