/* Copyright (C) 2019  Yi-Shin Lin
 *  ------------------------------------------------------------------------
 *  The algorithm is based on the paper:
 *  Blurton, S., P., Kesselmeier, M., & Gondan, M. (2017). The first-passage
 *  time distribution for the diffusion model with variable drift. Journal of
 *  Mathematical Psychology, doi: 10.1016/j.jmp.2016.11.003
 *  ------------------------------------------------------------------------
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

static double integrate_over_zr (double x, void * params)
{
  CDF *cdf = (CDF *) params;
  double out = cdf->DT < 0 ? 0 : cdf->G0_over_zr(cdf->DT, x) *
    R::dunif(x, cdf->lower, cdf->upper, 0);
  return out;
}

static double integrate_over_t (double x, void * params)
{
  CDF *cdf = (CDF *) params;
  double out;

  if (x < 0) {
    out = 0;
  } else {
    double term0 = cdf->G0(x);
    double term1 = R::dunif(x, cdf->lower, cdf->upper, 0);
    out = term0 * term1;
  }

  return out;
}

static double integrate_over_t_and_zr (double x, void * params)
{
  CDF *cdf0 = (CDF *) params;
  CDF cdf1(*cdf0);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double over_zr, error, out;
  gsl_function F;

  cdf1.lower = cdf1.mean_zr - .5*cdf1.szr; // lower
  cdf1.upper = cdf1.mean_zr + .5*cdf1.szr; // upper
  cdf1.DT    = x;

  F.function = &integrate_over_zr;
  F.params   = &cdf1;

  gsl_integration_qags (&F, cdf1.lower, cdf1.upper, 0, 1e-10, 1000, w,
                        &over_zr, &error);

  out = x < 0 ? 0 : over_zr * R::dunif(x, cdf0->lower, cdf0->upper, 0);

  gsl_integration_workspace_free (w);
  return out;
}

static void distribution (std::vector<double> rts, Parameters * params,
                          std::vector<double> & out)
{
  unsigned int nrt = rts.size();
  CDF * obj0 = new CDF (params);

  gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
  double result, error;
  gsl_function F;

  if (params->szr != 0 && params->st0 != 0)  {
    F.function = &integrate_over_t_and_zr;
    F.params = obj0; // F.params is void type

    for(size_t i=0; i<nrt; i++)
    {
      obj0->DT    = rts[i] - params->t0 - .5*params->d;
      obj0->lower = obj0->DT - .5*params->st0;
      obj0->upper = obj0->DT + .5*params->st0;

      gsl_integration_qags (&F, obj0->lower, obj0->upper, 0, 1e-10, 1000, w,
                            &result, &error);
      out[i] = result;
    }

  } else if (params->szr != 0 && params->st0 == 0) {
    obj0->lower = params->zr - params->szr/2; // lower
    obj0->upper = params->zr + params->szr/2; // upper

    F.function = &integrate_over_zr;
    F.params = obj0;

    for(size_t i=0; i<nrt; i++)
    {
      obj0->DT = rts[i] - params->t0 - .5*params->d;
      gsl_integration_qags (&F, obj0->lower, obj0->upper, 0, 1e-10, 1000, w,
                            &result, &error);
      out[i] = result;
    }


  } else if (params->szr == 0 && params->st0 != 0) {
    F.function = &integrate_over_t;
    F.params = obj0;

    for(size_t i=0; i<nrt; i++)
    {
      obj0->DT = rts[i] - params->t0 - 0.5*params->d;
      obj0->lower = obj0->DT - .5*params->st0;
      obj0->upper = obj0->DT + .5*params->st0;

      gsl_integration_qags (&F, obj0->lower, obj0->upper, 0, 1e-10, 1000, w,
                            &result, &error);
      out[i] = result;
    }

  } else if (params->st0 == 0 && params->szr == 0) {

    for(size_t i=0; i<nrt; i++)
    {
      obj0->DT = rts[i] - params->t0 - .5*params->d;
      out[i] = obj0->G0(obj0->DT);
    }

  } else {
    Rcout << "Unexpected condition \n";
  }

  gsl_integration_workspace_free (w);
  delete obj0;
}

// [[Rcpp::export]]
std::vector<double> prd (std::vector<double> rts,
                         std::vector<double> params,
                         double precision = 3,
                         unsigned int b = 2)
{
  // - CDF which finds the left-hand area (from 0 to RT)
  // - pass boundary to retrieve (1 = lower, 2 = upper)
  unsigned int nrt = rts.size();
  std::vector<double> out(nrt);

  if (nrt > MAX_INPUT_VALUES) Rcpp::stop("Allows only %d RTs.\n", MAX_INPUT_VALUES);
  if ((b < 1) || (b > 2)) Rcpp::stop ("Boundary must be either 2 (upper) or 1 (lower)\n");

  // Boundary selection has been done in the Parameters constructor
  // precision is to be compatible with F_calculator.
  Parameters * g_Params = new Parameters (params, precision, b-1);

  if (!g_Params->ValidateParams(true)) std::fill(out.begin(), out.end(), 0);
  else distribution (rts, g_Params, out);

  delete g_Params;
  return out;
}

/* Expose these functions for unit tests */

// std::vector<double> logMill(std::vector<double> x)
// {
//   std::vector<double> out(x.size());
//   for (size_t i=0; i<x.size(); i++)
//   {
//     if (x[i] >= 1e4) { out[i] = -std::log(x[i]); }
//     else             { out[i] = R::pnorm(x[i], 0, 1, 0, 1) - R::dnorm(x[i], 0, 1, 1); }
//   }
//   return out;
// }
// std::vector<double> test_G0(std::vector<double> params,
//                             std::vector<double> t)
// {
//   unsigned int nt = t.size();
//   std::vector<double> out(nt);
//   CDF * obj0 = new CDF (params);
//
//   for(size_t i=0; i<nt; i++) out[i] = obj0->G0(t[i] - params[6]);
//
//   delete obj0;
//   return out;
// }
// std::vector<double> test_Ga(std::vector<double> params,
//                             std::vector<double> t)
// {
//   unsigned int nt = t.size();
//   std::vector<double> out(nt);
//   std::vector<double> params_ = params;
//   params_[1] = -params[1];
//   params_[2] = 1-params[2];
//
//   CDF * obj0 = new CDF (params_);
//
//   for(size_t i=0; i<nt; i++) out[i] = obj0->G0(t[i] - params[6]);
//
//   delete obj0;
//   return out;
// }
