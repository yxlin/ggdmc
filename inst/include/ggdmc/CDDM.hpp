// CDDM.hpp
// version 03
#ifndef CDDM_HPP
#define CDDM_HPP

#include <RcppArmadillo.h>

using namespace Rcpp;

// GSL constants
#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define ROOT_EIGHT (2.0*M_SQRT2)

/* GSL routines: gsl_sf_bessel_zero_J0 & gsl_sf_bessel_J1 adapted to CDDM
 * in order to distribute the package for the Windows user
 *
 * The GSL routines license:
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* The following cheb_series is from specfunc/chebyshev.h
 * data for a Chebyshev series over a given interval
 */
struct cheb_series_struct {
  double * c;   /* coefficients                */
  int order;    /* order of expansion          */
  double a;     /* lower interval point        */
  double b;     /* upper interval point        */
  int order_sp; /* effective single precision order */
};
typedef struct cheb_series_struct cheb_series;


/*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

/* From specfunc/bessel_J1.c */

/* based on SLATEC besj1, 1983 version, w. fullerton */

/* chebyshev expansions

 series for bj1        on the interval  0.          to  1.60000d+01
 with weighted error   4.48e-17
 log weighted error  16.35
 significant figures required  15.77
 decimal places required  16.89

 */
static double bj1_data[12] = {
  -0.11726141513332787,
  -0.25361521830790640,
  0.050127080984469569,
  -0.004631514809625081,
  0.000247996229415914,
  -0.000008678948686278,
  0.000000214293917143,
  -0.000000003936093079,
  0.000000000055911823,
  -0.000000000000632761,
  0.000000000000005840,
  -0.000000000000000044,
};
static cheb_series bj1_cs = {
  bj1_data,
  11,
  -1, 1,
  8
};

/* The following two chebyshev series are from specfunc/bessel_amp_phase.c */
static double bm1_data[21] = {
  0.1047362510931285,
  0.00442443893702345,
  -0.00005661639504035,
  0.00000231349417339,
  -0.00000017377182007,
  0.00000001893209930,
  -0.00000000265416023,
  0.00000000044740209,
  -0.00000000008691795,
  0.00000000001891492,
  -0.00000000000451884,
  0.00000000000116765,
  -0.00000000000032265,
  0.00000000000009450,
  -0.00000000000002913,
  0.00000000000000939,
  -0.00000000000000315,
  0.00000000000000109,
  -0.00000000000000039,
  0.00000000000000014,
  -0.00000000000000005,
};

const cheb_series _gsl_sf_bessel_amp_phase_bm1_cs = {
  bm1_data,
  20,
  -1, 1,
  10
};

static double bth1_data[24] = {
  0.74060141026313850,
  -0.004571755659637690,
  0.000119818510964326,
  -0.000006964561891648,
  0.000000655495621447,
  -0.000000084066228945,
  0.000000013376886564,
  -0.000000002499565654,
  0.000000000529495100,
  -0.000000000124135944,
  0.000000000031656485,
  -0.000000000008668640,
  0.000000000002523758,
  -0.000000000000775085,
  0.000000000000249527,
  -0.000000000000083773,
  0.000000000000029205,
  -0.000000000000010534,
  0.000000000000003919,
  -0.000000000000001500,
  0.000000000000000589,
  -0.000000000000000237,
  0.000000000000000097,
  -0.000000000000000040,
};
const cheb_series _gsl_sf_bessel_amp_phase_bth1_cs = {
  bth1_data,
  23,
  -1, 1,
  12
};



class cddm
{
  public:
    double m_v1, m_v2, m_a, m_t0, m_sigma1, m_sigma2, m_eta1, m_eta2, m_tmax;
    double sigmasq, m_a2, m_twoa2, sigma2overa2, m_h, m_w;
    double m_eta1onsigma2, m_eta2onsigma2;
    unsigned int m_kmax, m_nw, m_sz;
    double m_sigma;
    // m_tmax,  /* Define the maximum time   */
    // m_nw,    /* Number of steps on circle */
    // m_sz     /* Number of time steps */
    arma::vec m_thetai;

    cddm(double v1, double v2, double a, double t0, double sigma1, double sigma2,
         double eta1, double eta2, double tmax, unsigned int kmax,
         unsigned int nw, unsigned int sz) :
      m_v1(v1), m_v2(v2), m_a(a), m_t0(t0), m_sigma1(sigma1), m_sigma2(sigma2),
      m_eta1(eta1), m_eta2(eta2), m_tmax(tmax), m_kmax(kmax), m_nw(nw), m_sz(sz)
    // dcircle
    {
      // Variables for d0_GSL
      sigmasq = m_sigma1*m_sigma1;
      m_a2    = m_a*m_a;
      m_twoa2 = 2.0*m_a2;
      sigma2overa2 = sigmasq / m_a2;
      m_h          = m_tmax / m_sz;    // time step size for d0_GSL

      // Variables for d function
      m_w = M_2PI / m_nw;
      if (m_eta1 < 1e-5) m_eta1 = 0.01; // guard agains NAN resulting from 0/0
      if (m_eta2 < 1e-5) m_eta2 = 0.01;

      /* Assume same diffusion in both directions */
      m_eta1onsigma2 = (m_eta1*m_eta1) / sigmasq;
      m_eta2onsigma2 = (m_eta2*m_eta2) / sigmasq;
    }

    cddm(double v1, double v2, double a, double t0, double sigma1, double sigma2,
         double tmax, double h, unsigned int nw) :
      m_v1(v1), m_v2(v2), m_a(a), m_t0(t0), m_sigma1(sigma1), m_sigma2(sigma2),
      m_tmax(tmax), m_h(h), m_nw(nw)
    // rcircle
    {
      if ( m_h <= 0  )  Rcpp::stop("h must be greater than 0.");
      if ( m_h >= 10)   Rcpp::stop("The time step is in second unit not in millisecond.");
      if (m_tmax <= 0)  Rcpp::stop("tmax must be greater than 0.");
      if (m_tmax < 1.0) Rcpp::Rcout << "tmax less than 1. Check tmax and h sequence\n";

      if (m_sigma1 < 0) Rcpp::stop("Drift rate SD for x axis must be greater than 0.");
      if (m_sigma2 < 0) Rcpp::stop("Drift rate SD for y axis must be greater than 0.");
      if ( m_t0 < 0  )  Rcpp::stop("t0 must be greater than 0.");

      m_sz = m_tmax/m_h; // nmax
      m_a2 = m_a*m_a;
      m_w  = M_2PI/m_nw; // m_w become inf when m_nw == 0 (NA will be casted to 0)
      m_thetai = arma::linspace(-M_PI, M_PI, m_nw+1); // +1 to circle back to starting

    }

    void d0_GSL(arma::vec & DT, arma::vec & Gt0)
    // First-passage-time density for Bessel process (GSL).
    // The function used GSL Library and generate PDF table.
    {
      double J0k[m_kmax], J0k2[m_kmax], J1k[m_kmax];
      unsigned int k, i;
      for (k = 0; k < m_kmax; k++)
      {
        J0k[k]  = gsl_sf_bessel_zero_J0(k+1);
        J0k2[k] = J0k[k] * J0k[k];
        J1k[k]  = gsl_sf_bessel_J1(J0k[k]);
      }

      DT[0] = 0;
      Gt0[0] = 0;
      for (i = 1; i < m_sz; i++)
      {
        DT[i] = i * m_h;
        Gt0[i] = 0;
        for (k = 0; k < m_kmax; k++)
        {
          Gt0[i] += J0k[k]*std::exp(-J0k2[k]*sigmasq*DT[i] / m_twoa2) / J1k[k];
        }
        Gt0[i] *= sigma2overa2;
      }
    }

    void d(arma::vec & DT, arma::vec & R, arma::vec & Gt0, arma::mat & Gt,
           arma::vec & Ptheta, arma::vec & Mt)
    /* ---------------------------------------------------------------------
        - Generate 2-D diffusion process PDF tables as well as Ptheta & Mt.
        - Calculate first-passage-time density and response probabilities for
          circular diffusion process; badix usually is 0.
       --------------------------------------------------------------------- */
    {
        double tscale, G11, G21, G12, G22, Girs1, Girs2;
        unsigned int i, j;
        double totalmass = 0;

        /* 1. Density of zero-drift process */
        d0_GSL(DT, Gt0);

        /* 2. Generate a series of angles; Response circle (1 x nw) */
        R[0] = -M_PI;
        for (i = 1; i < m_nw; i++) R[i] = R[i-1] + m_w;

        /* 3. Joint RT distribution/PDF (nw * sz) */
        for (i = 0; i < m_sz; i++)
        {
          // When eta1 = eta2 = 0, tscale becomes 1
          // Otherwise it depends on DT.
            tscale = std::sqrt(1/(1 + m_eta1onsigma2 * DT[i])) *
                     std::sqrt(1/(1 + m_eta2onsigma2 * DT[i]));
            // Rprintf("[m_eta1onsigma2, m_eta2onsigma2, DT[i]] [%6.4f, %6.4f, %6.4f]\n",
            //         m_eta1onsigma2, m_eta2onsigma2, DT[i]);

            G11 = 2.0*m_eta1*m_eta1 * (1 + m_eta1onsigma2 * DT[i]);
            G21 = 2.0*m_eta2*m_eta2 * (1 + m_eta2onsigma2 * DT[i]);
            // Rprintf("[tscale, G11, G21] [%6.4f, %6.4f, %6.4f]\n",tscale, G11, G21);

            for (j = 0; j < m_nw; j++)
            {
              G12 = m_v1 + m_a * m_eta1onsigma2 * std::cos(R[j]);
              G22 = m_v2 + m_a * m_eta2onsigma2 * std::sin(R[j]);
              // Rprintf("[G12, G22] [%6.4f, %6.4f]\n", G12, G22);

              Girs1 = std::exp( (G12*G12) / G11 - (m_v1*m_v1) / (m_eta1*m_eta1) / 2.);
              Girs2 = std::exp( (G22*G22) / G21 - (m_v2*m_v2) / (m_eta2*m_eta2) / 2.);

              // Rprintf("[Girs1, Girs2] [%6.4f, %6.4f]\n", Girs1, Girs2);

              // Identical
              // Gt[m_nw * i + j] = tscale*Girs1*Girs2*Gt0[i] / M_2PI;
              Gt(j, i) = tscale*Girs1*Girs2*Gt0[i] / M_2PI;
            }
        }

        /* 4. Total mass */
        for (j = 0; j < m_nw; j++)    // start from row
        {
          for (i = 1; i < m_sz; i++)  // Note start from 1
          {
            totalmass += (Gt[m_nw * i + j] + Gt[m_nw * (i - 1) + j]) / 2.0;
          }
        }
        totalmass *= m_w * m_h;
        // Rprintf("totalmass = %6.4f\n", totalmass); // .9990744

        //  Ptheta = a distribution of report outcomes,
        //  Mt     = distribution of (constant!) means.
        /* 5. Integrate joint densities to get means hitting probabilities */
        for (j = 0; j < m_nw; j++)
        {
          Ptheta[j] = 0;  // Each discrete theta integrates separately
          Mt[j] = 0;

          for (i = 1; i < m_sz; i++)
          {
              Ptheta[j] += (Gt[m_nw*i + j] + Gt[m_nw*(i - 1) + j]) /2.0;
              Mt[j] += (DT[i] * Gt[m_nw*i + j] + DT[i - 1] * Gt[m_nw*(i - 1)+ j]) / 2.0;
          }
          Ptheta[j] *= m_h / totalmass;
          Mt[j]     *= m_h / Ptheta[j] / totalmass;
        }
    }

    void r(unsigned int n, arma::mat & out)
    {
      arma::rowvec sigma_wt1, sigma_wt2;
      arma::mat Sigma_Wt, Xt(2, m_sz);
      double theta, Dt2;      // Distance_t^2
      unsigned int j;  // time step counter
      arma::rowvec mut1 = m_v1*arma::ones<arma::rowvec>(m_sz);
      arma::rowvec mut2 = m_v2*arma::ones<arma::rowvec>(m_sz);
      arma::mat Mut = arma::join_vert(mut1, mut2); // 2 x m_sz

      for(size_t i = 0; i<n; i++)
      {
        j   = 1;    // start from 0, move to the 1st, so 1
        Dt2 = 0;    // starting rPos = 0
        Xt.zeros(); // re-position back to the origin

        sigma_wt1 = m_sigma1*arma::randn<arma::rowvec>(m_sz);
        sigma_wt2 = m_sigma2*arma::randn<arma::rowvec>(m_sz);
        Sigma_Wt  = arma::join_vert(sigma_wt1, sigma_wt2); // 2 x m_sz

        while (Dt2 < m_a2 && j < m_sz)
        {
          Xt.col(j) = Xt.col(j-1) + Mut.col(j)*m_h + Sigma_Wt.col(j)*std::sqrt(m_h);
          Dt2 = Xt(0, j)*Xt(0, j) + Xt(1, j)*Xt(1, j);
          j++;
        }
        theta = std::atan2(Xt(1, j-1), Xt(0, j-1));
        if ( m_nw != 0 ) { out(i, 0) = divider(theta); } // R
        out(i, 1) = j*m_h + m_t0;   // RT
        out(i, 2) = theta;          // A
      }
    }

    bool ValidateParams (bool print)
    {
      using namespace Rcpp;
      bool valid = true;

      if (m_a < 0)      { valid = false; if (print) Rcout << "invalid parameter a = " << m_a << std::endl; }
      if (m_t0 < 0)     { valid = false; if (print) Rcout << "invalid parameter t0 = " << m_t0 << std::endl; }
      if (m_sigma1 < 0) { valid = false; if (print) Rcout << "invalid parameter sigma1 = " << m_sigma1 << std::endl; }
      if (m_sigma2 < 0) { valid = false; if (print) Rcout << "invalid parameter sigma2 = " << m_sigma2 << std::endl; }

      return valid;
    }

    double divider(double theta)
    {
      double out = NA_REAL;

      for(size_t i=0; i<(m_nw+1); i++)
      {
        if ( theta >= (m_thetai[i] - m_w/2.) && theta < (m_thetai[i] + m_w/2.) )
        {
          out = (i == m_nw) ? m_thetai[0] : m_thetai[i];
          break;
        }
      }
      return out;
    }

    /* The following 4 functions are from GSL by by G. Jungman. */
    /* Extracted from GNU GSL specfunc/bessel_J1.c. They remove error checking.
     * Use with caution.
     */
    double gsl_sf_bessel_zero_J0(unsigned int s)
    {
      /* CHECK_POINTER(result) */
      double out;

      if(s == 0)
      {
        out = 0.0;
        Rcpp::stop("bessel_zero_J0 error");
      }
      else
      {
        /* See [F. Lether, J. Comp. Appl .Math. 67, 167 (1996)]. */
        static const double P[] = { 1567450796.0/12539606369.0,
                                    8903660.0/2365861.0,
                                    10747040.0/536751.0,
                                    17590991.0/1696654.0 };
        static const double Q[] = { 1.0,
                                    29354255.0/954518.0,
                                    76900001.0/431847.0,
                                    67237052.0/442411.0 };

        const double beta = (s - 0.25) * M_PI;
        const double bi2  = 1.0/(beta*beta);
        const double R33num = P[0] + bi2 * (P[1] + bi2 * (P[2] + P[3] * bi2));
        const double R33den = Q[0] + bi2 * (Q[1] + bi2 * (Q[2] + Q[3] * bi2));
        const double R33 = R33num/R33den;
        out = beta + R33/beta;
      }
      return out;
    }

    /* From GNU GSL specfunc/cheb_eval.c */
    double cheb_eval_e(const cheb_series * cs, const double x)
    {
        int j;
        double d  = 0.0;
        double dd = 0.0;

        double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
        double y2 = 2.0 * y;

        double e = 0.0;

        for(j = cs->order; j>=1; j--) {
          double temp = d;
          d = y2*d - dd + cs->c[j];
          e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
          dd = temp;
        }

        {
          double temp = d;
          d = y*d - dd + 0.5 * cs->c[0];
          e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
        }

        return d;
    }

    /* From GNU GSL specfunc/bessel.c */
    double gsl_sf_bessel_sin_pi4_e(double y, double eps)
    {
      double out;

      const double sy = sin(y);
      const double cy = cos(y);
      const double s = sy + cy;
      const double d = sy - cy;
      // const double abs_sum = fabs(cy) + fabs(sy);
      double seps;
      double ceps;
      if(fabs(eps) < GSL_ROOT5_DBL_EPSILON) {
        const double e2 = eps*eps;
        seps = eps * (1.0 - e2/6.0 * (1.0 - e2/20.0));
        ceps = 1.0 - e2/2.0 * (1.0 - e2/12.0);
      }
      else {
        seps = sin(eps);
        ceps = cos(eps);
      }
      out = (ceps * d + seps * s)/ M_SQRT2;
      // result->err = 2.0 * GSL_DBL_EPSILON * (fabs(ceps) + fabs(seps)) * abs_sum / M_SQRT2;

      /* The following comment is from gsl-2.6
       * Try to account for error in evaluation of sin(y), cos(y).
       * See above.
       * FIXME ?
       */
      if(y > 1.0/GSL_DBL_EPSILON) {
        Rcpp::Rcout << "GSL error\n";
      }
      else if(y > 1.0/GSL_SQRT_DBL_EPSILON) {
        Rcpp::Rcout << "GSL error\n";
      }

      return out;
    }

    /* From GNU GSL specfunc/bessel_J1.c */
    double gsl_sf_bessel_J1(const double x)
    {
      double y = fabs(x), out;

      /* CHECK_POINTER(result) */

      if(y == 0.0) {
        out = 0.0;
      }
      else if(y < 2.0*GSL_DBL_MIN) {
        Rcpp::stop("UNDERFLOW_ERROR");
      }
      else if(y < ROOT_EIGHT * GSL_SQRT_DBL_EPSILON) {
        out = 0.5*x;
        return out;
      }
      else if(y < 4.0) {
        double c = cheb_eval_e(&bj1_cs, 0.125*y*y-1.0);
        out = x * (0.25 + c);
      }
      else {
        /* Because the leading term in the phase is y,
         * which we assume is exactly known, the error
         * in the cos() evaluation is bounded.
         */
        const double z  = 32.0/(y*y) - 1.0;
        double ca = cheb_eval_e(&_gsl_sf_bessel_amp_phase_bm1_cs,  z);
        double ct = cheb_eval_e(&_gsl_sf_bessel_amp_phase_bth1_cs, z);
        double sp = gsl_sf_bessel_sin_pi4_e(y, ct/y);

        const double sqrty = sqrt(y);
        const double ampl  = (0.75 + ca) / sqrty;
        out  = (x < 0.0 ? -ampl : ampl) * sp;
      }
      return out;
    }

};

arma::vec interp2(arma::vec DT_, arma::vec R_, arma::vec DT, arma::vec R,
                  arma::mat Gt, unsigned int sz, unsigned int nw);

arma::vec dcircle(arma::vec RT, arma::vec A, arma::vec P, double tmax,
                  unsigned int kmax, unsigned int sz, unsigned int nw);
#endif
