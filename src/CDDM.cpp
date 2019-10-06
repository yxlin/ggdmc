//    Copyright (C) <2019>  <Yi-Shin Lin>
//
//    This program is free software; you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation; either version 2 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License along
//    with this program; if not, write to the Free Software Foundation, Inc.,
//    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#include <ggdmc.hpp>

using namespace Rcpp;

//' The 2-dimension Diffusion Model
//'
//' Density, random generation for the 2-D diffusion model.
//'
//' The model has the main parameters, v1, v2, eta1, eta2, a, sigma, and t0.
//' tmax, kmax, sz and nw are tuning parameters for determining the set.
//' dcircle300 produces PDF table and others.
//'
//' @param x, n x is a data.matrix with first column of RT and a second column
//' of R. n is the numbers of observation.
//' @param P is a parameter vector, c(v1, v2, a, t0, sigma1, sigma2, eta1, eta2).
//' The sequence is important. v1 is the x-axis mean drift rate. v2 is the
//' y-axis mean drift rate. sigma1 is the x-axis within-trial drift rate SD.
//' sigma2 is the y-axix within-trial drift rate SD. a is decision threshold.
//' sigma1 and sigma2 must be 1 and identical, because this is what has been
//' thoroughly tested so far. Other values may return unknown results. t0
//' non-decision time.
//' @param tmax maximum time of the model
//' @param kmax the tuning parameter for Bessel function. Mostly 50.
//' @param h, sz the number of time steps (h = tmax / sz). h is time step.
//' Mostly .1 ms.
//' @param nw the number of theta steps (w = 2 * pi / nw)
//'
//' @return rcircle returns a n x 2 matrix. Each row is an [RT R] trial.
//' dcircle returns a n vector.
//' @examples
//' ## TODO examples
//' @export
// [[Rcpp::export]]
arma::vec dcircle(arma::vec RT, arma::vec A, arma::vec P, double tmax,
                  unsigned int kmax, unsigned int sz, unsigned int nw)
{
  // P[0] = v1; P[1] = v2; P[2] = a; P[3] = t0; P[4] = sigma1; P[5] = sigma2;
  // P[6] = eta1; P[7] = eta2;
  if (P[4] != P[5]) Rcpp::stop("Drift rate SD must be the same in both directions.");
  cddm * obj = new cddm(P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], tmax,
                        kmax, nw, sz);
  bool isvalid = obj->ValidateParams(false); // Not to print invalid parameter guesses

  arma::vec DT_tab(sz), Gt0(sz), A_tab(nw), Ptheta(nw), Mt(nw), DT;
  arma::mat Gt(nw,sz);
  Gt0.zeros(sz);
  arma::vec out(RT.n_elem);

  if (!isvalid)
  {
    out.fill(1e-10);
  }
  else
  {
    obj->d(DT_tab, A_tab, Gt0, Gt, Ptheta, Mt);
    DT = RT - P[3];
    out = interp2(DT, A, DT_tab, A_tab, Gt, sz, nw);
  }

  delete obj;
  return out;
}

//' @rdname dcircle
//' @export
// [[Rcpp::export]]
Rcpp::List dcircle300(arma::vec P, double tmax, unsigned int kmax,
                      unsigned int sz, unsigned int nw)
{
  if (P[4] != P[5]) Rcpp::stop("Drift rate SD must be the same in both directions.");
  cddm * obj = new cddm(P[0], P[1], P[2], P[3], P[4], P[5], P[6], P[7], tmax,
                        kmax, nw, sz);
  arma::vec DT(sz), Gt0(sz), R(nw), Ptheta(nw), Mt(nw);
  arma::mat Gt(nw, sz);
  Gt0.zeros(sz);
  obj->d(DT, R, Gt0, Gt, Ptheta, Mt);
  delete obj;

  return Rcpp::List::create(Rcpp::Named("DT")  = DT,
                            Rcpp::Named("R")   = R,
                            Rcpp::Named("Gt0") = Gt0,
                            Rcpp::Named("Gt")  = Gt,
                            Rcpp::Named("Ptheta") = Ptheta,
                            Rcpp::Named("Mt")     = Mt);
}

//' @rdname dcircle
//' @export
// [[Rcpp::export]]
arma::mat rcircle(unsigned int n, arma::vec P, double tmax, double h,
                  unsigned int nw)
{
  // P[0] = v1; P[1] = v2; P[2] = a; P[3] = t0; P[4] = sigma1; P[5] = sigma2;
  // P[6] = eta1; P[7] = eta2;
  cddm * obj = new cddm(P[0], P[1], P[2], P[3], P[4], P[5], tmax, h, nw);
  arma::mat out(n, 3); out.fill(NA_REAL);
  obj->r(n, out);
  delete obj;
  return out;   // [R RT A]
}

//' @rdname dcircle
//' @export
// [[Rcpp::export]]
Rcpp::List rcircle_process(arma::vec P, double tmax, double h)
{
  // P[0] = v1; P[1] = v2; P[2] = a; P[3] = t0; P[4] = sigma1; P[5] = sigma2;
  // P[6] = eta1; P[7] = eta2;
  if ( h <= 0  ) Rcpp::stop("h must be greater than 0.");
  if (tmax <= 0) Rcpp::stop("tmax must be greater than 0.");
  if (tmax < 1 ) Rcpp::Rcout << "tmax less than 1.\n";

  arma::vec T, sigma_wt1, sigma_wt2, mut1, mut2, out(3);
  arma::mat Sigma_Wt, Mut, tmp0, tmp1;
  double a2=P[2]*P[2], Dt2=0; // Distance_t^2

  T = arma::regspace(0, h, tmax);   // h must > 0
  unsigned int nmax=T.n_elem, i=1;

  mut1      = P[0]*arma::ones<arma::vec>(nmax);
  mut2      = P[1]*arma::ones<arma::vec>(nmax);
  sigma_wt1 = P[4]*arma::randn(nmax); // column vector
  sigma_wt2 = P[5]*arma::randn(nmax);
  tmp0      = arma::join_horiz(sigma_wt1, sigma_wt2); // nmax x 2
  tmp1      = arma::join_horiz(mut1, mut2); // nmax x 2
  Sigma_Wt  = tmp0.t(); // 2 x nmax
  Mut       = tmp1.t(); // 2 x nmax

  arma::vec xPos(nmax); xPos.fill(NA_REAL);
  arma::vec yPos(nmax); yPos.fill(NA_REAL);
  xPos(0) = 0;
  yPos(0) = 0;
  // arma::mat rPos(maxIter); rPos.fill(NA_REAL);
  arma::mat Xt(2, nmax); Xt.zeros();

  // change from i <= nmax to i < nmax, bc the last element in C++ is (nmax-1)
  while (Dt2 < a2 && i < nmax)
  {
    Xt.col(i) = Xt.col(i-1) + Mut.col(i)*h + Sigma_Wt.col(i)*std::sqrt(h);
    Dt2 = Xt(0, i)*Xt(0, i) + Xt(1, i)*Xt(1, i); // 2nd element; rPos
    xPos(i) = Xt(0, i-1);
    yPos(i) = Xt(1, i-1);
    i++;
  }
  // No need to minus 1, bc the 2nd element in C++ index is 1 not 2.
  out(0) = i * h + P[3];
  // Still need minus 1, bc the last time step is still (i-1)'th index.
  out(1) = std::atan2(Xt(1, i-1), Xt(0, i-1));

  // if this trial is a non-terminated trial. That is, if this trial has not
  // hit the boundary after tmax.
  out(2) = i > nmax; // auto-upcast
  return Rcpp::List::create(Rcpp::Named("out")  = out,
                            Rcpp::Named("xPos") = xPos,
                            Rcpp::Named("yPos") = yPos);

}

//' @rdname dcircle
//' @export
// [[Rcpp::export]]
Rcpp::List r1d(arma::vec P, double tmax, double h)
{
  // P[0] = v; P[1] = a; P[2] = z; P[3] = t0; P[4] = s;
  if ( h <= 0  ) Rcpp::stop("h must be greater than 0.");
  if (tmax <= 0) Rcpp::stop("tmax must be greater than 0.");
  if (tmax < 1 ) Rcpp::Rcout << "tmax less than 1.\n";
  if (P[2] > P[1]) Rcpp::stop("z > a");
  if (P[3] < 0) Rcpp::stop("t0 > 0");

  arma::vec T, sigma_wt, mut, out(2);
  double Dt;

  T = arma::regspace(0, h, tmax);   // h must > 0
  unsigned int nmax=T.n_elem, i=1;

  mut      = P[0]*arma::ones<arma::vec>(nmax);
  sigma_wt = P[4]*arma::randn(nmax);
  arma::vec  Xt(nmax); Xt.fill(NA_REAL);
  Xt(0) = P[2];
  Dt    = Xt(0);

  while (Dt < P[1] && Dt > 0 && i < nmax)
  {
    Xt(i) = Xt(i-1) + mut(i)*h + sigma_wt(i)*std::sqrt(h);
    Dt = Xt(i);
    i++;
  }
  out(0) = i * h + P[3];
  out(1) = i > nmax; // auto-upcast

  return Rcpp::List::create(Rcpp::Named("T")  = T,
                            Rcpp::Named("out")= out,
                            Rcpp::Named("Xt") = Xt);

}

arma::vec interp2(arma::vec DT_, arma::vec R_, arma::vec DT, arma::vec R,
                  arma::mat Gt, unsigned int sz, unsigned int nw)
{
  // ------------------------------------------------------------------
  // A workable but untidy version of circular 2-D linear interpolation
  // This function does not require monotonically increasing values.
  // Armadillo's interp2 won't work, bc it mandates monotonically
  // increasing XI and YI, which is not possible considering the real data of
  // c(RT, R) pairs will not increase/decrease monotonically, even we order
  // them.
  // - DT_ and R_ are from empirical data;
  // - DT and R are from PDF table.
  // ------------------------------------------------------------------
  unsigned int i, j, k;
  unsigned int n = DT_.n_elem;
  double Delta_x, Delta_y, delta_x, delta_y, fxy0, fxy1;
  arma::vec out(n), xcoor(2), ycoor(2);
  out.fill(1e-10);

  arma::vec tmp0(1); tmp0.fill(M_PI);
  arma::mat Gt_circular = arma::join_vert(Gt, Gt.row(0));
  arma::vec R_circular  = arma::join_vert(R, tmp0);

  for (i=0; i<n; i++)
  {
    xcoor.fill(NA_REAL); // Safeguard for out-of-range data
    ycoor.fill(NA_REAL);

    for(j=0; j<sz; j++)
    {
      if ( DT_[i] > DT[sz-1] )
      {
        // Rprintf("DT %6.4f greater than last DT %6.4f on the table\n", DT_[i], DT[sz-1]);
        break;
      }
      if ( (DT_[i] >= DT[j]) && (DT_[i] <  DT[j+1]))
      {
        xcoor(0) = j;   // lower point 0
        xcoor(1) = j+1; // higher point 1
        break;
      }
    }

    // for(k=0; k<(nw); k++)  // circular sheet
    for(k=0; k<(nw+1); k++)  // circular sheet
    {
      // if ( (R_[i] >= R[k]) && (R_[i] <  R[k+1]) )
      if ( (R_[i] >= R[k]) && (R_[i] <  R_circular[k+1]) )
      {
        ycoor(0) = k;   // lower point 0
        ycoor(1) = k+1; // higher point 1
        break;
      }

      // if (R_[i] == M_PI)
      // {
      //   // Rprintf("R_ %6.4f equals PI\n", R_[i]);
      //   ycoor(0) = 0;  // temporary fix for the unlikely exact match of the
      //   ycoor(1) = 1;  // end point
      //   break;
      // }

    }
    //   // Gt(k, j);       f(0,0)
    //   // Gt(k, j+1);     f(0,1)
    //   // Gt(k+1, j);     f(1,0)
    //   // Gt(k+1, j+1);   f(1,1)

    if ( xcoor.has_nan() || ycoor.has_nan() )
    {
      // xcoor.print("NA_REAL found. xcoor");
      // ycoor.print("NA_REAL found. ycoor");
      out[i] = 1e-10;
    } else
    {
      Delta_x = DT[xcoor(1)] - DT[xcoor(0)];
      delta_x = DT_[i] - DT[xcoor(0)];

      Delta_y = R_circular[ycoor(1)] - R_circular[ycoor(0)];
      delta_y = R_[i] - R_circular[ycoor(0)];
      // Delta_y = R[ycoor(1)] - R[ycoor(0)];
      // delta_y = R_[i] - R[ycoor(0)];

      // 2-D Interpolation equation
      fxy0 = (delta_x / Delta_x) * ( Gt_circular( ycoor(0), xcoor(1)) -
        Gt_circular(ycoor(0), xcoor(0)) ) + Gt_circular(ycoor(0), xcoor(0));
      fxy1 = (delta_x / Delta_x) * ( Gt_circular(ycoor(1), xcoor(1)) -
        Gt_circular(ycoor(0), xcoor(1)) ) + Gt_circular(ycoor(0), xcoor(1));
      // fxy0 = (delta_x / Delta_x) * ( Gt( ycoor(0), xcoor(1)) -
      //   Gt(ycoor(0), xcoor(0)) ) + Gt(ycoor(0), xcoor(0));
      // fxy1 = (delta_x / Delta_x) * ( Gt(ycoor(1), xcoor(1)) -
      //   Gt(ycoor(0), xcoor(1)) ) + Gt(ycoor(0), xcoor(1));

      out[i] = (delta_y/Delta_y) * (fxy1-fxy0) + fxy0;
    }
  }

  return out;
}

double test_gsl(double s, arma::vec P, double tmax, double h,
                unsigned int nw)
{
  cddm * obj = new cddm(P[0], P[1], P[2], P[3], P[4], P[5], tmax, h, nw);
  // double out = obj->gsl_sf_bessel_zero_J0(s);
  double out = obj->gsl_sf_bessel_J1(s);
  delete obj;
  return out;
}
