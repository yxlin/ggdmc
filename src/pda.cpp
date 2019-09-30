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

//' @export
// [[Rcpp::export]]
arma::vec spdf(arma::vec x, arma::vec RT, int n, double h_in,
               bool debug) {
  unsigned int nx  = x.n_elem;
  unsigned int nRT = RT.n_elem; // if defective densities, nRT != n
  double h, z0, z1, z1minusz0, dt, fil0_constant, minRT, maxRT;
  unsigned int ngrid=1024, half_ngrid=512;
  arma::vec out(nx);

  minRT = RT.min();
  maxRT = RT.max();
  h  = (h_in == 0) ? (0.8*arma::stddev(RT)*std::pow(nRT, -0.2)) : h_in;
  z0 = (minRT - 3.0*h < 0) ? 1e-10 : (minRT - 3.0*h);
  if (z0 < 0) { z0 = 0; if (debug) Rcout <<"z0 in SPDF is less than 0\n"; }
  z1 = maxRT > 10.0 ? 10.0 : maxRT + 3.0*h;

  if (nRT <= 10) {
    out.fill(1e-10);
  } else if (z1 <= z0) {

    if (debug)
    {
      Rcpp::Rcout << "[minRT maxRT nsRT] " << minRT <<  " " << maxRT << " " <<
        nRT << "\n";
    }
    out.fill(1e-10);

  } else {



    arma::vec z = arma::linspace<arma::vec>(z0, z1, ngrid);
    dt = z[1] - z[0];

    z1minusz0 = z1 - z0;
    fil0_constant = (-2.0*h*h*M_PI*M_PI) / (z1minusz0*z1minusz0);

    arma::vec filter0(ngrid);
    arma::vec h_binedge0(ngrid + 1);
    arma::vec signal0(ngrid);

    // Get binedge (1025), filter (1024) and histogram (1024) at one go -----------------
    for(size_t i=0; i < ngrid; i++)
    {
      h_binedge0[i] = z0 + dt*((double)i - 0.5); // Binedge

      if (i < (1 + half_ngrid)) {                // Filter
        filter0[i] = std::exp(fil0_constant * (double)(i*i));
      } else {
        int j = 2*(i - half_ngrid); // flipping
        filter0[i] = filter0[i-j];
      }
    }

    h_binedge0[ngrid] = (z0 + ((double)(ngrid - 1))*dt);

    arma::vec h_hist0 = arma::conv_to<arma::vec>::from(arma::histc(RT, h_binedge0)); // 1025
    signal0 = h_hist0.rows(0, ngrid-1) / (dt * (double)(n));

    arma::vec sPDF = arma::real(arma::ifft(filter0 % arma::fft(signal0))) ;
    arma::vec eDen; // a container for estiamted densities

    if (z.has_nan()) {
      eDen.fill(1e-10);
      Rcpp::stop("z has nan");
    } else {
      arma::vec uniquez = arma::unique(z);
      if (uniquez.size() <= 1) {
        eDen.fill(1e-10);
        Rcpp::stop("z has only 1 or 0 element");
      } else {
        arma::interp1(z, sPDF, x, eDen);
      }
    }

    for(size_t i=0; i < nx; i++)
    {
      out[i] = (eDen[i] < 1e-10 || std::isnan(eDen[i])) ? 1e-10 : eDen[i];
    }
  }

  return out;
}

