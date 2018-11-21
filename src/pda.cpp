#include <ggdmc.hpp>

using namespace Rcpp;

// [[Rcpp::plugins(openmp)]]
// Protect against compilers without OpenMP; eg OS X  clang
#ifdef _OPENMP
#include <omp.h>
#endif

#define OMP_MIN_VALUE 1e6

arma::mat check_finite(arma::vec RT, arma::vec R, bool debug) {

  arma::mat out;

  if (RT.has_inf()) {
    arma::vec finite_RT, finite_R;
    arma::uvec finite_idx = arma::find_finite(RT);
    finite_RT = RT.elem(finite_idx);
    finite_R  = R.elem(finite_idx);
    if (debug) Rcout << "Inf/NaN found. Only " << finite_idx.n_elem << " RTs are returned\n";
    out = arma::join_horiz(finite_RT, finite_R);
  } else {
    out = arma::join_horiz(RT, R);
  }

  return out;
}


//' Simulated likelihood functions for RT models
//'
//' CPU-based PDA
//'
//' @param x empirical data
//' @param RT simulated resposne times
//' @param n number of model simulations
//' @param h_in kernel bandwidth
//' @param debug debugging?
//' @return a vector
//' @export
// [[Rcpp::export]]
arma::vec spdf(arma::vec x, arma::vec RT, unsigned int n, double h_in,
  bool debug = false) {
  // do not change "int n" to "unsigned int n"
  unsigned int nx  = x.n_elem;
  unsigned int nRT = RT.n_elem; // if defective densities, nRT != n
  arma::vec out(nx);
  if (nRT <= 10) {
    out.fill(1e-10);
  } else {
    double minRT = RT.min();
    double maxRT = RT.max();

    double h  = (h_in == 0) ? (0.8*arma::stddev(RT)*std::pow(nRT, -0.2)) : h_in;
    double z0 = minRT <= 1e-10 ? 1e-10 : minRT - 3.0*h;
    if (z0 < 0) { z0 = 0; if (debug) Rcout <<"z0 in SPDF is less than 0\n"; }
    double z1 = maxRT > 10.0 ? 10.0 : maxRT + 3.0*h;
    unsigned int ngrid = 1024;
    unsigned int half_ngrid  = 0.5*ngrid;
    arma::vec z = arma::linspace<arma::vec>(z0, z1, ngrid);
    double dt = z[1] - z[0];
    double z1minusz0 = z1 - z0;
    double fil0_constant = (-2.0*h*h*M_PI*M_PI) / (z1minusz0*z1minusz0);

    arma::vec filter0(ngrid);
    arma::vec h_binedge0(ngrid + 1);
    arma::vec signal0(ngrid);

    // Get binedge (1025), filter (1024) and histogram (1024) at one go --------
    for(size_t i=0; i < ngrid; i++) {
      h_binedge0[i] = z0 + dt*((double)i - 0.5); // Binedge
      if (i < (1 + half_ngrid)) {                // Filter
        filter0[i] = std::exp(fil0_constant * (double)(i*i));
      } else {
        int j = 2*(i - half_ngrid); // flipping
        filter0[i] = filter0[i-j];
      }
    }

    h_binedge0[ngrid] = (z0 + ((double)(ngrid - 1))*dt); // 1025
    arma::vec h_hist0 = arma::conv_to<arma::vec>::from(arma::histc(RT, h_binedge0));
    signal0 = h_hist0.rows(0, ngrid-1) / (dt * (double)(n));

    arma::vec sPDF = arma::real(arma::ifft(filter0 % arma::fft(signal0))) ;
    arma::vec eDen; // a container for estiamted densities

    if(z.has_nan()) {
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

    for(size_t i = 0; i < nx; i++) {
      out[i] = (eDen[i] < 1e-10 || std::isnan(eDen[i])) ? 1e-10 : eDen[i];
    }
  }
  return out;
}

//' @rdname rplba1R
//' @export
// [[Rcpp::export]]
arma::mat rplba0(unsigned int n, arma::vec A, arma::vec b, double t0, arma::vec mean_v,
  arma::vec mean_w, arma::vec sd_v, double rD, double swt, int ncore = 1,
  bool debug = false) {

  int nmean_v = mean_v.n_rows;
  if (debug & (nmean_v < 2)) Rcpp::stop("Minimal 2 accumulators/responses!");
  if (debug & arma::any(sd_v < 0)) Rcout << "sd_v cannot be negative.\n";

  if (sd_v.n_elem == 1) sd_v = arma::repmat(sd_v, nmean_v, 1);
  if (A.n_elem == 1)  A = arma::repmat(A, nmean_v, 1);
  if (b.n_elem == 1)  b = arma::repmat(b, nmean_v, 1);

  double T0  = rD + swt;
  double dt0_stage1, dt0_stage2, dt1_stage1, dt1_stage2;
  arma::vec RT(n), R(n);
  size_t i;

    for (i = 0; i < n; i++) {
      double v0, v1, w0, w1, x0, x1, z0, z1, DT_tmp;
      int R_tmp;

      v0 = rtn_scalar(mean_v(0), sd_v(0), 0, INFINITY) ; // acc1 (X)
      v1 = rtn_scalar(mean_v(1), sd_v(1), 0, INFINITY) ; // acc2 (O)
      w0 = rtn_scalar(mean_v(0), sd_v(0), 0, INFINITY) ; // acc1 (X)
      w1 = rtn_scalar(mean_v(1), sd_v(1), 0, INFINITY) ; // acc2 (O)
      x0 = A(0) * R::runif(0, 1); // Stage 1 starting pos choice 0
      x1 = A(1) * R::runif(0, 1); // Stage 1 starting pos choice 1
      z0 = x0 + T0*v0; // Stage 2 starting pos choice 0
      z1 = x1 + T0*v1; // Stage 2 starting pos choice 1
      dt0_stage1 = (b(0) - x0) / v0;
      dt1_stage1 = (b(1) - x1) / v1;
      DT_tmp = dt0_stage1 < dt1_stage1 ? dt0_stage1 : dt1_stage1;
      R_tmp  = dt0_stage1 < dt1_stage1 ? 1 : 2;

      if (DT_tmp <= T0) {
        RT(i) = DT_tmp + t0;
        R(i)  = R_tmp;
      } else {
        dt0_stage2 = (b(0) - z0) / w0;
        dt1_stage2 = (b(1) - z1) / w1;
        RT(i) = dt0_stage2 < dt1_stage2 ? dt0_stage2 + T0 + t0 : dt1_stage2 + T0 + t0;
        R(i)  = dt0_stage2 < dt1_stage2 ? 1 : 2;
      }
    }

    arma::mat out = check_finite(RT, R, debug);
  return out;
}

//' @rdname rplba1R
//' @export
// [[Rcpp::export]]
arma::mat rplba1(unsigned int n, arma::vec A, arma::vec b, double t0,
  arma::vec mean_v, arma::vec mean_w, arma::vec sd_v, double rD, double swt,
  unsigned int ncore = 1, bool debug = false) {

  unsigned int nmean_v = mean_v.n_rows;
  if (debug & (nmean_v < 2)) Rcpp::stop("Minimal 2 accumulators/responses!");
  if (debug & arma::any(sd_v < 0)) Rcout << "sd_v cannot be negative.\n";

  if (sd_v.n_elem == 1) sd_v = arma::repmat(sd_v, nmean_v, 1);
  if (A.n_elem == 1)  A = arma::repmat(A, nmean_v, 1);
  if (b.n_elem == 1)  b = arma::repmat(b, nmean_v, 1);

  double T0  = rD + swt;
  double dt0_stage1, dt0_stage2, dt1_stage1, dt1_stage2;
  arma::vec RT(n), R(n);
  size_t i;

  for (i = 0; i < n; i++) {
    double v0, v1, w0, w1, x0, x1, z0, z1, DT_tmp;
    int R_tmp;

    v0 = rtn_scalar(mean_v(0), sd_v(0), 0, INFINITY) ; // acc1 (X)
    v1 = rtn_scalar(mean_v(1), sd_v(1), 0, INFINITY) ; // acc2 (O)
    w0 = rtn_scalar(mean_w(0), sd_v(0), 0, INFINITY) ; // acc1 (X)
    w1 = rtn_scalar(mean_w(1), sd_v(1), 0, INFINITY) ; // acc2 (O)
    x0 = A(0) * R::runif(0, 1); // Stage 1 starting pos choice 0
    x1 = A(1) * R::runif(0, 1); // Stage 1 starting pos choice 1
    z0 = x0 + T0*v0; // Stage 2 starting pos choice 0
    z1 = x1 + T0*v1; // Stage 2 starting pos choice 1
    dt0_stage1 = (b(0) - x0) / v0;
    dt1_stage1 = (b(1) - x1) / v1;
    DT_tmp = dt0_stage1 < dt1_stage1 ? dt0_stage1 : dt1_stage1;
    R_tmp  = dt0_stage1 < dt1_stage1 ? 1 : 2;

    if (DT_tmp <= T0) {
      RT(i) = DT_tmp + t0;
      R(i)  = R_tmp;
    } else {
      dt0_stage2 = (b(0) - z0) / w0;
      dt1_stage2 = (b(1) - z1) / w1;
      RT(i) = dt0_stage2 < dt1_stage2 ? dt0_stage2 + T0 + t0 : dt1_stage2 + T0 + t0;
      R(i)  = dt0_stage2 < dt1_stage2 ? 1 : 2;
    }
  }

  arma::mat out = check_finite(RT, R, debug);
  return out;
}


//' @rdname rplba1R
//' @export
// [[Rcpp::export]]
arma::mat rplba2(unsigned int n, arma::vec A, arma::vec b, double t0,
  arma::vec mean_v, arma::vec mean_w, arma::vec sd_v, arma::vec sd_w, double rD,
  double swt, unsigned int ncore = 1, bool debug = false) {
  // This is non-parallel version
  arma::vec x0(n), x1(n), v0(n), v1(n), dt0(n), dt1(n);
  arma::vec R(n), winDT(n), undone(n);
  double eswt = rD + swt;
  size_t i;

  for(i = 0; i < n; i++)
  {  // Stage 1 LBA
    v0[i] = rtn_scalar(mean_v[0], sd_v[0], 0, INFINITY) ; // acc1 (X)
    v1[i] = rtn_scalar(mean_v[1], sd_v[1], 0, INFINITY) ; // acc2 (O)
    x0[i] = R::runif(0, A[0]);
    x1[i] = R::runif(0, A[1]);
    dt0[i] = (b[0] - x0[i]) / v0[i];
    dt1[i] = (b[1] - x1[i]) / v1[i];
    R[i] = (dt0[i] < dt1[i]) ? 1 : 2;
    winDT[i]  = (dt0[i] < dt1[i]) ? dt0[i] : dt1[i];
    undone[i] = (winDT[i] <= eswt) ? false : true;
  }

  unsigned int n2 = arma::accu(undone);
  arma::uvec idx = find(undone);
  arma::vec z0 = x0.elem(idx);
  arma::vec z1 = x1.elem(idx);
  arma::vec w0 = v0.elem(idx);
  arma::vec w1 = v1.elem(idx);
  arma::vec dt0s2(n2), dt1s2(n2);

  // Stage 2 LBA
  for (size_t j = 0; j < n2; j++)
  {
    dt0s2[j] = (b[0] - (z0[j] + eswt * w0[j])) / (rtn_scalar(mean_w[0], sd_w[0], 0, INFINITY));
    dt1s2[j] = (b[1] - (z1[j] + eswt * w1[j])) / (rtn_scalar(mean_w[1], sd_w[1], 0, INFINITY));
    R[idx[j]]      = (dt0s2[j] < dt1s2[j]) ? 1 : 2;
    winDT[idx[j]]  = (dt0s2[j] < dt1s2[j]) ? (dt0s2[j] + eswt) : (dt1s2[j] + eswt);
  }

  arma::vec RT  = winDT + t0; // Add t0
  arma::mat out = check_finite(RT, R, debug);
  return out;
}

//' @rdname rplba1R
//' @export
// [[Rcpp::export]]
arma::mat rplba3(unsigned int n, arma::vec A, arma::vec B, arma::vec C,
  arma::vec mean_v, arma::vec mean_w, arma::vec sd_v, arma::vec sd_w,
  double rD, double tD, double swt, double t0) {

  double swt_r = rD + swt;  // rate delay + switch
  double swt_b = tD + swt; // thre delay + switch
  double b0    = A[0] + B[0];
  double b1    = A[1] + B[1];
  double c0    = b0 + C[0];           // changed threshold for acc 1
  double c1    = b1 + C[1];           // changed threshold for acc 2
  double swt1, swt2;                  // swt mutation
  bool a0=false, a1=false, a2=false;

  // Determine which switching condition occurs, depending on rate delay and
  // threshold delay. Both parameters are fed from a sampler.
  if (swt_r == swt_b) {       // condition 0: rate and thresold change co-occur
    a0 = true;
    swt1 = swt_r;
    swt2 = swt_r;
  } else if (swt_b < swt_r) { // condition 1: threshold change occurs early
    a1 = true;
    swt1 = swt_b;
    swt2 = swt_r;
  } else { // (swt_b > swt_r); condition 2: rate change occurs early
    a2 = true;
    swt1 = swt_r;
    swt2 = swt_b;
  }

  double swtD = swt2 - swt1;
  // --------------------------------------------------------------------------
  double u0, u1, v0, v1, w0, w1, Y0, Y1, Z0, Z1;
  double dt0_stage1, dt0_stage2, dt0_stage3, dt1_stage1, dt1_stage2, dt1_stage3;
  double DT_tmp1, DT_tmp2, DT_tmp3;
  double x0, x1;
  unsigned int R_tmp1,  R_tmp2,  R_tmp3;
  arma::vec RT(n), R(n);

  for(size_t i=0; i<n; i++)
  {
    // Stage 1
    u0 = rtn_scalar(mean_v[0], sd_v[0], 0, INFINITY) ;
    u1 = rtn_scalar(mean_v[1], sd_v[1], 0, INFINITY) ;
    // x0 = pVec[0] * Rf_runif(0, 1); // Stage 1 starting pos choice 0
    // x1 = pVec[1] * Rf_runif(0, 1); // Stage 1 starting pos choice 1
    x0 = Rf_runif(0, A[0]); // Stage 1 starting pos choice 0
    x1 = Rf_runif(0, A[1]); // Stage 1 starting pos choice 1
    dt0_stage1 = (b0 - x0) / u0;
    dt1_stage1 = (b1 - x1) / u1;
    DT_tmp1 = dt0_stage1 < dt1_stage1 ? dt0_stage1: dt1_stage1;
    R_tmp1  = dt0_stage1 < dt1_stage1 ? 1 : 2;

    // Stage 2
    // If rate changes earlier than threshold, threshold stays [b0 b1],
    // otherwise threshold changes eariler or equal to rate changes
    Y0 = a2 ? b0 - (x0 + swt1*u0) : c0 - (x0 + swt1*u0);
    Y1 = a2 ? b1 - (x1 + swt1*u1) : c1 - (x1 + swt1*u1);
    // If threshold changes earlier than rate, rate stays [u0 u1]
    // otherwise rate changes eariler or equal to threshold changes
    v0 = a1 ? u0 : rtn_scalar(mean_w[0], sd_w[0], 0, INFINITY);
    v1 = a1 ? u1 : rtn_scalar(mean_w[1], sd_w[1], 0, INFINITY);
    dt0_stage2 = Y0 / v0;
    dt1_stage2 = Y1 / v1;
    DT_tmp2 = dt0_stage2 < dt1_stage2 ? dt0_stage2 + swt1 : dt1_stage2 + swt1;
    R_tmp2  = dt0_stage2 < dt1_stage2 ? 1 : 2;

    // Stage 3
    Z0 = (a1) ? Y0 - swtD*v0 : (Y0 - swtD*v0) + C[0];
    Z1 = (a1) ? Y1 - swtD*v1 : (Y1 - swtD*v1) + C[1];
    w0 = (a1) ? rtn_scalar(mean_w[0], sd_w[0], 0, INFINITY) : v0;
    w1 = (a1) ? rtn_scalar(mean_w[1], sd_w[1], 0, INFINITY) : v1;
    dt0_stage3 = Z0 / w0;
    dt1_stage3 = Z1 / w1;
    DT_tmp3 = dt0_stage3 < dt1_stage3 ? dt0_stage3 + swt2 : dt1_stage3 + swt2;
    R_tmp3  = dt0_stage3 < dt1_stage3 ? 1 : 2;

    if (DT_tmp1 <= swt1) {
      RT[i] = DT_tmp1 + t0;
      R[i]  = R_tmp1;
    } else if (DT_tmp2 <= swt2 || a0) {
      RT[i] = DT_tmp2 + t0;
      R[i]  = R_tmp2;
    } else {
      RT[i] = DT_tmp3 + t0;
      R[i]  = R_tmp3;
    }
  }


  arma::mat out = arma::join_horiz(RT, R);
  return out;
}

//' @rdname rplba1R
//' @export
// [[Rcpp::export]]
arma::vec n1PDF_plba1(arma::vec x, arma::vec A, arma::vec b, arma::vec mean_v,
                      arma::vec sd_v, double t0, arma::vec mean_w,
                      double rD, double swt, unsigned int n, double h,
                      unsigned int ncore, bool debug) {

  arma::mat sim = rplba1(n, A, b, t0, mean_v, mean_w, sd_v, rD, swt, ncore, debug);
  arma::vec sRT = sim.col(0);
  arma::vec sR  = sim.col(1);
  arma::vec RT0 = sRT.elem(arma::find(sR==1));
  arma::vec out = spdf(x, RT0, n, h, debug);
  return out;
}

//' @rdname rplba1R
//' @export
// [[Rcpp::export]]
arma::vec n1PDF_plba2(arma::vec x, arma::vec A, arma::vec b, arma::vec mean_v,
                      arma::vec sd_v, double t0, arma::vec mean_w, arma::vec sd_w,
  double rD, double swt, unsigned int n, double h, unsigned int ncore, bool debug) {

  arma::mat sim = rplba2(n, A, b, t0, mean_v, mean_w, sd_v, sd_w, rD, swt, ncore, debug);
  arma::vec sRT = sim.col(0);
  arma::vec sR  = sim.col(1);
  arma::vec RT0 = sRT.elem(arma::find(sR==1));
  arma::vec out = spdf(x, RT0, n, h, debug);
  return out;
}

//' @rdname rplba1R
//' @export
// [[Rcpp::export]]
arma::vec n1PDF_plba3(arma::vec x, unsigned int n, arma::vec A, arma::vec B,
  arma::vec C, arma::vec mean_v, arma::vec sd_v, arma::vec mean_w,
  arma::vec sd_w, double rD, double tD, double swt, double t0, double h) {

  arma::mat sim = rplba3(n, A, B, C, mean_v, mean_w, sd_v, sd_w, rD, tD,
    swt, t0);
  arma::vec sRT = sim.col(0);
  arma::vec sR  = sim.col(1);
  arma::vec RT0 = sRT.elem(arma::find(sR == 1));
  arma::vec out = spdf(x, RT0, n, h, false);
  return out;

}


