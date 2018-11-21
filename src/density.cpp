#include <ggdmc.hpp>
#include <regex>

using namespace Rcpp;

arma::mat transform_ddm(arma::mat parmat, std::string cell,
  std::vector<std::string> dim1, arma::uvec isr1) {
  //   0  1  2   3  4   5   6  7    8   9
  //   a  v  zr  d  sz  sv  t0 st0  RT  precision (dmc uses names)
  for (size_t j = 0; j < dim1.size(); j++)
    if (dim1[j] == cell && isr1(j)) parmat.row(2) = 1 - parmat.row(2);
  return parmat;
}

arma::mat transform_norm(arma::mat parmat) {
  // A b  t0 mean_v sd_v st0
  parmat.row(1) += parmat.row(0); // + parmat.row(1); // calculate b = A + B
  return parmat;
}

arma::mat transform_plba0(arma::mat parmat) {
  parmat.row(1) = parmat.row(0) + parmat.row(1); // calculate b = A + B
  parmat.row(3).swap_cols(0, 1); // correct mean_w position
  parmat.swap_rows(3, 4);
  arma::mat out = parmat;
  return out;
}

arma::mat transform_plba1(arma::mat parmat) {
  // Rcout << "Before swap parmat \n" << arma::trans(parmat) << std::endl;
  //   0 A       0.75   0.75
  //   1 B       0.25   0.25
  //   2 mean_v  2.50   1.50
  //   3 mean_w  2.40   1.20
  //   4 sd_v    1.00   1.00
  //   5 rD      0.10   0.10
  //   6 t0      0.20   0.20
  //   7 swt   489.00 489.00

  parmat.row(1) = parmat.row(0) + parmat.row(1); // calculate b = A + B
  parmat.row(3).swap_cols(0, 1); // correct mean_w position
  parmat.swap_rows(3, 4);
  arma::mat out = parmat;
  return out;
}

arma::mat transform_plba2(arma::mat parmat) {
  // A b mean_v mean_w sd_v sd_w rD t0 swt
  parmat.row(1)  = parmat.row(0) + parmat.row(1);
  parmat.row(3).swap_cols(0, 1);
  return parmat;
}

arma::mat transform_plba3(arma::mat parmat) {
  parmat.row(4).swap_cols(0, 1);
  parmat.row(6).swap_cols(0, 1);
  return parmat;
}

arma::mat transform_cnorm(arma::mat parmat) {
  parmat.row(1) = parmat.row(0) + parmat.row(1);
  parmat.swap_rows(2, 5);
  parmat.swap_rows(3, 5);
  parmat.swap_rows(4, 5);
  for (size_t i = 0; i < parmat.row(5).n_elem; i++) {
    parmat(5, i) = 2 * R::plogis(parmat(5, i), 0, 1, true, false) - 1;
  }
  return parmat;
}


// arma::mat flipresp(arma::mat parmat, arma::mat n1order,
//   std::vector<std::string> dim1, std::vector<std::string> dim3,
//   std::vector<std::string> parnames, std::string cell) {
//
//   int nparname = parnames.size();
//   int nr       = dim3.size();
//   arma::mat out(nparname, nr);
//   arma::vec n1idx;
//
//   if(parmat.size() == 0) Rcpp::stop("parmat error in flipresp");
//   for (int i = 0; i < dim1.size(); i++) // s1.f1.r1, s2.f1.r1
//   {
//       if (dim1[i] == cell) {
//         n1idx = arma::vectorise(n1order.row(i));
//         for (int j = 0; j < nr; j++) {
//             out.col(j) = parmat.col(n1idx.at(j) - 1);
//         }
//       }
//   }
//
//   return out;
// }


// [[Rcpp::export]]
arma::mat FlipResponse_norm(arma::mat parmat, arma::umat n1mat,
                   std::vector<std::string> dim1,
                   std::vector<std::string> dim3,
                   std::vector<std::string> parnames, std::string cell) {

  unsigned int nparname = parnames.size();
  unsigned int nr       = dim3.size();
  arma::mat out(nparname, nr);
  arma::uvec n1idx;

  if (parmat.size() == 0) Rcpp::stop("parmat error in flipresp");
  for (size_t i = 0; i < dim1.size(); i++) // s1.f1.r1, s2.f1.r1
  {
    if (dim1[i] == cell) {
      n1idx = arma::vectorise(n1mat.row(i));
      for (size_t j = 0; j < nr; j++) {
        out.col(j) = parmat.col(n1idx.at(j) - 1);
      }
    }
  }

  return out;
}


//' Check parameter vector of DDM
//'
//' Check if all DDM parameters are within reasonable ragnes
//'
//' @param pVec a numeric vector storing the DDM parameters
//' @export
// [[Rcpp::export]]
bool checkddm2(std::vector<double> pVec) {
  // pVec is ordered as Voss & Voss's density.c wihtout RT.
  // ----------Name sequence correction----------
  // 0  1  2  3  4    5   6  7    8   9
  // a  v zr  d  szr  sv  t0 st0  RT  precision;
  // a  v  z  d  sz   sv  t0 st0  RT  precision (dmc uses names)
  // double aLimit = 2.0 ;  // Use Heathcote limits; Voss=10.0; Verdonck=1.0
  // double vLimit = 5.0 ;  // Voss = 50.0
  double zrLower = pVec[2] - 0.5*pVec[4];
  double zrUpper = pVec[2] + 0.5*pVec[4];
  double szLimit = 1.999999*std::min(pVec[2], 1 - pVec[2]) ;

  // 1. a cannot be less than z, or less than (or equal to) 0 or greater than
  // its limit
  bool c1 = pVec[0] < 0 ;
  // 2. If v's absolute value greater than vLimit (Voss's 50; Andrew's 5)
  //bool c2 = std::abs(pVec[1]) > vLimit ;
  // 3. Check zr. If zLower less than 1e-6 or zUpper greater than 1
  bool c3 = (zrLower < 1e-6 || zrUpper > 0.999999) ;  // 0 ~ 1
  // 4. If t0 - abs(d)/2 - st0/2 is less than 0
  bool c4 = (( pVec[6] - 0.5*std::abs(pVec[3]) - 0.5*pVec[7] ) < 0) ;
  // 5. If t0 is almost 0
  bool c5 = pVec[6] < 1e-6 ;
  // 6. If szr greater than 1 or less than 0 or greater than 2 x min(z)
  bool c6 = (pVec[4]  < 0) || pVec[4] > szLimit ;
  // 7. If sv greater than 2 or less than 0
  // bool c7 = pVec[5] > 2 || pVec[5] < 0;  // affect sv profile
  // 8. If st0 less than 0
  bool c8 = pVec[7] < 0 ;
  // 9. If t0 - st0/2 less than 0
  bool c9 = (pVec[6] - pVec[7]/2) < 0 ;

  bool out = c1 || c3 || c4 || c5 || c6 || c8 || c9;
  //bool out = c1 || c4 || c5 || c8 || c9 ;

  return out;
}


//' @rdname likelihood_norm
//' @export
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
arma::mat p_df(arma::vec pVec, std::string cell,
  std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::uvec isr1, arma::umat n1idx,
  bool n1order = true) {
  // should remove cpp11

  unsigned int nc = dim1.size() ;  // number of condition, eg s1.r1 etc.
  unsigned int np = dim2.size() ;  // number of parameter x condition, eg v.f1, v.f2
  unsigned int nr = dim3.size() ;  // number of accumualtor/response, eg r1, r2

  unsigned int nparname = parnames.size();
  arma::vec tmp(nparname); // user enter parameter names, no condition, but with constant
  arma::mat parmat(nparname, nr); // npar x nr

  // Iterate through accumulators, eg r1, r2, (if have any), r3, etc.
  // parnames is the parameter the user wishes to fit, which may/may not
  // differ from the parameter likelihood function expects
  // eg A, B, mean_v, sd_v, t0, st0 vs A, b, mean_v, sd_v, t0, st0
  // likelihood function expects to see b, but the user is asked to enter B.
  for (size_t i = 0; i < nr; i ++) // r1, r2, r3, etc.
  {
    for (size_t j = 0; j  < nc; j++) // eg s1.r1, s2.r1, s1.r2, s2.r2 etc.
    {
      if (dim1[j] == cell)
      {
        size_t ii = 0;
        for (size_t k = 0; k < np; k++) // eg allpars a, v.f1, v.f2, z, d, sz) ...
        {
          if (model(j, k, i))
          { // allpar is constant or NA, not in constant. allpar names == dim2 names
            tmp[ii] = allpar[k]; // constant values or NA
            for(size_t l = 0; l < pVec.size(); l++)
            { // replace NA with values in p.vector.
              if (pnames[l] == dim2[k] && std::isnan(tmp[ii])) tmp[ii] = pVec[l];
            }
            ii++ ;
          }
        }
      }
    }
    parmat.col(i) = tmp;
  }

  if (type == "rd") { parmat = transform_ddm(parmat, cell, dim1, isr1); }
  if (type == "norm" || type=="norm_pda" || type=="norm_pda_gpu") { parmat = transform_norm(parmat); }
  if (n1order) { parmat = FlipResponse_norm(parmat, n1idx, dim1, dim3, parnames, cell); }
  if (type == "plba0_gpu") { parmat = transform_plba0(parmat); }
  if (type == "plba1" || type == "plba1_gpu") { parmat = transform_plba1(parmat); }
  if (type == "plba2") { parmat = transform_plba2(parmat); }
  if (type == "plba3") { parmat = transform_plba3(parmat); }
  if (type == "cnorm") { parmat = transform_cnorm(parmat); }
  return arma::trans(parmat);
}


// [[Rcpp::export]]
arma::uvec getbounds(List data) {

  NumericVector modelAttr = data.attr("model");
  List matchMapOuter  = modelAttr.attr("match.map"); // Handle two-layer list
  List matchMapInner  = matchMapOuter["M"] ;      // "M" as match.map's name
  List modelDimension = modelAttr.attr("dimnames") ;
  std::vector<std::string> stimulusName   = matchMapInner.names() ;
  std::vector<std::string> dim1 = modelDimension[0] ;
  arma::uvec bound (dim1.size()); bound.zeros();
  unsigned int nmatch = matchMapInner.size();

  for (size_t i = 0; i < nmatch; i++) {
    std::string b = matchMapInner[i];  // *it == matchMapInner[i]; r1, r2
    std::regex matchPattern(stimulusName[i] + "(.*)" + b);

    for(size_t j = 0; j < dim1.size(); j++)
      if (bound(j) == 0) bound(j) = (std::regex_match(dim1[j], matchPattern)) ? 1 : 0;
  }
  return bound;
}

// [[Rcpp::export]]
arma::vec density_rd(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1, std::vector<std::string> dim2,
  std::vector<std::string> dim3, arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT, arma::uvec matchcell, arma::uvec isr1,
  double precision = 2.5) {

  std::vector<double> gVec(10);
  gVec.at(9) = precision;

  arma::mat pmat;
  arma::vec selectedRT, out(RT.n_elem);
  arma::uvec RTIdx;

  // a   v   zr   d  szr  sv  t0 st0 RT, precision;
  for (size_t i = 0; i < ise.n_elem; i++)
  {
    if (!ise[i])
    {
      pmat = p_df(pVec, dim1[i], pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, isr1, n1idx, true);
      for(size_t j = 0; j < 8; j++) gVec[j] = pmat(0, j);
      RTIdx = arma::find(cellidx.col(i) == 1);
      selectedRT = RT(RTIdx);

      if (matchcell(i))  // choose g_plus
      {
        for (size_t k = 0; k < selectedRT.n_elem; k++) {
          gVec.at(8) = selectedRT(k);
          out(RTIdx(k)) = checkddm2(gVec) ? 1e-10 : std::max(std::abs(g_plus(gVec)), 1e-10);
        }
      } else {
        for (size_t k = 0; k < selectedRT.n_elem; k++) {
          gVec.at(8) = selectedRT(k);
          out(RTIdx(k)) = checkddm2(gVec) ? 1e-10 : std::max(std::abs(g_minus(gVec)), 1e-10);
        }
      }
    }
  }

  return out;
}

// [[Rcpp::export]]
arma::vec density_norm(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift)
{
  arma::mat pmat;
  arma::vec out(RT.n_elem);
  arma::uvec RTIdx, valid_rows;

  // arma::vec n1PDFfixedt0(arma::vec rt, arma::vec A, arma::vec b, arma::vec mean_v,
  //                        arma::vec sd_v, arma::vec t0, bool posdrift) {

  // pmat: A   b  t0 mean_v sd_v st0
  for (size_t i = 0; i < ise.n_elem; i++) {
    if (!ise(i)) {
      pmat = p_df(pVec, dim1[i], pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, isr1, n1idx, true);
      RTIdx = arma::find(cellidx.col(i) == 1);
      valid_rows = arma::find_finite(pmat.col(1));
      arma::mat pmat_ = pmat.rows(valid_rows);
      out(RTIdx) = n1PDFfixedt0(RT(RTIdx), pmat_.col(0), pmat_.col(1),
          pmat_.col(3), pmat_.col(4), pmat_.col(2), posdrift);
    }
  }

  for (size_t j = 0; j < out.n_elem; j++) out(j) = std::max(out(j), 1e-10);
  return out;
}


// [[Rcpp::export]]
arma::vec density_cnorm_pda(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1,
  unsigned int nsim, double bw, bool debug = false) {
  arma::mat pmat;
  arma::vec out(RT.n_elem);
  arma::uvec RTIdx, valid_rows;

  for (size_t i = 0; i < ise.n_elem; i++) {
    if (!ise(i)) {
      pmat = p_df(pVec, dim1[i], pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, isr1, n1idx, true);
      RTIdx           = arma::find(cellidx.col(i) == 1);
      valid_rows      = arma::find_finite(pmat.col(1));
      arma::mat pmat_ = pmat.rows(valid_rows);

      //         A   b  t0 mean_v sd_v    corr_v st0
      // LEFT  1.5 2.7 0.2      3    1 0.3363755   0
      // RIGHT 1.5 2.7 0.2      2    1 0.3363755   0

      out(RTIdx) = n1PDF_cnorm(RT(RTIdx), pmat_.col(0), pmat_.col(1),
        pmat_.col(2), pmat_.col(3), pmat_.col(4), pmat_.col(6),
        pmat_(1,6), nsim, bw, debug);
    }
  }

  for (size_t j = 0; j < out.n_elem; j++) out(j) = std::max(out(j), 1e-10);
  return out;
}


// [[Rcpp::export]]
arma::vec density_norm_pda(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int npda = 16384,
  double bw = .01, bool debug = false)
{
  arma::mat pmat;
  arma::vec out(RT.n_elem);
  arma::uvec RTIdx;

  // A   b  t0 mean_v sd_v st0
  for (size_t i = 0; i < ise.n_elem; i++)
  {
      if (!ise(i))
      {
        pmat = p_df(pVec, dim1[i], pnames, allpar, parnames, model, type,
                      dim1, dim2, dim3, isr1, n1idx, true);

          if (arma::any(pmat.col(4) < 0)) {
            if (debug) Rcout << "sd_v negative\n";
            out.fill(1e-10);
          } else {
            RTIdx = arma::find(cellidx.col(i) == 1);
            out(RTIdx) = n1PDFfixedt0_pda(RT(RTIdx), pmat(0, 0), pmat(0, 1),
                  pmat.col(3), pmat.col(4), pmat(0, 2), npda, bw, debug);
          }
      }
  }

  for (size_t j = 0; j < out.n_elem; j++) out(j) = std::max(out(j), 1e-10);

  return out;
}

// [[Rcpp::export]]
arma::vec density_norm_gpu(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int nsim = 16384,
  double bw = .01, unsigned int gpuid = 0, unsigned int nthread = 32,
  bool debug = false)
{
  Rcpp::NumericVector tmp;
  arma::mat pmat;
  arma::vec out(RT.n_elem);
  arma::uvec RTIdx;
  // A   b  t0 mean_v sd_v st0
  for (size_t i = 0; i < ise.n_elem; i++)
  {
    if (!ise(i))
    {
      pmat = p_df(pVec, dim1[i], pnames, allpar, parnames, model, type,
         dim1, dim2, dim3, isr1, n1idx, true);

      if (arma::any(pmat.col(4) < 0)) {
         if (debug) Rcout << "sd_v negative\n"; out.fill(1e-10);
      } else {
         RTIdx = arma::find(cellidx.col(i) == 1);
         tmp   = n1PDF_gpu(RT(RTIdx), pmat(0, 0), pmat(0, 1),
          pmat.col(3), pmat.col(4), pmat(0, 2), nsim, nthread, gpuid, bw, debug);
         arma::vec armatmp = as<arma::vec>(tmp);
         out(RTIdx) = armatmp;
      }
    }
  }
  return out;
}


// [[Rcpp::export]]
arma::vec density_plba1(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int nsim = 16384,
  double bw = .01, unsigned int ncore = 1, bool debug = false)
{
  //   A    b mean_v sd_v mean_w  rD  t0 swt
  arma::mat pmat;
  arma::vec out(RT.n_elem);
  arma::uvec RTIdx;

  for (size_t i = 0; i < ise.n_elem; i++)
  {
    if (!ise(i))
    {
      pmat = p_df(pVec, dim1[i], pnames, allpar, parnames, model, type,
        dim1, dim2, dim3, isr1, n1idx, true);
      if (arma::any(pmat.col(4) < 0)) {
        if (debug) Rcout << "sd_v negative\n"; out.fill(1e-10);
      } else {
        RTIdx = arma::find(cellidx.col(i) == 1);
        out(RTIdx) = n1PDF_plba1(RT(RTIdx), pmat.col(0), pmat.col(1),
          pmat.col(2), pmat.col(3), pmat(0,6), pmat.col(4), pmat(0,5),
          pmat(0,7), nsim, bw, ncore, debug);
      }
    }
  }

  for (size_t j = 0; j < out.n_elem; j++) out(j) = std::max(out(j), 1e-10);
  return out;
}

// [[Rcpp::export]]
arma::vec density_plba0_gpu(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int nsim = 16384,
  double bw = .01, unsigned int ncore = 1, unsigned int gpuid = 0,
  unsigned int nthread = 32, bool debug = false)
{

  Rcpp::NumericVector tmp;
  arma::vec out(RT.n_elem);
  arma::uvec RTIdx;

  for (size_t i = 0; i < ise.n_elem; i++)
  {
    if (!ise(i))
    {
      //if (debug) Rcout << "In density pVec:\n " << pVec << std::endl;
      arma::mat pmat = p_df(pVec, dim1[i], pnames, allpar, parnames, model, type,
        dim1, dim2, dim3, isr1, n1idx, true);
      //if (debug) Rcout << "after pmat:\n " << pmat << std::endl;


      if (arma::any(pmat.col(0) < 0) || arma::any(pmat.col(1) < 0) ||
        arma::any(pmat.col(3) < 0) || arma::any(pmat.col(5) < 0) ||
        arma::any(pmat.col(6) < 0) || arma::any(pmat.col(7) < 0)) {

        if (debug && arma::any(pmat.col(3) < 0)) Rcout << "sd_v negative\n";
        if (debug && arma::any(pmat.col(0) < 0)) Rcout << "A negative\n";
        if (debug && arma::any(pmat.col(1) < 0)) Rcout << "b negative\n";
        if (debug && arma::any(pmat.col(5) < 0)) Rcout << "rD negative\n";
        if (debug && arma::any(pmat.col(6) < 0)) Rcout << "t0 negative\n";
        if (debug && arma::any(pmat.col(7) < 0)) Rcout << "swt negative\n";
        out.fill(1e-10);
      } else {
        RTIdx = arma::find(cellidx.col(i) == 1);
        //  A  b mean_v sd_v mean_w  rD  t0 swt
        //  0  1      2    3      4   5   6   7
        tmp = n1PDF_plba0_gpu(RT(RTIdx), pmat(0,0), pmat(0,1),
          pmat.col(2), pmat.col(3), pmat(0,6), pmat.col(4), pmat(0,5),
          pmat(0,7), nsim, nthread, gpuid, bw, debug);

        // Rcpp::NumericVector n1PDF_plba1_gpu(arma::vec x, double A, double b,
        // arma::vec mean_v,
        //   arma::vec sd_v, double t0, arma::vec mean_w, double rD, double swt,
        //   int n, int nthread, int gpuid, double bw, bool debug)
        arma::vec armatmp = as<arma::vec>(tmp);
        out(RTIdx) = armatmp;

      }
    }
  }

  for (size_t j = 0; j < out.n_elem; j++) out(j) = std::max(out(j), 1e-10);
  return out;
}


// [[Rcpp::export]]
arma::vec density_plba1_gpu(arma::vec pVec, std::vector<std::string> pnames,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, unsigned int nsim = 16384,
  double bw = .01, unsigned int ncore = 1, unsigned int gpuid = 0,
  unsigned int nthread = 32, bool debug = false)
{

  Rcpp::NumericVector tmp;
  arma::vec out(RT.n_elem);
  arma::uvec RTIdx;

  for (size_t i = 0; i < ise.n_elem; i++)
  {
    if (!ise(i))
    {
      //if (debug) Rcout << "In density pVec:\n " << pVec << std::endl;
      arma::mat pmat = p_df(pVec, dim1[i], pnames, allpar, parnames, model, type,
        dim1, dim2, dim3, isr1, n1idx, true);
      //if (debug) Rcout << "after pmat:\n " << pmat << std::endl;


      if (arma::any(pmat.col(0) < 0) || arma::any(pmat.col(1) < 0) ||
          arma::any(pmat.col(3) < 0) || arma::any(pmat.col(5) < 0) ||
          arma::any(pmat.col(6) < 0) || arma::any(pmat.col(7) < 0)) {

        if (debug && arma::any(pmat.col(3) < 0)) Rcout << "sd_v negative\n";
        if (debug && arma::any(pmat.col(0) < 0)) Rcout << "A negative\n";
        if (debug && arma::any(pmat.col(1) < 0)) Rcout << "b negative\n";
        if (debug && arma::any(pmat.col(5) < 0)) Rcout << "rD negative\n";
        if (debug && arma::any(pmat.col(6) < 0)) Rcout << "t0 negative\n";
        if (debug && arma::any(pmat.col(7) < 0)) Rcout << "swt negative\n";
        out.fill(1e-10);
      } else {
        RTIdx = arma::find(cellidx.col(i) == 1);
        //  A  b mean_v sd_v mean_w  rD  t0 swt
        //  0  1      2    3      4   5   6   7
        tmp = n1PDF_plba1_gpu(RT(RTIdx), pmat(0,0), pmat(0,1),
          pmat.col(2), pmat.col(3), pmat(0,6), pmat.col(4), pmat(0,5),
          pmat(0,7), nsim, nthread, gpuid, bw, debug);

        // Rcpp::NumericVector n1PDF_plba1_gpu(arma::vec x, double A, double b,
        // arma::vec mean_v,
        //   arma::vec sd_v, double t0, arma::vec mean_w, double rD, double swt,
        //   int n, int nthread, int gpuid, double bw, bool debug)
        arma::vec armatmp = as<arma::vec>(tmp);
        out(RTIdx) = armatmp;

      }
    }
  }

  for (size_t j = 0; j < out.n_elem; j++) out(j) = std::max(out(j), 1e-10);
  return out;
}


//' Calculate Summed, Log-likelihood of a Cognitive Model
//'
//' The function calculates log-likelihood for every trial.  The input must
//' be a data model instance.
//'
//' @param pVec a parameter vector
//' @param pnames a string vector storing the name of a parameter vector
//' @param allpar all parameters
//' @param parnames parameter names
//' @param model a model specification
//' @param type model type
//' @param dim1 first dimension of a model
//' @param dim2 second dimension of a model
//' @param dim3 third dimension of a model
//' @param n1idx n1 order index
//' @param ise an index vector storing if a cell is empty.
//' @param cellidx cell index
//' @param RT a RT vector
//' @param matchcell an index vector storing is the cell is a match response
//' @param isr1 is r1 index
//' @param posdrift a Boolean switch to enforce postive drift rate correction
//' @param nsim number of simulation
//' @param bw bandwidth
//' @param ncore number of parallel core
//' @param gpuid GPU card index on a multiple GPU machine
//' @param debug whether to print debugging information
//' for assuming drift rates are drawn from a normal distribution.
//' @return a double scalar
//' @examples
//' m1 <- BuildModel(
//'   p.map     = list(a = "1", v = "1", z = "1", d = "1", sz = "1", sv = "1",
//'               t0 = "1", st0 = "1"),
//'   constants = c(st0 = 0, d = 0),
//'   match.map = list(M = list(s1 = "r1", s2 = "r2")),
//'   factors   = list(S = c("s1", "s2")),
//'   responses = c("r1", "r2"),
//'   type      = "rd")
//'
//' p.vector <- c(a = 1, v = 1, z = 0.5, sz = 0.25, sv = 0.2, t0 = .15)
//'
//' ## Set up a model-data instance
//' dat <- simulate(m1, 128, ps = p.vector)
//' dmi <- BuildDMI(dat, m1)
//' ## sumloglike(p.vector, dmi)
//' ## [1] 0.3796048
//' @export
// [[Rcpp::export]]
double sumloglike(arma::vec pVec,
                std::vector<std::string> pnames,
                arma::vec allpar,
                std::vector<std::string> parnames,
                arma::ucube model,
                std::string type,
                std::vector<std::string> dim1,
                std::vector<std::string> dim2,
                std::vector<std::string> dim3,
                arma::umat n1idx, arma::uvec ise,
                arma::umat cellidx, arma::vec RT,
                arma::uvec matchcell, arma::uvec isr1, bool posdrift,
                unsigned int nsim, double bw, unsigned int ncore,
                unsigned int gpuid, bool debug) {
  std::string type1 ("norm");      // Available model type in ggdmc
  std::string type2 ("norm_pda");
  std::string type3 ("plba1");
  std::string type4 ("plba1_gpu");
  std::string type5 ("rd");
  std::string type6 ("plba0_gpu");
  std::string type7 ("norm_pda_gpu");
  std::string type8 ("cnorm");

  arma::vec density; density.fill(NA_REAL);

  if (type == type1) {         // LBA norm
     density = density_norm(pVec, pnames, allpar, parnames, model, type,
                            dim1, dim2, dim3, n1idx, ise, cellidx, RT,
                            matchcell, isr1, posdrift);
  } else if (type == type2) { // LBA norm_pda
     density = density_norm_pda(pVec, pnames, allpar, parnames, model, type,
                                dim1, dim2, dim3, n1idx, ise, cellidx, RT,
                                matchcell, isr1, nsim, bw, debug);
  } else if (type == type7) { // LBA norm_pda_gpu call gpda package
    unsigned int nthread = 32;
    density = density_norm_gpu(pVec, pnames, allpar, parnames, model, type,
      dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim, bw,
      gpuid, nthread, debug);

  } else if (type == type3) {
    density = density_plba1(pVec, pnames, allpar, parnames, model, type,
      dim1, dim2, dim3, n1idx, ise, cellidx, RT,
      matchcell, isr1, nsim, bw, ncore, debug);
  } else if (type == type4) {
    unsigned int nthread = 32;
    density = density_plba1_gpu(pVec, pnames, allpar, parnames, model, type,
      dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim, bw,
      ncore, gpuid, nthread, debug);

  } else if (type == type5) {
    density = density_rd(pVec, pnames, allpar, parnames, model, type,
      dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1);

  } else if (type == type6) {

    unsigned int nthread = 32;
    density = density_plba0_gpu(pVec, pnames, allpar, parnames, model, type,
      dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim, bw,
      ncore, gpuid, nthread, debug);

  } else if (type == type8) {
    density = density_cnorm_pda(pVec, pnames, allpar, parnames, model, type,
      dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, nsim, bw, debug);
  } else {
    Rcout << "Type not defined.\n";
    density.fill(1e-10);
  }

  double out = arma::accu(arma::log(density));
  if (std::isnan(out)) out = -INFINITY;
  return out;
}

// [[Rcpp::export]]
NumericVector profile_rd(arma::vec pVec,
  std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx,
  arma::uvec ise,
  arma::umat cellidx,
  arma::vec RT,
  arma::uvec matchcell,
  arma::uvec isr1,
  std::string pname,
  arma::vec ps)
{
  unsigned int npoint = ps.n_elem;
  arma::vec pVec_tmp, density;
  NumericVector out(npoint);
  pVec_tmp = pVec;
  //arma::ucube umodel = arma::conv_to<arma::ucube>::from(model);

  for (size_t i = 0; i < npoint; i++) {
    for (size_t j = 0; j < pVec.n_elem; j++) {if (pnames[j] == pname) pVec_tmp(j) = ps(i);}
    density = density_rd(pVec_tmp, pnames, allpar, parnames, model,
      type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1);
    out(i) = arma::accu(arma::log(density));
  }

  return out;
}


// [[Rcpp::export]]
NumericVector profile_norm(arma::vec pVec,
                       std::vector<std::string> pnames,
                       arma::vec allpar,
                       std::vector<std::string> parnames,
                       arma::ucube model,
                       std::string type,
                       std::vector<std::string> dim1,
                       std::vector<std::string> dim2,
                       std::vector<std::string> dim3,
                       arma::umat n1idx,
                       arma::uvec ise,
                       arma::umat cellidx,
                       arma::vec RT,
                       arma::uvec matchcell,
                       arma::uvec isr1,
                       std::string pname,
                       arma::vec ps,
                       bool posdrift)
{
  arma::vec pVec_tmp, density;
  NumericVector out(ps.n_elem);
  pVec_tmp = pVec;

  for (size_t i = 0; i < ps.n_elem; i++) {
    for (size_t j = 0; j < pVec.n_elem; j++) {if (pnames[j] == pname) pVec_tmp(j) = ps(i);}
      density = density_norm(pVec_tmp, pnames, allpar, parnames, model,
                             type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
                             posdrift);
      out(i) = arma::accu(arma::log(density));
  }

  return out;
}

// [[Rcpp::export]]
NumericVector profile_norm_pda(arma::vec pVec,
                           std::vector<std::string> pnames,
                           arma::vec allpar,
                           std::vector<std::string> parnames,
                           arma::ucube model,
                           std::string type,
                           std::vector<std::string> dim1,
                           std::vector<std::string> dim2,
                           std::vector<std::string> dim3,
                           arma::umat n1idx,
                           arma::uvec ise,
                           arma::umat cellidx,
                           arma::vec RT,
                           arma::uvec matchcell,
                           arma::uvec isr1,
                           std::string pname,
                           arma::vec ps, int nsim, double bw) {
  arma::vec pVec_tmp, density;
  NumericVector out(ps.n_elem);
  pVec_tmp = pVec;

  for (size_t i = 0; i < ps.n_elem; i++) {
    for (size_t j = 0; j < pVec.n_elem; j++) {if (pnames[j] == pname) pVec_tmp(j) = ps(i);}
    density = density_norm_pda(pVec_tmp, pnames, allpar, parnames, model,
                           type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
                           nsim, bw, false);
    out(i) = arma::accu(arma::log(density));
  }

  return out;
}

// [[Rcpp::export]]
NumericVector profile_cnorm_pda(arma::vec pVec,
  std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx,
  arma::uvec ise,
  arma::umat cellidx,
  arma::vec RT,
  arma::uvec matchcell,
  arma::uvec isr1,
  std::string pname,
  arma::vec ps, int nsim, double bw) {

  arma::vec density;
  NumericVector out(ps.n_elem);
  arma::vec pVec_tmp = pVec;

  for (size_t i = 0; i < ps.n_elem; i++) {
    for (size_t j = 0; j < pVec.n_elem; j++) {
      if (pnames[j] == pname) pVec_tmp(j) = ps(i);
    }
    density = density_cnorm_pda(pVec_tmp, pnames, allpar, parnames, model,
      type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
      nsim, bw, false);
    out(i) = arma::accu(arma::log(density));
  }

  return out;
}

// [[Rcpp::export]]
NumericVector profile_norm_gpu(arma::vec pVec,
  std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx,
  arma::uvec ise,
  arma::umat cellidx,
  arma::vec RT,
  arma::uvec matchcell,
  arma::uvec isr1,
  std::string pname,
  arma::vec ps, int nsim, double bw)
{
  arma::vec pVec_tmp, density;
  NumericVector out(ps.n_elem);
  pVec_tmp = pVec;

  for (size_t i = 0; i < ps.n_elem; i++) {
    for (size_t j = 0; j < pVec.n_elem; j++) {if (pnames[j] == pname) pVec_tmp(j) = ps(i);}
    density = density_norm_gpu(pVec_tmp, pnames, allpar, parnames, model,
      type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
      nsim, bw, false);
    out(i) = arma::accu(arma::log(density));
  }
  return out;
}

// [[Rcpp::export]]
NumericVector profile_plba1(arma::vec pVec,
  std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx,
  arma::uvec ise,
  arma::umat cellidx,
  arma::vec RT,
  arma::uvec matchcell,
  arma::uvec isr1,
  std::string pname,
  arma::vec ps, int nsim, double bw)
{
  arma::vec pVec_tmp, density;
  NumericVector out(ps.n_elem);
  pVec_tmp = pVec;
  int ncore = 1;

  for (size_t i = 0; i < ps.n_elem; i++) {
    for (size_t j = 0; j < pVec.n_elem; j++) {if (pnames[j] == pname) pVec_tmp(j) = ps(i);}
    density = density_plba1(pVec_tmp, pnames, allpar, parnames, model,
      type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
      nsim, bw, ncore, false);
    out(i) = arma::accu(arma::log(density));
  }

  return out;
}

// [[Rcpp::export]]
NumericVector profile_plba1_gpu(arma::vec pVec,
  std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx,
  arma::uvec ise,
  arma::umat cellidx,
  arma::vec RT,
  arma::uvec matchcell,
  arma::uvec isr1,
  std::string pname,
  arma::vec ps, int nsim, double bw, int gpuid, int nthread, bool debug)
{
  arma::vec pVec_tmp, density;
  NumericVector out(ps.n_elem);
  pVec_tmp = pVec;
  int ncore = 1;

  for (size_t i = 0; i < ps.n_elem; i++) {
    for (size_t j = 0; j < pVec.n_elem; j++) {if (pnames[j] == pname) pVec_tmp(j) = ps(i);}
    density = density_plba1_gpu(pVec_tmp, pnames, allpar, parnames, model,
      type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
      nsim, bw, ncore, gpuid, nthread, debug);
    out(i) = arma::accu(arma::log(density));
  }

  return out;
}

// [[Rcpp::export]]
NumericVector profile_plba0_gpu(arma::vec pVec,
  std::vector<std::string> pnames,
  arma::vec allpar,
  std::vector<std::string> parnames,
  arma::ucube model,
  std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx,
  arma::uvec ise,
  arma::umat cellidx,
  arma::vec RT,
  arma::uvec matchcell,
  arma::uvec isr1,
  std::string pname,
  arma::vec ps, int nsim, double bw, int gpuid, int nthread, bool debug)
{
  arma::vec pVec_tmp, density;
  NumericVector out(ps.n_elem);
  pVec_tmp = pVec;
  int ncore = 1;

  for (size_t i = 0; i < ps.n_elem; i++) {
    for (size_t j = 0; j < pVec.n_elem; j++) {if (pnames[j] == pname) pVec_tmp(j) = ps(i);}
    density = density_plba0_gpu(pVec_tmp, pnames, allpar, parnames, model,
      type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
      nsim, bw, ncore, gpuid, nthread, debug);
    out(i) = arma::accu(arma::log(density));
  }

  return out;
}

