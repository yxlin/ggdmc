#include <ggdmc.hpp>

using namespace Rcpp;

//' Generate Random Choice RT Data from LNR Model
//'
//' Race among \code{nacc} accumulators, using log-normal race model
//'
//' @param n numbers of observation
//' @param meanlog a n_acc length vector or a n_acc x n matrix. mean
//' of the distribution on the log scale without default value
//' @param sdlog a n_acc length vector or a n_acc x n matrix. Standard
//' deviation of the distribution on the log scale without default value.
//' @param t0 a scalar, a vector of length number of accumulators or a
//' matrix with 1 row per accumulator, when start time differ on each trial
//' @param st0 range of non-decision time variability, must be a scalar, as the same
//' variability is assumed in a common encoding/production stage
//' @return a matrix
//' @examples
//' ## A simple demo
//' pmat <- matrix(c(-1, 0, 1, 1, .2, .2, 0, 0), 2)
//' set.seed(123)
//' dat0 <- rlnr(4, pmat[,1], pmat[,2], pmat[,3], pmat[1,4])
//' ##           [,1] [,2]
//' ## [1,] 0.4100361    0
//' ## [2,] 0.4922407    0
//' ## [3,] 1.7855260    1
//' ## [4,] 0.4822220    1
//' ##
//' ## Three accumulators
//' n <- 1e5
//' meanlog <- c(.5, .75, 1);
//' sdlog <- c(1,1,1)
//' t0 <- c(.2,1,1)
//' set.seed(123)
//' dat1 <- rlnr(n, meanlog, sdlog, t0)
//' table(dat1[,2])
//' hist(dat1[,1], breaks = "fd", main = "", xlab ="")
//'
//' ## t0 has one element only
//' t0 <- .2
//' set.seed(123)
//' dat2 <- rlnr(n, meanlog, sdlog, t0)
//' table(dat2[,2])
//' hist(dat2[,1], breaks = "fd", freq = FALSE, main = "", xlab ="")
//' ## check t0 noise
//' st0 <- 1
//' set.seed(123)
//' dat3 <- rlnr(n, meanlog, sdlog, t0, st0)
//' table(dat3[,2])
//' hist(dat3[,1], breaks = "fd", freq = FALSE, main = "", xlab ="")
//'
//' @export
// [[Rcpp::export]]
arma::mat rlnr(unsigned int n, arma::vec meanlog, arma::vec sdlog, arma::vec t0,
  double st0 = 0) {

  unsigned int n_acc = meanlog.n_elem;
  arma::mat RT(n_acc, n), out;
  if (t0.n_elem == 1) t0 = arma::repmat(t0, n_acc, 1);

  for (size_t i = 0; i < n_acc; i++) {
      for(size_t j = 0; j < n; j++) RT(i, j) = R::rlnorm(meanlog(i), sdlog(i)) + t0(i);
  }

  arma::urowvec idx = arma::index_min(RT, 0);
  arma::colvec winner = arma::conv_to< arma::colvec >::from(idx);
  arma::vec fastRT(n);
  for (size_t k = 0; k < n; k++) fastRT(k) = RT(idx(k), k);

  if (st0 != 0) {
    arma::vec RTst0 = fastRT + st0 * arma::randu(n);
    out = arma::join_horiz(RTst0, winner);
  } else {
    out = arma::join_horiz(fastRT, winner);
  }

  return out;
}

// [[Rcpp::export]]
DataFrame rlnrDF(unsigned int n, arma::vec meanlog, arma::vec sdlog,
  arma::vec t0, double st0 = 0) {

    unsigned int n_acc = meanlog.n_elem;
    arma::mat RT(n_acc, n);
    DataFrame out;

    if (t0.n_elem == 1) t0 = arma::repmat(t0, n_acc, 1);

    for (size_t i = 0; i < n_acc; i++) {
        for(size_t j = 0; j < n; j++) RT(i, j) = R::rlnorm(meanlog(i), sdlog(i)) + t0(i);
    }

    arma::urowvec idx = arma::index_min(RT, 0);
    arma::colvec winner = arma::conv_to< arma::colvec >::from(idx);
    arma::vec fastRT(n);
    for (size_t k = 0; k < n; k++) fastRT(k) = RT(idx(k), k);

    if (st0 != 0) {
        arma::vec RTst0 = fastRT + st0 * arma::randu(n);
        out = DataFrame::create(Named("RT") = RTst0,
                                Named("R")  = 1 + winner);
    } else {
        out = DataFrame::create(Named("RT") = fastRT,
                                Named("R")  = 1 + winner);
    }
    return out;
}


// [[Rcpp::export]]
arma::vec n1PDFfixedt0_lnr1(arma::mat x, arma::vec meanlog, arma::vec sdlog)
{
  int n     = x.n_cols;
  int n_acc = x.n_rows;
  arma::vec density(n), out(n);

  // The first row must be the shortest RT
  for (int i = 0; i < n; i++) out(i) = R::dlnorm(x(0, i), meanlog(0), sdlog(0), false);

  // other accumulators
  if (n_acc > 1) {
    for (int j = 1; j < n_acc; j++) {
      for (int k = 0; k < n; k++) out(k) *= R::plnorm(x(j, k), meanlog(j), sdlog(j), false, false);
    }
  }
  return out;
}

// [[Rcpp::export]]
arma::vec n1PDFfixedt0_lnr2(arma::mat x, arma::mat meanlog, arma::mat sdlog)
{
  int n     = x.n_cols;
  int n_acc = x.n_rows;
  arma::vec density(n), out(n);

  // The first row must be the shortest RT
  for (int i = 0; i < n; i++) out(i) = R::dlnorm(x(0, i), meanlog(0, i), sdlog(0, i), false);

  // other accumulators
  if (n_acc > 1) {
    for (int j = 1; j < n_acc; j++) {
      for (int k = 0; k < n; k++) out(k) *= R::plnorm(x(j, k), meanlog(j, k), sdlog(j, k), false, false);
    }
  }
  return out;
}
