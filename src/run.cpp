#include <ggdmc.hpp>
#include <chrono>
#include <random>
using namespace Rcpp;

//' Generate a Gamma Vector
//'
//' This is part of DE-MCMC algorithm. This function generates a gamma vector
//' for element-wise computation in Armadillo C++. This function is based on
//' p242 ter Braak (2006) who cited Roberts and Rosenthal (2001)
//'
//' @param npar number of parameters.
//' @param gammamult a tuning parameter stands for for gamma mutation. Default
//' value is 2.38.
//' @param hyper a boolean switch, indicating to calculate hyper gamma
//' @return a vector
//' @examples
//' pVec <- c(A = 1.51, b = 2.7, muv1 = 3.32, muv2 = 2.24, t_ND = 0.08,
//'           muw1 = 1.51, muw2 = 3.69, t_delay = 0.31, sv = 1, swt = 0.5)
//' gamma <- GetGamma(length(pVec), 2.38)
//' @export
// [[Rcpp::export]]
arma::vec GetGamma(unsigned int npar, double gammamult, bool hyper) {
  arma::vec out(npar);
  double divisor = hyper ? std::sqrt(4.0*(double)npar) : std::sqrt(2.0*(double)npar);
  for (size_t i = 0; i < npar; i++) {
    out(i) = std::isnan(gammamult) ? R::runif(0.5, 1.0) : gammamult/divisor;
  }
  return out;
}

//' Draw n other chains and shuffle them
//'
//' This is part of DE-MCMC algorithm. \code{PickChains} draws \code{n}
//' chains out of \code{length(chains)} chains, excluding the kth chain.
//' \code{GetSubchains} is used in \code{migration} operator. It draws a subset
//' of chains in \code{nchain} chains.
//'
//' \code{Getsubchains} is part of the Migration algroithms. It does two-step
//' shuffling. In step 1, it selects a number l (integer) uniformly between 1
//' and k to be the number of subpopulations for migration.  In step 2, it
//' generates the chain index (0 to nchain - 1) and lastely it shuffles them
//'
//'
//' @param k the kth processed chain. Must be an integer within the range of 0
//' to \code{nchain - 1}. No check for errorly using R index.
//' @param nchain number of chains to draw.
//' @param chains an integer vector, indicating chain index, e.g., 0:23
//' @param ngroup number of distributed groups
//' @return a column vector
//' @keywords PickChains, getsubchains
//' @export
//' @examples
//' chains <- 0:23
//'
//' ## Presuming current processing chain is the 1st chain (C index = 0)
//' ## pick 2 chains out of 24 chains, excluding current chain.
//' PickChains(0, 2, chains)
//'
//' ## Example outputs
//' ##      [,1]
//' ## [1,]   17
//' ## [2,]   12
//' ##      [,1]
//' ## [1,]    2
//' ## [2,]    5
//' ##      [,1]
//' ## [1,]    5
//' ## [2,]    3
//' ##      [,1]
//' ## [1,]   10
//' ## [2,]    8
//' ##      [,1]
//' ## [1,]   15
//' ## [2,]    8
//'
//' ## get a random number of subchains
//' GetSubchains(24)
//' ##       [,1]
//' ##  [1,]    0
//' ##  [2,]    3
//' ##  [3,]    5
//' ##  [4,]    9
//' ##  [5,]   10
//' ##  [6,]   12
//' ##  [7,]   14
//' ##  [8,]   15
//' ##  [9,]   18
//' ## [10,]   20
//' ## [11,]   21
//' ## [12,]   22
//'
//' @export
// [[Rcpp::export]]
arma::uvec PickChains(unsigned int k, unsigned int nchain, arma::uvec chains) {
  chains.shed_row(k);
  arma::uvec rchains = arma::shuffle(chains);
  return rchains.rows(0, nchain - 1);
}

//' @rdname PickChains
//' @export
// [[Rcpp::export]]
arma::uvec GetSubchains(unsigned int nchain, bool debug) {

  arma::uvec out;
  unsigned nsubchain = (unsigned int)std::ceil((double)nchain * R::runif(0.0, 1.0));

  // Step 2. Generate the chain index (0 to nchain - 1) and shuffle them
  arma::uvec chains  = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  arma::uvec rchains = arma::shuffle(chains);

  // From the shuffled chains, take the first nsubchain out of them.
  arma::uvec subchains = rchains.rows(0, nsubchain - 1);
  if (debug) {
    out = subchains;
  } else {
    out = arma::sort(subchains);
  }
  return out;   // Return shuffled subchains
}

//' @rdname PickChains
//' @export
// [[Rcpp::export]]
arma::uvec SelectEmigrants(unsigned int ngroup, unsigned int k) {

  arma::uvec groups, groups_, rgroups, emigrants;
  // randomly determine how many groups to emigrate
  unsigned int l = (unsigned int)std::ceil((double)(ngroup-1) *
    R::runif(0.0, 1.0));

  // remembering where the current group is in the original order
  groups = arma::linspace<arma::uvec>(0, ngroup - 1, ngroup);
  groups_ = groups;

  groups.shed_row(k); // remove current group
  rgroups = arma::shuffle(groups); // shuffle the rest of group members
  emigrants = rgroups.rows(0, l - 1); // pick the first to the l groups

  // glue the current group back to the selected emigrating groups
  // arma::uvec out = arma::join_cols(groups_.row(k), emigrants);
  // return arma::sort(out); // sorting them. E.g. 0, 2, 3, 5 within 0, 1,..., 6
  return arma::sort(emigrants);
}


void MutateDGMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, arma::vec dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift, unsigned int ngroup,
  double rp, unsigned int npda, double bw, unsigned int ncore,
  unsigned int gpuid, arma::uvec& rj) {

  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::mat theta  = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int start, end, k;
  unsigned int m = nchain / ngroup;
  arma::vec tmp(npar);
  arma::uvec subchains;

  for (size_t i = 0; i < ngroup; i++) {
    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t ii = 0; ii < m; ii++) {
      k = subchains(ii);

      for (size_t j = 0; j < npar; j++) {
        tmp(j) = theta(j, k) + R::runif(-rp, rp);
      }

      tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift, npda, bw, ncore,
        gpuid, false);

      tmp_logpos = tmp_lp + tmp_ll;
      cur_logpos = usell(k) + uselp(k);

      mh = std::exp(tmp_logpos - cur_logpos);
      if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
        theta.col(k) = tmp;
        uselp(k) = tmp_lp;
        usell(k) = tmp_ll;
        rj(k)    = 5;
      } else {
        rj(k)    = 6;
      }
    }
  }
  usetheta = arma::trans(theta);
}


void MutateDGMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube thetas,
  arma::vec pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, arma::vec ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  arma::vec sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, unsigned int ngroup,
  double rp, arma::uvec& rj) {

  unsigned int k, start, end;
  double tmp_logpos, cur_logpos, tmp_hll, tmp_hlp, mh;
  arma::mat useloc = arma::trans(usephi[0]); //npar x nchain
  arma::mat usesca = arma::trans(usephi[1]);
  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  unsigned int m = nchain / ngroup;
  arma::uvec subchains;
  arma::vec tmp_loc(npar), tmp_sca(npar);

  for (size_t i = 0; i < ngroup; i++) {
    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t j = 0; j < m; j++) {
      k = subchains(i);

      /* !!! This is critical !!!*/
      // Update usehll for new theta; nsub x npar x nchain
      usehll(k) = sumloghlike(thetas.slice(k), pdists, useloc.col(k),
        usesca.col(k), plower, pupper, plog);
      cur_logpos = usehlp(k) + usehll(k);

      for (size_t jj = 0; jj < npar; jj++) {
        tmp_loc(jj) = useloc(jj, k) + R::runif(-rp, rp);
        tmp_sca(jj) = usesca(jj, k) + R::runif(-rp, rp);
      }

      tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
                             sp2, llower, slower, lupper, supper, llog, slog);
      tmp_hll = sumloghlike(thetas.slice(k), pdists, tmp_loc, tmp_sca,
                            plower, pupper, plog);
      tmp_logpos = tmp_hlp + tmp_hll;
      // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      mh = std::exp(tmp_logpos - cur_logpos);
      if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
        useloc.col(k) = tmp_loc; // npar x nchain
        usesca.col(k) = tmp_sca;
        usehlp(k)  = tmp_hlp;
        usehll(k)  = tmp_hll;
        rj(k)    = 5;
      } else {
        rj(k)    = 6;
      }
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}

void MutateDGMCDataChains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  arma::vec dists,
  arma::mat usephi0, arma::mat usephi1,
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  unsigned int ngroup, unsigned int npda, double bw, unsigned int ncore,
  unsigned int gpuid, double rp, arma::uvec& rj) {

  unsigned int start, end, k;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::mat p1_ = arma::trans(usephi0); // npar x nchain
  arma::mat p2_ = arma::trans(usephi1); // npar x nchain
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int m = nchain / ngroup;
  arma::vec tmp;
  arma::uvec subchains;

  for (size_t i = 0; i < ngroup; i++) {
    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t j = 0; j < m; j++) {
      k = subchains(j);
      cur_logpos = usell(k) + uselp(k);

      for(size_t jj = 0; jj < npar; jj++) {
        tmp(j) = theta(j, k) + R::runif(-rp, rp);
      }

      tmp_lp = sumlogprior(tmp, dists, p1_.col(k), p2_.col(k), lower, upper,
                           islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
                          dim2, dim3, n1idx, ise, cellidx, RT, matchcell,
                          isr1, posdrift,
                          npda, bw, ncore, gpuid, false);
      tmp_logpos = tmp_lp + tmp_ll;
      mh = std::exp(tmp_logpos - cur_logpos);
      if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
        theta.col(k) = tmp;
        uselp(k) = tmp_lp;
        usell(k) = tmp_ll;
        rj(k) = 5;
      } else {
        rj(k) = 6;
      }
    }
  }
  usetheta = arma::trans(theta);
}


void CrossoverDGMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
  arma::vec pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, arma::vec ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  arma::vec sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, unsigned int ngroup,
  double rp, double gammaMult, arma::uvec& hrj) {

  unsigned int k0, k1, k2, start, end, m;
  double tmp_hll, tmp_hlp, tmp_logpos, cur_logpos, mh;
  arma::uvec dchains, subchains; // direction chains
  arma::vec hgamma, tmp_loc, tmp_sca;
  arma::mat useloc = arma::trans(usephi(0)); //npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int npar   = useloc.n_rows;
  unsigned int nchain = useloc.n_cols;
  hgamma = GetGamma(npar, gammaMult, true);
  m = nchain / ngroup;  // m == 7; nchain = 21
  arma::vec noise(npar);

  for (size_t i = 0; i < ngroup; i++) {  // G0, G1, G1

    start = i * m;                // 0,  7, 14
    end   = ((1 + i) * m - 1);    // 6, 13, 20
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t ii = 0; ii < m; ii++) { // 0, 1, ..., 6

      for(size_t iii = 0; iii < npar; iii++) noise(iii) = R::runif(-rp, rp);

      dchains = PickChains(ii, 2, subchains); // eg 7-13; 8 11
      k0 = subchains(ii);  // 7
      k1 = dchains(0);     // 8
      k2 = dchains(1);     // 11

      tmp_loc = useloc.col(k0) + (hgamma % (useloc.col(k1) - useloc.col(k2))) +
        noise;
      tmp_sca = usesca.col(k0) + (hgamma % (usesca.col(k1) - usesca.col(k2))) +
        noise;
      tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
        sp2, llower, slower, lupper, supper, llog, slog);
      tmp_hll = sumloghlike(theta.slice(k0), pdists, tmp_loc, tmp_sca,
        plower, pupper, plog);
      tmp_logpos = tmp_hlp + tmp_hll;

      /* !!! This is critical !!!*/
      // Update usehll for new theta; nsub x npar x nchain
      usehll(k0) = sumloghlike(theta.slice(k0), pdists, useloc.col(k0),
        usesca.col(k0), plower, pupper, plog);
      cur_logpos = usehlp(k0) + usehll(k0);

      mh = std::exp(tmp_logpos - cur_logpos);
      // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {
        useloc.col(k0) = tmp_loc; // npar x nchain
        usesca.col(k0) = tmp_sca;
        usehlp(k0)  = tmp_hlp;
        usehll(k0)  = tmp_hll;
        hrj(k0) = 3;
      } else {
        hrj(k0) = 4;
      }
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
}

void CrossoverDGMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  arma::vec dists,
  arma::mat usephi0, arma::mat usephi1,
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  unsigned int ngroup,
  unsigned int npda, double bw, unsigned int ncore, unsigned int gpuid,
  double rp, double gammamult, arma::uvec& rj) {

  unsigned int k0, k1, k2, start, end, m;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::uvec dchains, subchains;
  arma::vec gamma, tmp;

  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain matrix
  arma::mat phi0  = arma::trans(usephi0); // npar x nchain
  arma::mat phi1  = arma::trans(usephi1); // npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma  = GetGamma(npar, gammamult);
  m = nchain / ngroup;  // 7
  arma::vec noise(npar);

  for (size_t i = 0; i < ngroup; i++) {  // G0, G1, G2
    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t ii = 0; ii < m; ii++) {

      for(size_t iii = 0; iii < npar; iii++) noise(iii) = R::runif(-rp, rp);

      dchains = PickChains(ii, 2, subchains); // Pick 2 within 0, 1, 2, 3, e.g.
      k0 = subchains(ii);
      k1 = dchains(0);
      k2 = dchains(1);

      tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + noise;
      tmp_lp = sumlogprior(tmp, dists, phi0.col(k0), phi1.col(k0), lower, upper,
                           islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
                          dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
                          posdrift,
                          npda, bw, ncore, gpuid, false);
      tmp_logpos = tmp_ll + tmp_lp;
      // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;

      cur_logpos = uselp(k0) + usell(k0);
      mh = std::exp(tmp_logpos - cur_logpos);
      if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {
        theta.col(k0) = tmp;
        uselp(k0) = tmp_lp;
        usell(k0) = tmp_ll;
        rj(k0) = 3;
      } else {
        rj(k0) = 4;
      }
    }
  }
  usetheta = arma::trans(theta);
}


void CrossoverDMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
  arma::vec pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, arma::vec ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  arma::vec sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, double rp,
  double gammaMult, arma::uvec& rj) {

  unsigned int k0, k1, k2;
  // unsigned int rk0, rk1, rk2;
  double tmp_hll, tmp_hlp, tmp_logpos, cur_logpos, mh;
  arma::uvec chains, rchains, subchains, rsubchains;
  arma::vec hgamma, tmp_loc, tmp_sca;
  arma::mat useloc = arma::trans(usephi(0)); // npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int npar   = useloc.n_rows;
  unsigned int nchain = useloc.n_cols;
  hgamma = GetGamma(npar, gammaMult, true); // true = extras *2 as p1 and p2
  //chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));
  arma::vec noise(npar);
  // double a, d;

  for (size_t i = 0; i < nchain; i++) {
    subchains = PickChains(chains(i), 2, chains); // (b-a) * R::runif(1) + a;
    // rsubchains = PickChains(rchains(i), 2, rchains); // (b-a) * R::runif(1) + a;
    k0 = chains(i);
    k1 = subchains(0);
    k2 = subchains(1);

    // rk0 = rchains(i);
    // rk1 = rsubchains(0);
    // rk2 = rsubchains(1);

    for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);

    // b = usehll(k0);
    // c = usehlp(k0);
    /* CRITICAL. Update usehll for new theta; nsub x npar x nchain */
    // usehll(k0) = sumloghlike(theta.slice(k0), pdists, useloc.col(k0),
    //    usesca.col(k0), plower, pupper, plog);
    // usehlp(k0) = sumloghprior(useloc.col(k0), usesca.col(k0), ldists, sdists,
    //   lp1, sp1, lp2, sp2, llower, slower, lupper, supper, llog, slog);
    cur_logpos = usehlp(k0) + usehll(k0);
    // a = usehll(k0);
    // d = usehlp(k0);
    // if (b != a) Rcout << "Chain: " << k0 << " Hyper b and a " << b << " "
    //    << a << std::endl;
    // if (d != c) Rcout << "Chain: " << k0 << " Hyper c and d " << c << " "
    //                   << d << std::endl;

    // update mu and sigma
    tmp_loc = useloc.col(k0) + noise +
              ( hgamma % (useloc.col(k1) - useloc.col(k2)) );
    tmp_sca = usesca.col(k0) + noise +
              ( hgamma % (usesca.col(k1) - usesca.col(k2)) );

    // theta: nsub x npar x nchain == ps: nchain x nsub x npar
    tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
      sp2, llower, slower, lupper, supper, llog, slog);

    tmp_hll = sumloghlike(theta.slice(k0), pdists, tmp_loc, tmp_sca, plower,
      pupper, plog);
    tmp_logpos = tmp_hlp + tmp_hll;

    mh = std::exp(tmp_logpos - cur_logpos);
    if(std::isnan(mh)) Rcout << "mh " << mh << "\n";
    // if ( (R::runif(0, 1) < mh) )  {
    if ( (!std::isnan(mh)) && (R::runif(0, 1) < mh) )  {
      usephi(0).row(k0) = arma::trans(tmp_loc); // npar x nchain
      usephi(1).row(k0) = arma::trans(tmp_sca);

      // useloc.col(k0) = tmp_loc; // npar x nchain
      // usesca.col(k0) = tmp_sca;
      usehlp(k0) = tmp_hlp;
      usehll(k0) = tmp_hll;
      rj(k0) = 3;
    } else {
      rj(k0) = 4;
    }

  }
  // usephi(0) = arma::trans(useloc);
  // usephi(1) = arma::trans(usesca);
}

void CrossoverDMCHyperchains_blocked(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
  arma::vec pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, arma::vec ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  arma::vec sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, double rp,
  double gammaMult, arma::uvec& rj, unsigned int j) {

  double cur_logpos, tmp_logpos;
  unsigned int nchain = usephi(0).n_rows;
  // unsigned int npar   = usephi(0).n_cols;
  double hgamma = 1.19;
  arma::uvec subchains;

  // arma::uvec chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  arma::uvec chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));
  arma::mat useloc = arma::trans(usephi(0)); // npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  unsigned int k0;
  // unsigned int j = block;
  arma::vec useloc_, usesca_, diff1, diff2, tmp_loc, tmp_sca;

  for (size_t i = 0; i < nchain; i++) {
    k0 = chains(i);
    useloc_ = useloc.col(k0);
    usesca_ = usesca.col(k0);
    tmp_loc = useloc_;
    tmp_sca = usesca_;

    usehll(k0) = sumloghlike(theta.slice(chains(i)), pdists, useloc_, usesca_,
      plower, pupper, plog);
    cur_logpos = usehlp(k0) + usehll(k0);

    subchains = PickChains(k0, 2, chains); // (b-a) * R::runif(1) + a;
    diff1 = useloc.col(subchains(0)) - useloc.col(subchains(1));
    diff2 = usesca.col(subchains(0)) - usesca.col(subchains(1));

    tmp_loc(j) = useloc_(j) + R::runif(-rp, rp) + hgamma * diff1(j);
    tmp_sca(j) = usesca_(j) + R::runif(-rp, rp) + hgamma * diff2(j);

    double tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1,
        lp2, sp2, llower, slower, lupper, supper, llog, slog);
    double tmp_hll = sumloghlike(theta.slice(chains(i)), pdists, tmp_loc,
                                 tmp_sca, plower, pupper, plog);
    tmp_logpos = tmp_hlp + tmp_hll;

    double mh = std::exp(tmp_logpos - cur_logpos);
    // if(std::isnan(mh)) Rcout << "mh " << mh << "\n";
    if (!std::isnan(mh) && (R::runif(0, 1) < mh) )  {
        usephi(0).row(k0).col(j) = tmp_loc(j);
        usephi(1).row(k0).col(j) = tmp_sca(j);
        // usephi(0).row(i) = arma::trans(tmp_loc); // npar x nchain
        // usephi(1).row(i) = arma::trans(tmp_sca);
        usehlp(k0) = tmp_hlp;
        usehll(k0) = tmp_hll;
        rj(k0) = 3;
      } else { rj(k0) = 4; }
  }    // chain loop
}


void CrossoverDMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, arma::vec dists,
  arma::mat usephi0, arma::mat usephi1, arma::vec lower,
  arma::vec upper, arma::uvec islog, arma::vec allpar,
  std::vector<std::string> parnames, arma::ucube model, std::string type,
  std::vector<std::string> dim1, std::vector<std::string> dim2,
  std::vector<std::string> dim3, arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT, arma::uvec matchcell, arma::uvec isr1,
  bool posdrift,
  unsigned int npda, double bw, unsigned int ncore,
  unsigned int gpuid, double rp, double gammamult, arma::uvec& rj) {

  unsigned int k0, k1, k2;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::uvec chains, subchains;
  arma::vec gamma, tmp;

  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain
  arma::mat phi0  = arma::trans(usephi0);  // npar x nchain
  arma::mat phi1  = arma::trans(usephi1);  // npar x nchain

  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma  = GetGamma(npar, gammamult);
  // chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));
  arma::vec noise(npar);
  double b, a;

  for (size_t i = 0; i < nchain; i++) {

    for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);
    subchains = PickChains(chains(i), 2, chains);
    k0 = chains(i);
    k1 = subchains(0);
    k2 = subchains(1);

    b = uselp(k0);
    uselp(k0) = sumlogprior(theta.col(k0), dists, phi0.col(k0), phi1.col(k0),
      lower, upper, islog);
    usell(k0) = sumloglike(theta.col(k0), pnames, allpar, parnames, model,
      type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift,
      npda, bw, ncore, gpuid, false);
    cur_logpos = uselp(k0) + usell(k0);
    a = uselp(k0);
    if (b != a) Rcout << "Data b and a " << b << " " << a << std::endl;

    tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + noise;
    tmp_lp = sumlogprior(tmp, dists, phi0.col(k0), phi1.col(k0), lower, upper,
      islog);
    tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1, dim2,
      dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift,
      npda, bw, ncore, gpuid, false);
    tmp_logpos = tmp_ll + tmp_lp;

    // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
    mh = std::exp(tmp_logpos - cur_logpos);
    if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
       theta.col(k0) = tmp;
       uselp(k0)     = tmp_lp;
       usell(k0)     = tmp_ll;
       rj(k0) =3;
    } else {
       rj(k0) =4;
    }
  }
  usetheta = arma::trans(theta);
}

void CrossoverDMCDatachains_blocked(arma::mat& usetheta,    // nchain x npar
                            arma::vec& uselp, arma::vec& usell,
                            std::vector<std::string> pnames,
                            arma::vec dists,
                            arma::mat usephi0, arma::mat usephi1, arma::vec lower,
                            arma::vec upper, arma::uvec islog, arma::vec allpar,
                            std::vector<std::string> parnames, arma::ucube model,
                            std::string type,
                            std::vector<std::string> dim1,
                            std::vector<std::string> dim2,
                            std::vector<std::string> dim3, arma::umat n1idx,
                            arma::uvec ise,
                            arma::umat cellidx, arma::vec RT, arma::uvec matchcell,
                            arma::uvec isr1, bool posdrift,
                            unsigned int npda, double bw, unsigned int ncore,
                            unsigned int gpuid, double rp, double gammamult,
                            arma::uvec& rj, unsigned int block) {

  unsigned int k0;
  // , k1, k2, rk;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::uvec chains, subchains;
  arma::vec tmp;
  unsigned int j = block;

  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain
  arma::mat phi0  = arma::trans(usephi0);  // npar x nchain
  arma::mat phi1  = arma::trans(usephi1);  // npar x nchain

  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;

  double gamma = 2.38 / std::sqrt(2.0);
  // gamma  = GetGamma(npar, gammamult);
  // chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));
  arma::vec noise(npar);
  // double b, a;

  for (size_t i = 0; i < nchain; i++) {

    subchains = PickChains(chains(i), 2, chains);
    k0 = chains(i);
    // k1 = subchains(0);
    // k2 = subchains(1);

    // rk = rchains(i);

    arma::vec theta_ = theta.col(k0);
    arma::vec tmp_theta = theta_;

    // uselp(k0) = sumlogprior(theta.col(k0), dists, phi0.col(k0), phi1.col(k0),
    //       lower, upper, islog);

    // b = uselp(k0);
    // usell(k0) = sumloglike(theta.col(k0), pnames, allpar, parnames, model,
    //       type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift, npda,
    //       bw, ncore, gpuid, false);

    cur_logpos = uselp(k0) + usell(k0);
    // a = uselp(k0);
    // if (b != a) Rcout << "Data b and a " << b << " " << a << std::endl;

    arma::vec theta0 = theta.col(subchains(0));
    arma::vec theta1 = theta.col(subchains(1));

    tmp_theta(j) = theta_(j) + R::runif(-rp, rp) +
      gamma * ( theta0(j) - theta1(j) );

    // tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + noise;
    tmp_lp = sumlogprior(tmp_theta, dists, phi0.col(k0), phi1.col(k0), lower,
      upper, islog);
    tmp_ll = sumloglike(tmp_theta, pnames, allpar, parnames, model, type, dim1,
      dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift, npda, bw,
      ncore, gpuid, false);
    tmp_logpos = tmp_ll + tmp_lp;

    // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
    mh = std::exp(tmp_logpos - cur_logpos);
    if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
      theta.col(k0) = tmp_theta;
      uselp(k0)     = tmp_lp;
      usell(k0)     = tmp_ll;
      rj(k0) =3;
    } else {
      rj(k0) =4;
    }
  }
  usetheta = arma::trans(theta);
}

void CrossoverDGMCChains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,  // nchain x 1
  std::vector<std::string> pnames, arma::vec dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper,
  arma::uvec islog, arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  unsigned int ngroup,
  unsigned int npda, double bw, unsigned int ncore, unsigned int gpuid,
  double rp, double gammamult, arma::uvec& rj) {

  // Rcout << "CrossoverDGMCChains" << std::endl;
  unsigned int k0, k1, k2, start, end, m;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::uvec dchains, subchains; // direction chains
  arma::vec gamma, tmp;

  arma::mat theta  = arma::trans(usetheta); // npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma = GetGamma(npar, gammamult); // .5 * arma::randu(npar) + .5;
  m = nchain / ngroup;
  arma::vec noise(npar);

  for (size_t i = 0; i < ngroup; i++) {

    start = i * m;
    end   = ((1 + i) * m - 1);
    subchains = arma::linspace<arma::uvec>(start, end, m);

    for (size_t ii = 0; ii < m; ii++) {

      for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);

      dchains = PickChains(ii, 2, subchains); //  dchains are on 0, 1, ... m index
      k0 = subchains(ii);
      k1 = dchains(0);
      k2 = dchains(1);
      cur_logpos = usell(k0) + uselp(k0);
      tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + noise;

      tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
      tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
        dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift,
        npda, bw, ncore, gpuid, false);
      tmp_logpos = tmp_lp + tmp_ll;
      mh = std::exp(tmp_logpos - cur_logpos);
      if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
        theta.col(k0) = tmp;
        uselp(k0) = tmp_lp;
        usell(k0) = tmp_ll;
        rj(k0)    = 3;
      } else {
        rj(k0)    = 4;
      }
    }

  }
  usetheta = arma::trans(theta);
}


void CrossoverDMCChains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,  // nchain x 1
  std::vector<std::string> pnames,
  arma::vec dists, arma::vec p1, arma::vec p2,
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  double rp, double gammamult, bool force, unsigned int npda,
  double bw, unsigned int ncore, unsigned int gpuid, bool debug, arma::uvec& rj)
{
  // Ter Braak's crossover (2006), expect I randomly select a current chain
  unsigned int k0, k1, k2;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::uvec chains, subchains;
  arma::vec gamma, tmp;
  arma::mat theta = arma::trans(usetheta);   // theta: npar x nchain;
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  gamma = GetGamma(npar, gammamult); // ter Braak's step size
  if (debug) {
    // Rcout << "No shuffle\n";
    chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  } else {
    // Rcout << "Shuffling\n";
    chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));
  }

  arma::vec noise(npar);

  for (size_t i = 0; i < nchain; i++) {
    // Among the rest of the subpopulation xi \ {x^i_j}, sample two
    // chain/chromosomes, (xi_s and xi_t)
    subchains = PickChains(chains(i), 2, chains);
    k0 = chains(i);
    k1 = subchains(0);
    k2 = subchains(1);

    // noise = 2.0*rp*arma::randu<arma::vec>(npar) + rp;
    // Hu & Tsai and Liang & Wong use norm, but ter Braak does not
    // Uniform(-rp, rp); (b-a) * R::runif(1) + a;
    // tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + noise;

    for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);

    tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + noise;

    // PDA re-calculation
    if (force) {
      uselp(k0) = sumlogprior(theta.col(k0), dists, p1, p2, lower, upper, islog);
      usell(k0) = sumloglike (theta.col(k0), pnames, allpar, parnames,
        model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
        posdrift, npda, bw, ncore, gpuid, debug);
    }
    cur_logpos = usell(k0) + uselp(k0);

    tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
    tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1, dim2,
      dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift, npda, bw, ncore,
      gpuid, debug);
    tmp_logpos = tmp_lp + tmp_ll;
    // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
    mh = std::exp(tmp_logpos - cur_logpos);
    if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
      theta.col(k0) = tmp;
      uselp(k0) = tmp_lp;
      usell(k0) = tmp_ll;
      rj(k0) = 3;
    } else {
      rj(k0) = 4;
    }
  }
  usetheta = arma::trans(theta);
}

void CrossoverDMCChains_blocked(arma::mat& usetheta,    // nchain x npar
                                arma::vec& uselp, arma::vec& usell,  // nchain x 1
                                std::vector<std::string> pnames,
                                arma::vec dists, arma::vec p1, arma::vec p2,
                                arma::vec lower, arma::vec upper, arma::uvec islog,
                                arma::vec allpar, std::vector<std::string> parnames,
                                arma::ucube model, std::string type,
                                std::vector<std::string> dim1,
                                std::vector<std::string> dim2,
                                std::vector<std::string> dim3,
                                arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
                                arma::uvec matchcell, arma::uvec isr1, bool posdrift,
                                double rp, double gammamult, bool force, unsigned int npda,
                                double bw, unsigned int ncore, unsigned int gpuid, bool debug,
                                arma::uvec& rj, unsigned int j)
{
  // Ter Braak's crossover (2006), expect I randomly select a current chain
  unsigned int k0, k1, k2;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::uvec chains, subchains;
  arma::mat theta = arma::trans(usetheta);   // theta: npar x nchain;
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  double gamma = gammamult / std::sqrt(2.0*(double)npar);

  // chains = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));

  for (size_t i = 0; i < nchain; i++) {
    subchains = PickChains(chains(i), 2, chains);
    k0 = chains(i);
    k1 = subchains(0);
    k2 = subchains(1);

    // PDA re-calculation
    if (force) {
      uselp(k0) = sumlogprior(theta.col(k0), dists, p1, p2, lower, upper, islog);
      usell(k0) = sumloglike (theta.col(k0), pnames, allpar, parnames,
            model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1,
            posdrift, npda, bw, ncore, gpuid, debug);
    }
    cur_logpos = usell(k0) + uselp(k0);

    arma::vec tmp_theta = theta.col(k0);
    arma::vec theta0 = theta.col(k1);
    arma::vec theta1 = theta.col(k2);

    tmp_theta(j) = theta(j, k0) + R::runif(-rp, rp) +
      gamma * ( theta0(j) - theta1(j) );

    tmp_lp = sumlogprior(tmp_theta, dists, p1, p2, lower, upper, islog);
    tmp_ll = sumloglike(tmp_theta, pnames, allpar, parnames, model, type, dim1, dim2,
                        dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift, npda, bw, ncore,
                        gpuid, debug);
    tmp_logpos = tmp_lp + tmp_ll;
    // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
    mh = std::exp(tmp_logpos - cur_logpos);
    if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
      theta.col(k0) = tmp_theta;
      uselp(k0) = tmp_lp;
      usell(k0) = tmp_ll;
      rj(k0) = 3;
    } else {
      rj(k0) = 4;
    }
  }
  usetheta = arma::trans(theta);
}


void MigrateDMCHyperchains_old(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll,  arma::cube theta,
  arma::vec pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, arma::vec ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  arma::vec sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, double rp,
  arma::uvec& rj, bool debug) {

  arma::mat useloc = arma::trans(usephi(0)); // useloc: npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  arma::vec usehlp_ = usehlp;
  arma::vec usehll_ = usehll;

  unsigned int npar   = useloc.n_rows;
  unsigned int nchain = useloc.n_cols;
  arma::uvec subchains   = GetSubchains(nchain, debug);
  unsigned int nsubchain = subchains.n_elem;
  unsigned int k;

  double tmp_logpos, cur_logpos, mh;
  arma::mat tmp_loc(npar, nsubchain), tmp_sca(npar, nsubchain);
  arma::vec cur_hlp(nsubchain), cur_hll(nsubchain);
  arma::vec tmp_hlp(nsubchain), tmp_hll(nsubchain), noise(npar);

  for(size_t i = 0; i < nsubchain; i++) { // 0, 1, 2, 5, ...
    k = subchains(i);
    // Uniform(-rp, rp); (b-a) * R::runif(1) + a;
    for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);
    tmp_loc.col(i) = useloc.col(k) + noise; // proposal
    tmp_sca.col(i) = usesca.col(k) + noise; // proposal

    /* CRITICAL. Update usehll for new theta; nsub x npar x nchain */
    usehll_(k) = sumloghlike(theta.slice(k), pdists, useloc.col(k),
      usesca.col(k), plower, pupper, plog); // this changes too
    cur_hlp(i) = usehlp_(k);
    cur_hll(i) = usehll_(k);

    // theta: nsub x npar x nchain
    tmp_hlp(i) = sumloghprior(tmp_loc.col(i), tmp_sca.col(i), ldists,
      sdists, lp1, sp1, lp2, sp2, llower, slower, lupper, supper, llog, slog);
    tmp_hll(i) = sumloghlike(theta.slice(k), pdists, tmp_loc.col(i),
      tmp_sca.col(i), plower, pupper, plog);
  }

  tmp_logpos = tmp_hll(nsubchain - 1) + tmp_hlp(nsubchain - 1);
  // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
  cur_logpos = cur_hll(0) + cur_hlp(0);   // migrate to the first subchain
  mh = std::exp(tmp_logpos - cur_logpos);
  if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {
    useloc.col(subchains(0)) = tmp_loc.col(nsubchain - 1);
    usesca.col(subchains(0)) = tmp_sca.col(nsubchain - 1);
    usehlp_(subchains(0)) = tmp_hlp(nsubchain - 1);
    usehll_(subchains(0)) = tmp_hll(nsubchain - 1);
    rj(subchains(0)) = 1;
  } else {
    rj(subchains(0)) = 2;
  }

  // Continue migration, if nsubchain > 1
  if (nsubchain != 1) {
    for(size_t k = 1; k < (nsubchain - 1); k++) {
      tmp_logpos = tmp_hll(k) + tmp_hlp(k);
      // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      cur_logpos = cur_hll(k + 1) + cur_hlp(k + 1);
      mh = std::exp(tmp_logpos - cur_logpos);

      if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {
        useloc.col(subchains(k + 1))   = tmp_loc.col(k);
        usesca.col(subchains(k + 1))   = tmp_sca.col(k);
        usehlp_(subchains(k + 1)) = tmp_hlp(k);
        usehll_(subchains(k + 1)) = tmp_hll(k);
        rj(subchains(k + 1)) = 1;
      } else {
        rj(subchains(k + 1)) = 2;
      }
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
  usehlp    = usehlp_;
  usehll    = usehll_;
}

void MigrateDMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
  arma::vec pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, arma::vec ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  arma::vec sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog, double rp,
  arma::uvec& rj, bool debug) {

  arma::mat useloc = arma::trans(usephi(0));
  arma::mat usesca = arma::trans(usephi(1));
  arma::vec usehlp_ = usehlp;
  arma::vec usehll_ = usehll;

  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  arma::uvec subchains = GetSubchains(nchain, debug);
  unsigned int nsubchain = subchains.n_elem;
  arma::vec cur_hlp(nsubchain), cur_hll(nsubchain), cur_loc, cur_sca,
  tmp_loc(npar), tmp_sca(npar);

  unsigned int next_chain, k;
  double tmp_hlp, tmp_hll, tmp_logpos, cur_logpos, mh;

  for(size_t i = 0; i < nsubchain; i++) { // 0, 1, 2, 5, ...

    next_chain = ((i+1) == nsubchain) ? subchains(0) : subchains(i+1);
    k = subchains(i);

    // for (size_t j = 0; j < npar; j++) {
    //   tmp_loc(j) = useloc(j, k) + R::runif(-rp, rp);
    //   tmp_sca(j) = usesca(j, k) + R::runif(-rp, rp);
    // }


    for (size_t j = 0; j < npar; j++) {
      tmp_loc(j) = useloc(j, k) + R::rnorm(useloc(j, k), rp);
      tmp_sca(j) = usesca(j, k) + R::rnorm(usesca(j, k), rp);
    }

    tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists,
      lp1, sp1, lp2, sp2, llower, slower, lupper, supper, llog, slog);
    tmp_hll = sumloghlike(theta.slice(k), pdists, tmp_loc, tmp_sca, plower,
      pupper, plog);
    tmp_logpos = tmp_hlp + tmp_hll;
    // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;


    /* CRITICAL. Update usehll for new theta; nsub x npar x nchain */
    usehll_(next_chain) = sumloghlike(theta.slice(next_chain), pdists,
    useloc.col(next_chain), usesca.col(next_chain), plower, pupper, plog);
    cur_logpos = usehlp(next_chain) + usehll(next_chain);
    mh = std::exp(tmp_logpos - cur_logpos);

    if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {
      useloc.col(next_chain) = tmp_loc;
      usesca.col(next_chain) = tmp_sca;
      usehlp(next_chain) = tmp_hlp;
      usehll(next_chain) = tmp_hll;
      rj(next_chain) = 1;
    } else {
      rj(next_chain) = 2;
    }
  }
  usephi(0) = arma::trans(useloc);
  usephi(1) = arma::trans(usesca);
  usehlp    = usehlp_;
  usehll    = usehll_;
}

void MigrateDGMCHyperchains(arma::field<arma::mat>& usephi,
  arma::vec& usehlp, arma::vec& usehll,  arma::cube& thetas, // nsub x npar x nchain
  arma::vec pdists, arma::vec plower, arma::vec pupper,
  arma::uvec plog, arma::vec ldists, arma::vec lp1,
  arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
  arma::vec sdists, arma::vec sp1, arma::vec sp2,
  arma::vec slower, arma::vec supper, arma::uvec slog,
  unsigned int ngroup, double rp, arma::uvec& hrj) {
  // Rcout << "MigrateDGMC Hyperchains" << std::endl;
  unsigned int start, end, next_chain, k, l;
  double tmp_hlp, tmp_hll, tmp_logpos, mh, cur_logpos;
  // cur_hlp, cur_hll, ;

  arma::mat useloc = arma::trans(usephi(0)); // nchain x npar to npar x nchain
  arma::mat usesca = arma::trans(usephi(1));

  arma::mat useloc_ = arma::trans(usephi(0)); // nchain x npar to npar x nchain
  arma::mat usesca_ = arma::trans(usephi(1));


  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  unsigned int m = nchain / ngroup;    // number of chains in each group
  arma::uvec subchains, subgroups;
  arma::vec tmp_loc(npar), tmp_sca(npar);


  for (size_t i = 0; i < ngroup; i++) {
    subgroups = SelectEmigrants(ngroup, i);
    l = subgroups.n_elem;

    for (size_t j = 0; j < l; j++) {
      start = subgroups(j) * m;              // 0, 10, 15, 25
      end   = ((1 + subgroups(j)) * m - 1);  // 4, 14, 19, 29
      subchains = arma::linspace<arma::uvec>(start, end, m);

      for (size_t jj = 0; jj < m; jj++) {
        next_chain = (j == (l-1)) ? subgroups(0)*m + jj : subgroups(j+1) * m + jj;
        k = subchains(jj);

        for (size_t jjj = 0; jjj < npar; jjj++) {
          tmp_loc(jjj) = useloc(jjj, k) + R::rnorm(useloc(jjj, k), rp);
          tmp_sca(jjj) = usesca(jjj, k) + R::rnorm(usesca(jjj, k), rp);
        }

        tmp_hlp = sumloghprior(tmp_loc, tmp_sca, ldists, sdists, lp1, sp1, lp2,
          sp2, llower, slower, lupper, supper, llog, slog);
        tmp_hll = sumloghlike(thetas.slice(k), pdists, tmp_loc, tmp_sca, plower,
          pupper, plog);
        tmp_logpos = tmp_hlp + tmp_hll;

        /* !!! This is critical !!!*/
        // Update usehll for new theta; nsub x npar x nchain
        usehll(next_chain) = sumloghlike(thetas.slice(next_chain), pdists,
          useloc.col(next_chain),
          usesca.col(next_chain), plower, pupper, plog);
        cur_logpos = usehlp(next_chain) + usehll(next_chain);

        // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
        mh = std::exp(tmp_logpos - cur_logpos);
        if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {

          useloc_.col(next_chain) = tmp_loc;
          usesca_.col(next_chain) = tmp_sca;
          usehlp(next_chain)    = tmp_hlp;
          usehll(next_chain)    = tmp_hll;
          hrj(next_chain) = 1;
        } else {
          hrj(next_chain) = 2;
        }
      }
    }
  }

  usephi(0) = arma::trans(useloc_);
  usephi(1) = arma::trans(usesca_);
}

void MigrateDGMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  arma::vec dists,
  arma::mat usephi0, arma::mat usephi1, // nchain x npar
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  unsigned int ngroup, double rp, unsigned int npda,
  double bw, unsigned int ncore, unsigned int gpuid, arma::uvec& rj)
{
  // Rcout << "MigrateDGMC Data chains" << std::endl;
  unsigned int start, end, next, k, l;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  // cur_lp, cur_ll
  arma::mat phi0 = arma::trans(usephi0);
  arma::mat phi1 = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  arma::mat theta_ = arma::trans(usetheta); // theta: npar x nchain;

  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int m      = nchain / ngroup;
  arma::uvec subchains, subgroups;
  arma::vec tmp(npar);

  for (size_t i = 0; i < ngroup; i++) {
    subgroups = SelectEmigrants(ngroup, i);
    l = subgroups.n_elem;

    for (size_t j = 0; j < l; j++) {
      start = subgroups(j) * m;              // 0, 10, 15, 25
      end   = ((1 + subgroups(j)) * m - 1);  // 4, 14, 19, 29
      subchains = arma::linspace<arma::uvec>(start, end, m);
      for (size_t jj = 0; jj < m; jj++) {
        next  = (j == (l-1)) ? subgroups(0)*m + jj : subgroups(j + 1)*m + jj;
        k = subchains(jj);

        for (size_t jjj = 0; jjj < npar; jjj++) {
          tmp(jjj) = theta(jjj, k) + R::rnorm(theta(jjj, k), rp);
        }

        tmp_lp = sumlogprior(tmp, dists, phi0.col(k), phi1.col(k), lower, upper,
          islog);
        tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift,
          npda, bw, ncore, gpuid, false);
        tmp_logpos = tmp_lp + tmp_ll;
        cur_logpos = uselp(next) + usell(next);
        // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
        mh = std::exp(tmp_logpos - cur_logpos);
        if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {
          theta_.col(next) = tmp;
          uselp(next) = tmp_lp;
          usell(next) = tmp_ll;
          rj(next) = 1;
        } else {
          rj(next) = 2;
        }
      }
    }
  }
  usetheta = arma::trans(theta_);
}


void MigrateDMCDatachains(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  arma::vec dists,
  arma::mat usephi0, arma::mat usephi1, // nchain x npar
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  unsigned int nsim, double bw, unsigned int ncore, unsigned int gpuid,
  double rp, double gammamult, arma::uvec& rj, bool debug) {
  arma::mat phi0  = arma::trans(usephi0);
  arma::mat phi1  = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  unsigned int npar    = theta.n_rows;
  unsigned int nchain  = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain, debug); // eg. 0, 1, 3, 4, 8
  unsigned int nsubchain = subchains.n_elem;

  unsigned int next_chain, k;
  arma::vec theta_cur, theta_star(npar);
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;

  for (size_t i = 0; i < nsubchain; i++) {

    next_chain = ((i+1) == nsubchain) ? subchains(0) : subchains(i+1);
    k = subchains(i);

    theta_cur  = theta.col(next_chain);

    // for (size_t jj = 0; jj < npar; jj++) {
    //   theta_star(jj) = theta(jj, k) + R::runif(-rp, rp);
    // }

    for (size_t jj = 0; jj < npar; jj++) {
      theta_star(jj) = theta(jj, k) + + R::rnorm(theta(jj, k), rp);
    }

    tmp_lp = sumlogprior(theta_star, dists, phi0.col(k), phi1.col(k), lower,
      upper, islog);
    tmp_ll = sumloglike(theta_star, pnames, allpar, parnames, model, type, dim1,
      dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift, nsim, bw, ncore,
      gpuid, false);

    // if (std::isinf(tmp_lp) && tmp_lp > 0.0) tmp_lp = 1e-10;
    // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
    tmp_logpos = tmp_lp + tmp_ll;
    cur_logpos = uselp(next_chain) + usell(next_chain);

    mh = std::exp(tmp_logpos - cur_logpos);
    if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {
      theta.col(next_chain) = theta_star;
      uselp(next_chain) = tmp_lp;
      usell(next_chain) = tmp_ll;
      rj(next_chain) = 1;
    } else {
      rj(next_chain) = 2;
    }
  }
  usetheta = arma::trans(theta);
}

void MigrateDMCDatachains_old(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames,
  arma::vec dists,
  arma::mat usephi0, arma::mat usephi1, // nchain x npar
  arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames,
  arma::ucube model, std::string type,
  std::vector<std::string> dim1,
  std::vector<std::string> dim2,
  std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise,
  arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  unsigned int nsim, double bw, unsigned int ncore, unsigned int gpuid,
  double rp, double gammamult, arma::uvec& rj, bool debug) {

  arma::mat phi0  = arma::trans(usephi0);
  arma::mat phi1  = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  arma::vec uselp_ = uselp;
  arma::vec usell_ = usell;
  unsigned int npar    = theta.n_rows;
  unsigned int nchain  = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain, debug); // eg. 0, 1, 3, 4, 8
  unsigned int nsubchain = subchains.n_elem;
  unsigned int k;

  arma::vec cur_lp(nsubchain), cur_ll(nsubchain), tmp_lp(nsubchain),
            tmp_ll(nsubchain), noise(npar);
  double tmp_logpos, cur_logpos, mh;
  arma::mat tmp(npar, nsubchain);

  for (size_t i = 0; i < nsubchain; i++) {
    k = subchains(i);

    for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);
    tmp.col(i) = theta.col(k) + noise; // proposal
    cur_lp(i) = uselp_(k);
    cur_ll(i) = usell_(k);
    tmp_lp(i) = sumlogprior(tmp.col(i), dists, phi0.col(k), phi1.col(k), lower,
      upper, islog);
    tmp_ll(i) = sumloglike(tmp.col(i), pnames, allpar, parnames, model, type,
      dim1, dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift,
      nsim, bw, ncore, gpuid, false);
  }

  tmp_logpos = tmp_ll(nsubchain - 1) + tmp_lp(nsubchain - 1);
  // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
  cur_logpos = cur_ll(0) + cur_lp(0);   // migrate to the first subchain

  mh = std::exp(tmp_logpos - cur_logpos);
  if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
    theta.col(subchains(0)) = tmp.col(nsubchain - 1);
    uselp_(subchains(0)) = tmp_lp(nsubchain - 1);
    usell_(subchains(0)) = tmp_ll(nsubchain - 1);
    rj(subchains(0)) = 1;
  } else {
    rj(subchains(0)) = 2;
  }

  if (nsubchain != 1) {
    for(size_t k = 1; k < (nsubchain - 1); k++) {
      tmp_logpos = tmp_ll(k) + tmp_lp(k);
      // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      cur_logpos = cur_ll(k + 1) + cur_lp(k + 1);
      mh = std::exp(tmp_logpos - cur_logpos);

      if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
        theta.col(subchains(k + 1)) = tmp.col(k);
        uselp_(subchains(k + 1)) = tmp_lp(k);
        usell_(subchains(k + 1)) = tmp_ll(k);
        rj(subchains(k + 1)) = 1;
      } else {
        rj(subchains(k + 1)) = 2;
      }
    }
  }
  usetheta = arma::trans(theta);
  uselp    = uselp;
  usell    = usell;
}

void MigrateDGMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, arma::vec dists,
  arma::vec p1, arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  unsigned int ngroup, double rp, unsigned int npda, double bw,
  unsigned int ncore, unsigned int gpuid, arma::uvec& rj) {

  unsigned int start, end, next_chain, k, l;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  // , cur_lp, cur_ll;
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  unsigned int m = nchain / ngroup; // m = 30 / 6 = 5
  arma::uvec subchains, subgroups;
  arma::vec tmp(npar);

  for (size_t i = 0; i < ngroup; i++) {
    subgroups = SelectEmigrants(ngroup, i); // 6 groups: 0, 2, 3, 5
    l = subgroups.n_elem; // l = 4

    for (size_t j = 0; j < l; j++) { // j = 0, 1, 2, 4
      start = subgroups(j) * m;              // 0, 10, 15, 25
      end   = ((1 + subgroups(j)) * m - 1);  // 4, 14, 19, 29
      subchains = arma::linspace<arma::uvec>(start, end, m);

      for (size_t jj = 0; jj < m; jj++) { // jj = 0, 1, 2, 3, 4
        // corresponding index in the next group
        next_chain = (j == (l - 1)) ? subgroups(0)*m + jj : subgroups(j+1)*m + jj;

        k = subchains(jj);
        for (size_t jjj = 0; jjj < npar; jjj++) {
          tmp(jjj) = theta(jjj, k) + R::rnorm(theta(jjj, k), rp);
        }

        tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, islog);
        tmp_ll = sumloglike(tmp, pnames, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT,matchcell, isr1, posdrift,
          npda, bw, ncore, gpuid, false);
        tmp_logpos = tmp_lp + tmp_ll;

        cur_logpos = uselp(next_chain) + usell(next_chain);
        mh = std::exp(tmp_logpos - cur_logpos);
        if ( !std::isnan(mh) && (R::runif(0, 1) < mh ) )  {
          theta.col(next_chain) = tmp;
          uselp(next_chain) = tmp_lp;
          usell(next_chain) = tmp_ll;
          rj(next_chain)    = 1;
        } else {
          rj(next_chain)    = 2;
        }
      }

    }

  }
  usetheta = arma::trans(theta);
}

void MigrateDMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, arma::vec dists, arma::vec p1,
  arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift, double rp, double gammamult,
  bool force, unsigned int nsim, double bw, unsigned int ncore,
  unsigned int gpuid, arma::uvec& rj, bool debug) {

  double tmp_lp, tmp_ll, tmp_logpos,cur_logpos, mh;
  arma::vec theta_cur;

  arma::mat theta    = arma::trans(usetheta);    // theta: npar x nchain
  arma::mat theta_in = theta;

  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain, debug); // eg, 0, 1, 3, 4, 8;
  unsigned int nsubchain = subchains.n_elem;   // could be just 1 chain
  unsigned int next_chain, k;
  arma::vec theta_star(npar);

  for(size_t i = 0; i < nsubchain; i++) {
    next_chain = ((i+1) == nsubchain) ? subchains(0) : subchains(i+1);

    k = subchains(i);
    theta_cur = theta.col(next_chain);

    for (size_t j = 0; j < npar; j++) {
      theta_star(j) = theta(j, k) + R::rnorm(theta(j, k), rp);
    }

    tmp_lp = sumlogprior(theta_star, dists, p1, p2, lower, upper, islog);
    tmp_ll = sumloglike(theta_star, pnames, allpar, parnames, model, type, dim1,
      dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift,
      nsim, bw, ncore, gpuid, false);
    tmp_logpos = tmp_lp + tmp_ll;
    // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
    cur_logpos = uselp(next_chain) + usell(next_chain);

    mh = std::exp(tmp_logpos - cur_logpos);
    if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
      theta_in.col(next_chain) = theta_star;
      uselp(next_chain) = tmp_lp;
      usell(next_chain) = tmp_ll;
      rj(next_chain) = 1;
    } else {
      rj(next_chain) = 2;
    }
  }
  usetheta = arma::trans(theta_in);
}

// void MigrateDMCChains(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
//   std::vector<std::string> pnames, arma::vec dists, arma::vec p1,
//   arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
//   arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
//   std::string type, std::vector<std::string> dim1,
//   std::vector<std::string> dim2, std::vector<std::string> dim3,
//   arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
//   arma::uvec matchcell, arma::uvec isr1, double rp, double gammamult,
//   bool force, unsigned int nsim, double bw, unsigned int ncore,
//   unsigned int gpuid, arma::uvec& rj) {
//
//   double tmp_lp, tmp_ll, tmp_logpos, cur_lp, cur_ll, cur_logpos, mh;
//   arma::vec theta_cur;
//
//   arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
//   unsigned int npar   = theta.n_rows;
//   unsigned int nchain = theta.n_cols;
//   arma::uvec subchains = GetSubchains(nchain); // eg, 0, 1, 3, 4, 8;
//   unsigned int nsubchain = subchains.n_elem;   // could be just 1 chain
//   unsigned int next_chain, k;
//   arma::vec theta_star(npar);
//
//   for(size_t i = 0; i < nsubchain; i++) {
//     next_chain = ((i+1) == nsubchain) ? subchains(0) : subchains(i+1);
//
//     k = subchains(i);
//     theta_cur = theta.col(next_chain);
//
//     for (size_t j = 0; j < npar; j++) {
//       theta_star(j) = theta(j, k) + R::rnorm(theta(j, k), rp);
//     }
//
//     tmp_lp = sumlogprior(theta_star, dists, p1, p2, lower, upper, islog);
//     tmp_ll = sumloglike(theta_star, pnames, allpar, parnames, model, type, dim1,
//       dim2, dim3, n1idx, ise, cellidx, RT, matchcell, isr1, posdrift, nsim, bw, ncore,
//       gpuid, false);
//     tmp_logpos = tmp_lp + tmp_ll;
//     // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
//     cur_logpos = uselp(next_chain) + usell(next_chain);
//
//     mh = std::exp(tmp_logpos - cur_logpos);
//     if ( !std::isnan(mh) && (R::runif(0, 1) < mh) )  {
//       theta.col(next_chain) = theta_star;
//       uselp(next_chain) = tmp_lp;
//       usell(next_chain) = tmp_ll;
//       rj(next_chain) = 1;
//     } else {
//       rj(next_chain) = 2;
//     }
//   }
//   usetheta = arma::trans(theta);
// }

void MigrateDMCChains_old(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, arma::vec dists, arma::vec p1,
  arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::vec RT,
  arma::uvec matchcell, arma::uvec isr1, bool posdrift,
  double rp, double gammamult,
  bool force, unsigned int nsim, double bw, unsigned int ncore,
  unsigned int gpuid, arma::uvec& rj, bool debug) {

  // Rcout << "Migrate old" << std::endl;
  double tmp_logpos, cur_logpos, mh;
  arma::mat theta = arma::trans(usetheta);    // theta: npar x nchain
  arma::vec uselp_ = uselp;
  arma::vec usell_ = usell;

  unsigned int npar   = theta.n_rows;
  unsigned int nchain = theta.n_cols;
  // arma::uvec subchains = GetSubchains(nchain); // eg, 0, 1, 3, 4, 8; could be just 1 chain
  arma::uvec subchains = GetSubchains(nchain, debug);
  unsigned int nsubchain = subchains.n_elem;
  arma::mat tmp(npar, nsubchain);
  arma::vec cur_lp(nsubchain), cur_ll(nsubchain), noise(npar);
  arma::vec tmp_lp(nsubchain), tmp_ll(nsubchain);

  for(size_t i = 0; i < nsubchain; i++) {
    for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);

    tmp.col(i) = theta.col(subchains(i)) + noise; // proposal
    cur_lp(i) = uselp_(subchains(i));
    cur_ll(i) = usell_(subchains(i));

    tmp_lp(i) = sumlogprior(tmp.col(i), dists, p1, p2, lower, upper, islog);
    tmp_ll(i) = sumloglike(tmp.col(i), pnames, allpar, parnames,
      model, type, dim1, dim2, dim3, n1idx, ise, cellidx, RT,
      matchcell, isr1, posdrift, nsim, bw, ncore, gpuid, false);
  }

  // Conduct topological migration within subchains;
  // starting from the last subchain;
  // each individual chain/chromosome in subchains is treated as a subgroup
  tmp_logpos = tmp_ll(nsubchain - 1) + tmp_lp(nsubchain - 1);
  cur_logpos = cur_ll(0) + cur_lp(0);   // migrate to the first subchain
  mh = std::exp(tmp_logpos - cur_logpos);
  // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
  if (!std::isnan(tmp_logpos) && R::runif(0, 1) < mh) {
    theta.col(subchains(0)) = tmp.col(nsubchain - 1);
    uselp_(subchains(0)) = tmp_lp(nsubchain - 1);
    usell_(subchains(0)) = tmp_ll(nsubchain - 1);
    rj(subchains(0)) = 1;
  } else {
    rj(subchains(0)) = 2;
  }

  // Continue migration, if nsubchain > 1
  if (nsubchain != 1) {
    for(size_t k = 1; k < (nsubchain - 1); k++) {
      tmp_logpos = tmp_ll(k) + tmp_lp(k);
      cur_logpos = cur_ll(k + 1) + cur_lp(k + 1);

      mh = std::exp(tmp_logpos - cur_logpos);
      // if (std::isnan(tmp_logpos)) tmp_logpos = -INFINITY;
      if (!std::isnan(tmp_logpos) && R::runif(0, 1) < mh) {
        theta.col(subchains(k + 1)) = tmp.col(k);
        uselp_(subchains(k + 1))     = tmp_lp(k);
        usell_(subchains(k + 1))     = tmp_ll(k);
        rj(subchains(k + 1)) = 1;
      } else {
        rj(subchains(k + 1)) = 2;
      }
    }
  }
  usetheta = arma::trans(theta);
  uselp    = uselp;
  usell    = usell;
}

// [[Rcpp::export]]
List run_dgmc(List samples, arma::uvec force, unsigned int report, double pm,
  double qm, double gammamult, unsigned int ncore, unsigned int ngroup) {

  List samples_in(clone(samples)); // so R' original samples stays
  CheckPnames(samples_in);

  List data          = samples_in["data"];
  double rp          = samples_in["rp"];
  arma::ucube model  = data.attr("model");
  unsigned int npda  = data.attr("n.pda");
  double bw          = data.attr("bw");
  unsigned int gpuid = data.attr("gpuid");
  arma::vec RT       = data["RT"];
  unsigned int nmc   = samples_in["nmc"];
  unsigned int start_R = samples_in["start"];
  unsigned int thin = samples_in["thin"];
  unsigned int store_i = start_R - 1;
  unsigned int nsamp = 1 + (nmc - start_R) * thin;
  arma::cube theta = samples_in["theta"] ; // nchain x npar x nmc
  arma::mat lp     = samples_in["summed_log_prior"] ; // nmc x nchain
  arma::mat ll     = samples_in["log_likelihoods"] ;  // nmc x nchain
  arma::mat usetheta = theta.slice(store_i) ; // nchain x npar
  arma::vec uselp = arma::trans(lp.row(store_i)) ; // nchain x 1
  arma::vec usell = arma::trans(ll.row(store_i)) ;  // nchain x 1

  // extract data-model arguments
  NumericVector modelAttr = data.attr("model");
  std::string type = modelAttr.attr("type");
  arma::vec allpar = modelAttr.attr("all.par");
  List modelDim    = modelAttr.attr("dimnames");
  arma::umat n1idx  = modelAttr.attr("n1.order");
  arma::uvec mc    = modelAttr.attr("match.cell");
  bool posdrift    = modelAttr.attr("posdrift");
  arma::uvec ise   = data.attr("cell.empty") ;
  arma::umat cellidx = cellIdx2Mat(data);
  std::vector<std::string> parnames = modelAttr.attr("par.names");
  std::vector<std::string> dim1 = modelDim[0] ; // row
  std::vector<std::string> dim2 = modelDim[1] ; // col; para
  std::vector<std::string> dim3 = modelDim[2] ; // slice; r1, r2
  arma::uvec isr1 = GetIsR1(modelAttr, type);
  NumericVector pvec = modelAttr.attr("p.vector");
  unsigned int npar = pvec.size();
  std::vector<std::string> pnames = pvec.attr("names");

  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar), pdists(npar);
  arma::uvec plog(npar);
  List pprior = samples_in["p.prior"];
  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);

  arma::umat rejection_rate;
  arma::uvec rj;
  InitializeOneSubject(samples_in, rejection_rate); // nchain x nsamp

  for (size_t i = 1; i < nsamp; i++) {

      rj = rejection_rate.col(i);

      if (R::runif(0.0, 1.0) < pm) {
        MigrateDGMCChains(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
          plower, pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3,
          n1idx, ise, cellidx, RT, mc, isr1, posdrift, ngroup, rp, npda, bw,
          ncore, gpuid, rj);
      } else if (R::runif(0.0, 1.0) <= qm) {
        MutateDGMCChains(usetheta, uselp, usell, pnames, pdists, pp1,
          pp2, plower, pupper, plog, allpar, parnames, model, type, dim1, dim2,
          dim3, n1idx, ise, cellidx, RT, mc, isr1, posdrift, ngroup, rp, npda,
          bw, ncore, gpuid, rj);
      } else {
        CrossoverDGMCChains(usetheta, uselp, usell, pnames, pdists,
          pp1, pp2, plower, pupper, plog, allpar, parnames, model, type, dim1,
          dim2, dim3, n1idx, ise, cellidx, RT, mc, isr1, posdrift, ngroup,
          npda, bw, ncore, gpuid, rp, gammamult, rj);
      }

      rejection_rate.col(i) = rj;
      if (i % thin == 0) {
         store_i++;
         if ((store_i + 1) % report == 0) Rcout << store_i + 1 << " ";
         lp.row(store_i) = uselp.t();  // nmc x nchain
         ll.row(store_i) = usell.t();  // nmc x nchain
         theta.slice(store_i) = usetheta ;
      }
  }
  /* ------------------Output ------------------------------------- */
  samples_in["summed_log_prior"] = lp;     // nmc x nchain
  samples_in["log_likelihoods"]  = ll;     // nmc x nchain
  samples_in["theta"]            = theta;  // nchain x npar x nmc
  samples_in["rejection_rate"]   = rejection_rate;
  Rcout << std::endl;
  return samples_in;
}

// [[Rcpp::export]]
List run_dmc(List samples, arma::uvec force, unsigned int report, double pm,
             double pm0, double gammamult, unsigned int ncore, bool slice) {

  List samples_in(clone(samples)); // so R' original samples stays
  CheckPnames(samples_in);

  List data          = samples_in["data"];
  double rp          = samples_in["rp"];
  arma::ucube model  = data.attr("model");
  unsigned int npda  = data.attr("n.pda");
  double bw          = data.attr("bw");
  bool debug         = data.attr("debug");
  unsigned int gpuid = data.attr("gpuid");
  arma::vec RT       = data["RT"];
  unsigned int nmc   = samples_in["nmc"];
  unsigned int start_R = samples_in["start"];
  unsigned int thin = samples_in["thin"];
  unsigned int store_i = start_R - 1;
  unsigned int nsamp = 1 + (nmc - start_R) * thin;
  unsigned int nchain  = samples_in["n.chains"];

  arma::cube theta = samples_in["theta"];    // nchain x npar x nmc
  arma::mat lp = samples_in["summed_log_prior"];  // nmc x nchain
  arma::mat ll = samples_in["log_likelihoods"];
  arma::mat usetheta = theta.slice(store_i);   // nchain x npar
  arma::vec uselp = arma::trans(lp.row(store_i)); // nchains x 1
  arma::vec usell = arma::trans(ll.row(store_i));

  // extract data-model options
  NumericVector modelAttr = data.attr("model");
  std::string type = modelAttr.attr("type");
  arma::vec allpar = modelAttr.attr("all.par");
  List modelDim    = modelAttr.attr("dimnames");
  arma::umat n1idx = modelAttr.attr("n1.order");
  arma::uvec mc    = modelAttr.attr("match.cell");
  bool posdrift    = modelAttr.attr("posdrift");
  arma::uvec ise   = data.attr("cell.empty") ;
  arma::umat cellidx = cellIdx2Mat(data);
  std::vector<std::string> parnames = modelAttr.attr("par.names");
  std::vector<std::string> dim1 = modelDim[0] ; // row
  std::vector<std::string> dim2 = modelDim[1] ; // col; para
  std::vector<std::string> dim3 = modelDim[2] ; // slice; r1, r2
  arma::uvec isr1 = GetIsR1(modelAttr, type);

  NumericVector pvec = modelAttr.attr("p.vector");
  std::vector<std::string> pnames = pvec.attr("names");
  unsigned int npar = pnames.size();

  List pprior = samples_in["p.prior"];
  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar), pdists(npar);
  arma::uvec plog(npar);
  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);

  arma::umat rejection_rate;
  arma::uvec rj;
  InitializeOneSubject(samples_in, rejection_rate); // nchain x nmc

  for (size_t i = 1; i < nsamp; i++) {// From 1 not 0

    rj = arma::zeros<arma::uvec>(nchain);

    if (R::runif(0.0, 1.0) < pm) {

	    MigrateDMCChains(usetheta, uselp, usell, pnames, pdists, pp1, pp2, plower,
            pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3, n1idx,
            ise, cellidx, RT, mc, isr1, posdrift, rp, gammamult, force(i), npda, bw, ncore,
            gpuid, rj, debug);

    } else if (R::runif(0.0, 1.0) < pm0) {
      // Rcout <<"MigrateDMCChains_old\n";
      MigrateDMCChains_old(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
                           plower, pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3,
                           n1idx, ise, cellidx, RT, mc, isr1, posdrift, rp, gammamult, force(i), npda, bw,
                           ncore, gpuid, rj, debug);
    } else {

      if (slice) {
        // Rcout << "Use slicing";
        for (size_t l = 0; l < npar; l++) {
          CrossoverDMCChains_blocked(usetheta, uselp, usell, pnames, pdists,
                pp1, pp2, plower, pupper, plog, allpar, parnames, model, type,
                dim1, dim2, dim3, n1idx, ise, cellidx, RT, mc, isr1, posdrift,
                rp, gammamult, force(i), npda, bw, ncore, gpuid, debug, rj, l);
        }
      } else {
        // Rcout << "Non blocked";
        CrossoverDMCChains(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
          plower, pupper, plog, allpar, parnames, model, type, dim1, dim2, dim3,
          n1idx, ise, cellidx, RT, mc, isr1, posdrift, rp, gammamult, force(i),
          npda, bw, ncore, gpuid, debug, rj);
      }

    }


    if (i % thin == 0) {
      store_i++;
      if ((store_i + 1) % report == 0) Rcout << store_i + 1 << " ";
      lp.row(store_i) = uselp.t();   // nmc x nchain
      ll.row(store_i) = usell.t();   // nmc x nchain
      theta.slice(store_i)  = usetheta;
      rejection_rate.col(store_i) = rj;
    }
  }

  /* ------------------Output ------------------------------------- */
  samples_in["summed_log_prior"] = lp;     // nmc x nchain
  samples_in["log_likelihoods"]  = ll;     // nmc x nchain
  samples_in["theta"]            = theta;  // nchain x npar x nmc
  samples_in["rejection_rate"]   = rejection_rate;
  Rcout << std::endl;
  return samples_in;
}

// [[Rcpp::export]]
List run_hyper_dmc(List samples, unsigned int report, double pm, double pm0, double hpm, double hpm0,
  double gammamult, unsigned int ncore, bool debug) {

  List samples_in(clone(samples));
  CheckHyperPnames(samples_in);

  List hyper = samples_in.attr("hyper");
  unsigned int npar = hyper["n.pars"];
  unsigned int nchain = hyper["n.chains"];
  unsigned int nmc  = hyper["nmc"];
  unsigned int thin = hyper["thin"];
  unsigned int start= hyper["start"]; // start_R == 1;
  double rp = hyper["rp"];            // rp is defined in initialise
  unsigned int nsub = samples.size();
  unsigned int start_C = start - 1;   // start_C == 0;
  unsigned int store_i = start_C;    // store_i == 0;
  unsigned int nsamp = 1 + (nmc - start) * thin;

  /* data_hyper/thetas (nsub x npar x nchain) == cps (nchain x nsub x npar) */
  std::vector<std::string> pnames = hyper("p.names");
  List phi      = hyper["phi"];
  arma::mat hlp = hyper["h_summed_log_prior"]; // nmc x nchain
  arma::mat hll = hyper["h_log_likelihoods"];
  arma::cube location = phi[0]; // nchain x npar x nmc
  arma::cube scale    = phi[1];

  arma::field<arma::vec> blocks(npar);  // reserved blocks variable
  arma::field<arma::mat> usephi(2);
  usephi(0) = location.slice(start_C); // nchain x npar
  usephi(1) = scale.slice(start_C);
  arma::vec usehlp = arma::trans(hlp.row(start_C)); // nchain
  arma::vec usehll = arma::trans(hll.row(start_C));

  List subject0 = samples_in[0];
  List pprior   = subject0["p.prior"];
  List ppprior  = hyper["pp.prior"];     /* Extract pprior & ppprior */
  List lprior   = ppprior[0];
  List sprior   = ppprior[1];

  std::vector<std::string> types(nsub);
  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar), lp1(npar),
     lp2(npar), llower(npar), lupper(npar), sp1(npar), sp2(npar), slower(npar),
     supper(npar), bws(nsub),  pdists(npar), ldists(npar), sdists(npar);
  arma::uvec llog(npar), slog(npar), plog(npar),substore(nsub), npdas(nsub),
             gpuids(nsub);
  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);
  GetPrior(lprior, ldists, lp1, lp2, llower, lupper, llog);
  GetPrior(sprior, sdists, sp1, sp2, slower, supper, slog);

  // Extract subject level data before entering the loop.
  arma::field<arma::cube> subtheta(nsub); // nchains x npar x nmc; nmc x nchain
  arma::field<arma::mat> usetheta(nsub), lp(nsub), ll(nsub);
  arma::field<arma::vec> uselp(nsub), usell(nsub), allpars(nsub), RTs(nsub);
  arma::field<arma::umat> n1idx(nsub), cellidx(nsub);
  arma::field<arma::uvec> matchcells(nsub), emptycells(nsub), isr1(nsub);
  arma::uvec posdrift(nsub);
  arma::field<std::vector<std::string>> parnames(nsub), dim1s(nsub), dim2s(nsub),
               dim3s(nsub);
  arma::field<arma::ucube> models(nsub);

  TransformSubjects(samples_in, subtheta, usetheta, lp, uselp, ll, usell,
    substore, types, allpars, n1idx, matchcells, emptycells, cellidx, parnames,
    dim1s, dim2s, dim3s, isr1, posdrift, models, npdas, bws, gpuids, RTs);


  arma::cube theta0   = GetTheta0(samples_in); // nsub x npar x nchain
  // Hyper-level data: theta0 == ps (nchain x nsub x npar)
  // extract first nmc thetas from each participants

  arma::umat hyper_rejection_rate(nchain, nmc, arma::fill::zeros);
  arma::field<arma::umat> rejection_rate(nsub);
  InitializeSubjectRJ(samples_in, rejection_rate);
  arma::uvec hrj, rj; // nchain x nsamp,
  arma::mat uselp_tmp;

  for (size_t k = 0; k < nchain; k++) { // usetheta is nchain x npar
    usehll(k) = sumloghlike(theta0.slice(k), pdists,
      arma::trans(usephi(0).row(k)), arma::trans(usephi(1).row(k)),
      plower, pupper, plog);
  }


   for (size_t i = 1; i < nsamp; i++) {
       hrj = arma::zeros<arma::uvec>(nchain);

       if (R::runif(0, 1) < hpm) {

           MigrateDMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
             plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
             sdists, sp1, sp2, slower, supper, slog, rp, hrj, debug);

       } else if (R::runif(0, 1) < hpm0) {

         MigrateDMCHyperchains_old(usephi, usehlp, usehll, theta0, pdists,
                               plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
                               sdists, sp1, sp2, slower, supper, slog, rp, hrj, debug);

       } else {
         for (size_t l = 0; l < npar; l++) {
           CrossoverDMCHyperchains_blocked(usephi, usehlp, usehll, theta0, pdists,
                                           plower, pupper, plog, ldists, lp1,
                                           lp2, llower, lupper, llog,
                                           sdists, sp1, sp2, slower, supper,
                                           slog, rp, gammamult, hrj, l);
         }

          // CrossoverDMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
          //   plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
          //   sdists, sp1, sp2, slower, supper, slog, rp, gammamult, hrj);
       }

       // arma::mat A = arma::shuffle(usephi(0), 1);  // nchain x npar
       // arma::mat B = arma::shuffle(usephi(1), 1);
       // usephi(0) = A;
       // usephi(1) = B;

       /* ----------------------------------------------------------
        *  Update data level and hyper data
        * ----------------------------------------------------------*/
      for (size_t j = 0; j < nsub; j++) {
          // usephi: nchain x npar; usethetas(j): nchain x npar
          uselp(j) = UpdatePriors(usetheta(j), pdists, usephi(0), usephi(1),
             plower, pupper, plog);
          rj = arma::zeros<arma::uvec>(nchain);

          if (R::runif(0, 1) < pm) {

              MigrateDMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
                pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
                parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
                n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j),
                isr1(j), posdrift(j), npdas(j), bws(j), ncore, 0, rp, gammamult, rj, debug);
          } else if (R::runif(0, 1) < pm) {

            MigrateDMCDatachains_old(usetheta(j), uselp(j), usell(j), pnames,
                                 pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
                                 parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
                                 n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j),
                                 isr1(j), posdrift(j), npdas(j), bws(j), ncore, 0, rp, gammamult,
                                 rj, debug);


          } else {

            for (size_t l = 0; l < npar; l++) {
              CrossoverDMCDatachains_blocked(usetheta(j), uselp(j), usell(j), pnames,
                         pdists, usephi(0), usephi(1), plower, pupper,
                         plog, allpars(j),
                         parnames(j), models(j), types[j], dim1s(j),
                         dim2s(j), dim3s(j),
                         n1idx(j), emptycells(j), cellidx(j), RTs(j),
                         matchcells(j),
                         isr1(j), posdrift(j), npdas(j), bws(j), ncore, 0, rp,
                         gammamult, rj, l);
            }
//
//             CrossoverDMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
//               pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
//               parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
//               n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j),
//               isr1(j), posdrift(j), npdas(j), bws(j), ncore, 0, rp, gammamult, rj);
          }

          // theta0s: nsub x npar x nchain == ps: nchain x nsub x npar
          for (size_t k = 0; k < nchain; k++) { // usetheta is nchain x npar
            theta0.slice(k).row(j) = usetheta(j).row(k);
          }


          if ( i % thin == 0 ) {
            substore(j)++;
            // nmc x nchain; nchain x 1
            lp(j).row(substore(j)) = uselp(j).t();
            ll(j).row(substore(j)) = usell(j).t();  // usetheta: nchain x npar
            subtheta(j).slice(substore(j)) = usetheta(j);
            rejection_rate(j).col(substore(j)) = rj;
          }
      } // end of subject loop

      // for (size_t k = 0; k < nchain; k++) { // usetheta is nchain x npar
      //   usehll(k) = sumloghlike(theta0.slice(k), pdists,
      //     arma::trans(usephi(0).row(k)), arma::trans(usephi(1).row(k)),
      //     plower, pupper, plog);
      // }

      if (i % thin == 0) {
        store_i++;
        if ((store_i+1) % report == 0) Rcout << store_i + 1 << " ";
        hlp.row(store_i)  = usehlp.t(); // nmc x nchain = nchain x 1
        hll.row(store_i)  = usehll.t();
        location.slice(store_i) = usephi(0); // nchain x npar x nmc = nchain x npar
        scale.slice(store_i)    = usephi(1);
        hyper_rejection_rate.col(store_i) = hrj;
      }
   }   // end nsamp

   /* ----------------------------------------------------------
    *  Reconstruct data-level and then hyper-level
    * ----------------------------------------------------------*/
   for (size_t ii = 0; ii < nsub; ii++) {
     List subject = samples_in[ii];
     subject["summed_log_prior"] = lp(ii);
     subject["log_likelihoods"]  = ll(ii);
     subject["theta"]  = subtheta(ii);
     subject["rejection_rate"]   = rejection_rate(ii);
     samples_in[ii] = subject;
   }

   Rcpp::List newphi      = Rcpp::List::create(
     Rcpp::Named("location")   = location,
     Rcpp::Named("scale")      = scale);

   hyper["h_log_likelihoods"]  = hll;
   hyper["h_summed_log_prior"] = hlp;
   hyper["phi"]                = newphi;
   hyper["rejection_rate"]     = hyper_rejection_rate;

   samples_in.attr("hyper") = hyper;
   return samples_in;
}

// [[Rcpp::export]]
List run_hyper_dgmc(List samples, unsigned int report, double pm,
  double hpm, double qm, double hqm, double gammamult, unsigned int ngroup,
  unsigned int ncore) {

  List samples_in(clone(samples));
  CheckHyperPnames(samples_in);

  List hyper = samples_in.attr("hyper");
  unsigned int npar  = hyper["n.pars"];
  unsigned int nmc   = hyper["nmc"];
  unsigned int thin  = hyper["thin"];
  unsigned int startR= hyper["start"]; // start_R == 1;
  unsigned int nchain= hyper["n.chains"];
  double rp = hyper["rp"];            // rp is defined in initialise
  unsigned int nsub = samples.size();
  unsigned int startC = startR - 1;   // start_C == 0;
  unsigned int store_i = startC;      // store_i == 0;
  unsigned int nsamp = 1 + (nmc - startR) * thin;

  /* Hyperparameters, hyper logprior and hyper loglikelihood */
  std::vector<std::string> pnames = hyper("p.names");
  List phi      = hyper["phi"];
  arma::mat hlp = hyper["h_summed_log_prior"]; // nmc x nchain
  arma::mat hll = hyper["h_log_likelihoods"];
  arma::cube location = phi[0]; // nchain x npar x nmc
  arma::cube scale    = phi[1];
  // extract first nmc thetas from each participants; nsub x npar x nchain
  arma::cube theta0  = GetTheta0(samples_in);

  arma::field<arma::vec> blocks(npar); // reserved blocks variable
  arma::field<arma::mat> usephi(2);
  usephi(0) = location.slice(startC); // nchain x npar
  usephi(1) = scale.slice(startC);
  arma::vec usehlp = arma::trans(hlp.row(startC)); // nchain
  arma::vec usehll = arma::trans(hll.row(startC));

  List subject0 = samples_in[0];
  List pprior   = subject0["p.prior"];
  List ppprior  = hyper["pp.prior"];
  List lprior   = ppprior[0];
  List sprior   = ppprior[1];
  std::vector<std::string> types(nsub);
  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar), lp1(npar),
     lp2(npar), llower(npar), lupper(npar), sp1(npar), sp2(npar), slower(npar),
     supper(npar), bws(nsub), pdists(npar), ldists(npar), sdists(npar);
  arma::uvec llog(npar), slog(npar), plog(npar), substore(nsub), npdas(nsub),
              gpuids(nsub);
  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);
  GetPrior(lprior, ldists, lp1, lp2, llower, lupper, llog);
  GetPrior(sprior, sdists, sp1, sp2, slower, supper, slog);

  // Extract subject level data before entering for loop.
  arma::field<arma::cube> subtheta(nsub); // nchains x npar x nmc
  arma::field<arma::mat> usetheta(nsub), lp(nsub), ll(nsub); // nmc x nchain
  arma::field<arma::vec> uselp(nsub), usell(nsub), allpars(nsub), RTs(nsub);
  arma::field<arma::umat> n1idx(nsub), cellidx(nsub);
  arma::field<arma::uvec> matchcells(nsub), emptycells(nsub), isr1(nsub);
  arma::field<std::vector<std::string>> parnames(nsub), dim1s(nsub),
     dim2s(nsub), dim3s(nsub);
  arma::field<arma::ucube> models(nsub);
  arma::uvec posdrift(nsub);


  TransformSubjects(samples_in, subtheta, usetheta, lp, uselp, ll, usell,
    substore, types, allpars, n1idx, matchcells, emptycells, cellidx, parnames,
    dim1s, dim2s, dim3s, isr1, posdrift, models, npdas, bws, gpuids, RTs);


  arma::umat hyper_rejection_rate(nchain, nmc, arma::fill::zeros);
  arma::field<arma::umat> rejection_rate(nsub);
  InitializeSubjectRJ(samples_in, rejection_rate);
  arma::uvec hrj, rj; // nchain x nsamp,

  for (size_t i = 1; i < nsamp; i++) {

      hrj = arma::zeros<arma::uvec>(nchain);

      if (R::runif(0, 1) < hpm) {
           MigrateDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
              plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
              sdists, sp1, sp2, slower, supper, slog, ngroup, rp, hrj);
      } else if (R::runif(0, 1) <= hqm) {
           MutateDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists, plower,
              pupper, plog, ldists, lp1, lp2, llower, lupper, llog, sdists, sp1,
              sp2, slower, supper, slog, ngroup, rp, hrj);
      } else {
           CrossoverDGMCHyperchains(usephi, usehlp, usehll, theta0, pdists,
              plower, pupper, plog, ldists, lp1, lp2, llower, lupper, llog,
              sdists, sp1, sp2, slower, supper, slog, ngroup, rp, gammamult,
              hrj);
      }

      /* ----------------------------------------------------------
       *  Update data level and hyper data
       * ----------------------------------------------------------*/
      for (size_t j = 0; j < nsub; j++) {
        uselp(j) = UpdatePriors(usetheta(j), pdists, usephi(0), usephi(1),
          plower, pupper, plog); // usephi: nchain x npar
        rj = arma::zeros<arma::uvec>(nchain);


        if (R::runif(0, 1) < pm) {
            MigrateDGMCDatachains(usetheta(j), uselp(j), usell(j), pnames,
              pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
              parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
              n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j),
              isr1(j), posdrift(j), ngroup, rp, npdas(j), bws(j), ncore, 0, rj);
        } else if (R::runif(0, 1) <= qm) {
            MutateDGMCDataChains(usetheta(j), uselp(j), usell(j), pnames,
              pdists, usephi(0), usephi(1), plower, pupper, plog, allpars(j),
              parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
              n1idx(j), emptycells(j), cellidx(j), RTs(j), matchcells(j),
              isr1(j), posdrift(j),
              ngroup, npdas(j), bws(j), ncore, 0, rp, rj);
        } else {
            CrossoverDGMCDatachains(usetheta(j), uselp(j), usell(j),
              pnames, pdists, usephi(0), usephi(1), plower, pupper, plog,
              allpars(j), parnames(j), models(j), types[j], dim1s(j), dim2s(j),
              dim3s(j), n1idx(j), emptycells(j), cellidx(j), RTs(j),
              matchcells(j), isr1(j), posdrift(j), ngroup, npdas(j), bws(j), ncore, 0, rp,
              gammamult, rj);
        }


        if ( i % thin == 0 ) {
          substore(j)++;
          lp(j).row(substore(j)) = uselp(j).t();  // nmc x nchain; nchain x 1
          ll(j).row(substore(j)) = usell(j).t();  // usetheta: nchain x npar
          subtheta(j).slice(substore(j)) = usetheta(j);
          rejection_rate(j).col(substore(j)) = rj;
        }

        for (size_t k = 0; k < nchain; k++) { // usetheta is nchain x npar
          theta0.slice(k).row(j) = usetheta(j).row(k);
        }
      }  // end of subject loop

      if (i % thin == 0) {
         store_i++;
         if ((store_i+1) % report == 0) Rcout << store_i + 1 << " ";
         hlp.row(store_i)  = usehlp.t(); // nmc x nchain = nchain x 1
         hll.row(store_i)  = usehll.t();
         location.slice(store_i) = usephi(0);
         scale.slice(store_i)    = usephi(1);
         hyper_rejection_rate.col(store_i) = hrj;
      }
  } // end of MCMC iteration


  /* ----------------------------------------------------------
   *  Reconstruct data-level and then hyper-level
   * ----------------------------------------------------------*/
  for (size_t ii = 0; ii < nsub; ii++) {
    List subject = samples_in[ii];
    subject["summed_log_prior"] = lp(ii);
    subject["log_likelihoods"]  = ll(ii);
    subject["theta"]            = subtheta(ii);
    subject["rejection_rate"]   = rejection_rate(ii);
    samples_in[ii] = subject;
  }

  List newphi      = Rcpp::List::create(
    Named("location")   = location,
    Named("scale")      = scale);
  hyper["h_log_likelihoods"]  = hll;
  hyper["h_summed_log_prior"] = hlp;
  hyper["phi"]                = newphi;
  hyper["rejection_rate"]     = hyper_rejection_rate;

  samples_in.attr("hyper") = hyper;
  return samples_in;
}

