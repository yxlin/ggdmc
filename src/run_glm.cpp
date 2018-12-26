#include <ggdmc.hpp>
#include <chrono>
#include <random>
using namespace Rcpp;

void Crossover_glm(arma::mat& usetheta,    // nchain x npar
  arma::vec& uselp, arma::vec& usell,  // nchain x 1
  std::vector<std::string> pnames,
  arma::vec dists, arma::vec p1, arma::vec p2,
  arma::vec lower, arma::vec upper, arma::uvec lg,
  std::string type, arma::mat Y,
  arma::mat X, double rp, double gammamult, bool force,
  arma::uvec& rj)
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
  chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));
  arma::vec noise(npar);

  for (size_t i = 0; i < nchain; i++) {
    // Among the rest of the subpopulation xi \ {x^i_j}, sample two
    // chain/chromosomes, (xi_s and xi_t)
    subchains = PickChains(chains(i), 2, chains);
    k0 = chains(i);
    k1 = subchains(0);
    k2 = subchains(1);

    for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);

    tmp = theta.col(k0) + (gamma % (theta.col(k1) - theta.col(k2))) + noise;

    // PDA re-calculation
    if (force) {
      uselp(k0) = sumlogprior(theta.col(k0), dists, p1, p2, lower, upper, lg);
      usell(k0) = sumloglike_glm(theta.col(k0), type, X, Y);
    }
    cur_logpos = usell(k0) + uselp(k0);

    tmp_lp = sumlogprior(tmp, dists, p1, p2, lower, upper, lg);
    tmp_ll = sumloglike_glm(tmp, type, X, Y);
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



void Migrate_glm(arma::mat& usetheta, arma::vec& uselp, arma::vec& usell,
  std::vector<std::string> pnames, arma::vec dists, arma::vec p1,
  arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::mat Y,
  arma::uvec matchcell, arma::uvec isr1, arma::mat X, bool posdrift, double rp,
  double gammamult, bool force, unsigned int nsim, double bw, unsigned int ncore,
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
    tmp_ll = sumloglike_glm(theta_star, type, X, Y);
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


void Migrate_old_glm(arma::mat& usetheta, arma::vec& uselp,
                              arma::vec& usell,
  std::vector<std::string> pnames, arma::vec dists, arma::vec p1,
  arma::vec p2, arma::vec lower, arma::vec upper, arma::uvec islog,
  arma::vec allpar, std::vector<std::string> parnames, arma::ucube model,
  std::string type, std::vector<std::string> dim1,
  std::vector<std::string> dim2, std::vector<std::string> dim3,
  arma::umat n1idx, arma::uvec ise, arma::umat cellidx, arma::mat Y,
  arma::uvec matchcell, arma::uvec isr1, arma::mat X, bool posdrift,
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
    tmp_ll(i) = sumloglike_glm(tmp.col(i), type, X, Y);
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
List run_glm(List samples, arma::uvec force, unsigned int report, double pm,
             double pm0, double gammamult, unsigned int ncore) {

  List samples_in(clone(samples)); // so R' original samples stays
  CheckPnames(samples_in);

  List data          = samples_in["data"];
  double rp          = samples_in["rp"];
  arma::ucube model  = data.attr("model");
  unsigned int npda  = data.attr("n.pda");
  double bw          = data.attr("bw");
  bool debug         = data.attr("debug");
  unsigned int gpuid = data.attr("gpuid");

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

  arma::mat Y      = data.attr("Y");
  arma::mat X_     = data.attr("X"); // glm independent / predictor variable
  arma::mat X = arma::join_horiz(arma::ones(X_.n_rows), X_);

  // arma::vec Y_     = data["Y"];
  // arma::vec X_     = data["X"]; // glm independent / predictor variable
  // arma::mat Y = Y_;
  // arma::mat X = arma::join_horiz(arma::ones(X_.n_rows), X_);


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

      Migrate_glm(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
                       plower, pupper, plog, allpar, parnames, model, type,
                       dim1, dim2, dim3, n1idx, ise, cellidx, Y, mc, isr1,
                       X, posdrift, rp, gammamult, force(i), npda, bw, ncore,
                       gpuid, rj, debug);

    } else if (R::runif(0.0, 1.0) < pm0) {
      Migrate_old_glm(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
                           plower, pupper, plog, allpar, parnames, model, type,
                           dim1, dim2, dim3, n1idx, ise, cellidx, Y, mc, isr1,
                           X, posdrift, rp, gammamult, force(i), npda, bw,
                           ncore, gpuid, rj, debug);
    } else {
      Crossover_glm(usetheta, uselp, usell, pnames, pdists, pp1, pp2,
                    plower, pupper, plog, type, Y, X, rp, gammamult,
                    force(i), rj);
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

void MigrateHyper_glm(arma::field<arma::mat>& usephi,
                           arma::vec& usehlp,
                           arma::vec& usehll,
                           arma::cube theta,
                           arma::vec pdists,
                           arma::vec lower, arma::vec upper, arma::uvec lg,
                           arma::vec ldists, arma::vec lp1, arma::vec lp2,
                           arma::vec llower, arma::vec lupper, arma::uvec llg,
                           arma::vec sdists, arma::vec sp1, arma::vec sp2,
                           arma::vec slower, arma::vec supper, arma::uvec slg,
                           double rp, arma::uvec& rj) {

  arma::mat useloc = arma::trans(usephi(0));
  arma::mat usesca = arma::trans(usephi(1));
  arma::vec usehlp_ = usehlp;
  arma::vec usehll_ = usehll;

  unsigned int nchain = useloc.n_cols;
  unsigned int npar   = useloc.n_rows;
  arma::uvec subchains = GetSubchains(nchain, false);
  unsigned int nsubchain = subchains.n_elem;
  arma::vec cur_hlp(nsubchain), cur_hll(nsubchain), cur_loc, cur_sca,
            tmp_loc(npar), tmp_sca(npar);

  unsigned int next_chain, k;
  double tmp_hlp, tmp_hll, tmp_logpos, cur_logpos, mh;

  arma::uvec idx1 = arma::find_finite(lp1);

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

    tmp_hlp = sumloghprior(tmp_loc(idx1), tmp_sca(idx1), ldists(idx1),
                           sdists(idx1), lp1(idx1), sp1(idx1), lp2(idx1),
                           sp2(idx1), llower(idx1), slower(idx1), lupper(idx1),
                           supper(idx1), llg(idx1), slg(idx1));
    tmp_hll = sumloghlike(theta.slice(k), pdists, tmp_loc, tmp_sca, lower,
                          upper, lg);
    tmp_logpos = tmp_hlp + tmp_hll;
    if (std::isnan(tmp_logpos)) tmp_logpos = R_NegInf;

    /* CRITICAL. Update usehll for new theta; nsub x npar x nchain */
    usehll_(next_chain) = sumloghlike(theta.slice(next_chain), pdists,
            useloc.col(next_chain), usesca.col(next_chain), lower, upper, lg);
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

void MigrateHyper_old_glm(arma::field<arma::mat>& usephi,
                          arma::vec& usehlp,
                          arma::vec& usehll,
                          arma::cube theta,
                          arma::vec pdists,
                          arma::vec lower, arma::vec upper, arma::uvec lg,
                          arma::vec ldists, arma::vec lp1, arma::vec lp2,
                          arma::vec llower, arma::vec lupper, arma::uvec llg,
                          arma::vec sdists, arma::vec sp1, arma::vec sp2,
                          arma::vec slower, arma::vec supper, arma::uvec slg,
                          double rp, arma::uvec& rj) {

  arma::mat useloc = arma::trans(usephi(0)); // useloc: npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  arma::vec usehlp_ = usehlp;
  arma::vec usehll_ = usehll;

  unsigned int npar    = useloc.n_rows;
  unsigned int nchain  = useloc.n_cols;
  arma::uvec subchains = GetSubchains(nchain, false);
  unsigned int nsubchain = subchains.n_elem;
  unsigned int k;

  double tmp_logpos, cur_logpos, mh;
  arma::mat tmp_loc(npar, nsubchain), tmp_sca(npar, nsubchain);
  arma::vec cur_hlp(nsubchain), cur_hll(nsubchain);
  arma::vec tmp_hlp(nsubchain), tmp_hll(nsubchain), noise(npar);

  arma::uvec idx1 = arma::find_finite(lp1);

  arma::vec tmp_loc_, tmp_sca_;

  for(size_t i = 0; i < nsubchain; i++) { // 0, 1, 2, 5, ...
    k = subchains(i);
    // Uniform(-rp, rp); (b-a) * R::runif(1) + a;
    for(size_t j = 0; j < npar; j++) noise(j) = R::runif(-rp, rp);
    tmp_loc.col(i) = useloc.col(k) + noise; // proposal
    tmp_sca.col(i) = usesca.col(k) + noise; // proposal

    /* CRITICAL. Update usehll for new theta; nsub x npar x nchain */
    usehll_(k) = sumloghlike(theta.slice(k), pdists, useloc.col(k),
            usesca.col(k), lower, upper, lg); // this changes too
    cur_hlp(i) = usehlp_(k);
    cur_hll(i) = usehll_(k);

    arma::vec tmp_loc_ = tmp_loc.col(i);
    arma::vec tmp_sca_ = tmp_sca.col(i);
    // theta: nsub x npar x nchain
    tmp_hlp(i) = sumloghprior(tmp_loc_(idx1), tmp_sca_(idx1), ldists(idx1),
            sdists(idx1), lp1(idx1), sp1(idx1), lp2(idx1), sp2(idx1),
            llower(idx1), slower(idx1), lupper(idx1), supper(idx1), llg(idx1),
            slg(idx1));

    tmp_hll(i) = sumloghlike(theta.slice(k), pdists, tmp_loc.col(i),
            tmp_sca.col(i), lower, upper, lg);
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

void CrossoverHyper_glm(arma::field<arma::mat>& usephi, arma::vec& usehlp,
               arma::vec& usehll, arma::cube theta,
               arma::vec pdist, arma::vec p1, arma::vec p2,
               arma::vec lower, arma::vec upper, arma::uvec lg,
               arma::vec ldist, arma::vec lp1, arma::vec lp2,
               arma::vec llower, arma::vec lupper, arma::uvec llg,
               arma::vec sdist, arma::vec sp1, arma::vec sp2,
               arma::vec slower, arma::vec supper, arma::uvec slg,
               double rp, double gammaMult, arma::uvec& rj, unsigned int j)
{
  double cur_logpos, tmp_logpos, tmp_hlp, tmp_llik, tmp_slik, tmp_hll;
  unsigned int nchain = usephi(0).n_rows;
  double hgamma = 1.19;
  arma::uvec subchains;
  unsigned int k0;
  arma::uvec chains_ = arma::linspace<arma::uvec>(0, nchain - 1, nchain);
  arma::uvec chains = arma::shuffle(chains_);
  arma::mat useloc = arma::trans(usephi(0)); // npar x nchain
  arma::mat usesca = arma::trans(usephi(1));

  arma::uvec pidx_fin = arma::find_finite(pdist);
  arma::uvec lidx_fin = arma::find_finite(ldist);
  arma::uvec sidx_fin = arma::find_finite(sdist);
  arma::uvec pidx_non = arma::find_nonfinite(pdist);
  arma::uvec lidx_non = arma::find_nonfinite(ldist);
  arma::uvec sidx_non = arma::find_nonfinite(sdist);
  unsigned int nppar = pidx_fin.n_elem;
  unsigned int nlpar = lidx_fin.n_elem;
  unsigned int nspar = sidx_fin.n_elem;

  arma::mat thetak0;

  arma::vec useloc_, usesca_, diff1, diff2, tmp_loc, tmp_sca;

  for (size_t i = 0; i < nchain; i++) {
    k0 = chains(i);
    useloc_ = useloc.col(k0); // three elements;
    usesca_ = usesca.col(k0); // two elements (+ NA)
    tmp_loc = useloc_;  // make a copy
    tmp_sca = usesca_;
    thetak0 = theta.slice(k0);
    // sum log h like; theta = nsub x npar x nchain
    usehll(k0) = sumloghlike(thetak0.cols(pidx_fin), pdist(pidx_fin), useloc_(pidx_fin),
           usesca_(pidx_fin), lower(pidx_fin), upper(pidx_fin), lg(pidx_fin));
    cur_logpos = usehlp(k0) + usehll(k0);

    subchains = PickChains(k0, 2, chains); // (b-a) * R::runif(1) + a;
    diff1 = useloc.col(subchains(0)) - useloc.col(subchains(1)); // 3
    diff2 = usesca.col(subchains(0)) - usesca.col(subchains(1)); // 2

    tmp_loc(j) = useloc_(j) + R::runif(-rp, rp) + hgamma * diff1(j);
    tmp_sca(j) = usesca_(j) + R::runif(-rp, rp) + hgamma * diff2(j);

    // sum log h prior
    tmp_llik = sumlogprior(tmp_loc(lidx_fin), ldist(lidx_fin), lp1(lidx_fin),
                           lp2(lidx_fin), llower(lidx_fin), lupper(lidx_fin),
                           llg(lidx_fin));
    tmp_slik = sumlogprior(tmp_sca(sidx_fin), sdist(sidx_fin), sp1(sidx_fin),
                           sp2(sidx_fin), slower(sidx_fin), supper(sidx_fin),
                           slg(sidx_fin));
    tmp_hlp = tmp_llik + tmp_slik;

    // tmp_sca(idx1).print("tmp_sca idx1");
    // sp1(idx1).print("sp1 idx1");
    // sp2(idx1).print("sp2 idx1");

    // tmp_loc(lidx_non) = p1(pidx_fin);
    // tmp_sca(sidx_non) = p2(pidx_fin);
    tmp_hll = sumloghlike(thetak0.cols(pidx_fin), pdist(pidx_fin),
                          tmp_loc(pidx_fin), tmp_sca(pidx_fin),
                          lower(pidx_fin), upper(pidx_fin), lg(pidx_fin));
    tmp_logpos = tmp_hlp + tmp_hll;
    double mh = std::exp(tmp_logpos - cur_logpos);

    if (!std::isnan(mh) && (R::runif(0, 1) < mh) )  {
         // only has two columns
          usephi(0).row(k0).col(j) = tmp_loc(j);
          usephi(1).row(k0).col(j) = tmp_sca(j);
          usehlp(k0) = tmp_hlp;
          usehll(k0) = tmp_hll;
       rj(k0) = 3;
    } else {
      rj(k0) = 4;
    }
  }    // chain loop
}      // function closing

void CrossoverData_glm(arma::mat& usetheta,    // nchain x npar
                       arma::vec& uselp, arma::vec& usell,
                       std::vector<std::string> pnames,
                       arma::vec dist, arma::vec p1, arma::vec p2,
                       arma::vec ldist, arma::vec sdist,
                       arma::mat usephi0, arma::mat usephi1,
                       arma::vec lower, arma::vec upper, arma::uvec lg,
                       arma::vec allpar,
                       std::vector<std::string> parnames, arma::ucube model,
                       std::string type,
                       std::vector<std::string> dim1,
                       std::vector<std::string> dim2,
                       std::vector<std::string> dim3, arma::umat n1idx,
                       arma::uvec ise,
                       arma::umat cellidx, arma::mat Y, arma::uvec matchcell,
                       arma::uvec isr1, arma::mat X, double rp,
                       double gammamult, arma::uvec& rj) {

  unsigned int k0, k1, k2;
  double tmp_lp, tmp_ll, tmp_logpos, cur_logpos, mh;
  arma::uvec chains, subchains;
  arma::vec tmp;
  unsigned int npar = usetheta.n_cols;

  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain
  arma::mat useloc = arma::trans(usephi0);  // npar x nchain
  arma::mat usesca = arma::trans(usephi1);  // npar x nchain

  unsigned int nchain = theta.n_cols;
  // double gamma = 2.38 / std::sqrt(2.0);
  arma::vec gamma  = GetGamma(npar, gammamult);
  chains = arma::shuffle(arma::linspace<arma::uvec>(0, nchain - 1, nchain));
  arma::vec noise(npar);

  arma::vec phi0k0, phi1k0;

  arma::uvec pidx_fin = arma::find_finite(dist);
  arma::uvec lidx_fin = arma::find_finite(ldist);
  arma::uvec sidx_fin = arma::find_finite(sdist);
  arma::uvec pidx_non = arma::find_nonfinite(dist);
  arma::uvec lidx_non = arma::find_nonfinite(ldist);
  arma::uvec sidx_non = arma::find_nonfinite(sdist);

  arma::vec thetak0, theta0, theta1, tmp_theta;
  for (size_t i = 0; i < nchain; i++) {
    subchains = PickChains(chains(i), 2, chains);
    k0 = chains(i);
    k1 = subchains(0);
    k2 = subchains(1);

    thetak0 = theta.col(k0); // npar x nchain
    tmp_theta = thetak0;
    cur_logpos = uselp(k0) + usell(k0);

    theta0 = theta.col(k1);
    theta1 = theta.col(k2);
    tmp_theta(pidx_fin) = thetak0(pidx_fin) + gamma(pidx_fin) % (theta0(pidx_fin) - theta1(pidx_fin)) +
      R::runif(-rp, rp) ;
    // tmp_theta = theta_ + gamma % (theta0 - theta1) +
    //   R::runif(-rp, rp) ;
    // phi0.col(i).print("phi0.col ");
    phi0k0 = useloc.col(i);
    phi1k0 = usesca.col(i);

    tmp_lp = sumlogprior(tmp_theta(pidx_fin), dist(pidx_fin),
                         phi0k0(pidx_fin),
                         phi1k0(pidx_fin), lower(pidx_fin), upper(pidx_fin),
                         lg(pidx_fin));
    // tmp_theta(pidx_non) = phi0k0(pidx_non);
    tmp_ll = sumloglike_glm(tmp_theta, type, X, Y);
    tmp_logpos = tmp_ll + tmp_lp;
    if (std::isnan(tmp_logpos)) tmp_logpos = R_NegInf;

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

void MigrateData_glm(arma::mat& usetheta,    // nchain x npar
                     arma::vec& uselp, arma::vec& usell,
                     std::vector<std::string> pnames,
                     arma::vec dists,
                     arma::mat usephi0, arma::mat usephi1, // nchain x npar
                     arma::vec lower, arma::vec upper, arma::uvec lg,
                     arma::vec allpar, std::vector<std::string> parnames,
                     arma::ucube model, std::string type,
                     std::vector<std::string> dim1,
                     std::vector<std::string> dim2,
                     std::vector<std::string> dim3,
                     arma::umat n1idx, arma::uvec ise, arma::umat cellidx,
                     arma::mat Y, arma::uvec matchcell, arma::uvec isr1,
                     arma::mat X, double rp, double gammamult, arma::uvec& rj)
{
  arma::mat phi0  = arma::trans(usephi0);
  arma::mat phi1  = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  unsigned int npar      = theta.n_rows;
  unsigned int nchain    = theta.n_cols;
  arma::uvec subchains   = GetSubchains(nchain, false); // eg. 0, 1, 3, 4, 8
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
                         upper, lg);
    tmp_ll = sumloglike_glm(theta_star, type, X, Y);
    // if (std::isinf(tmp_lp) && tmp_lp > 0.0) tmp_lp = 1e-10;

    tmp_logpos = tmp_lp + tmp_ll;
    cur_logpos = uselp(next_chain) + usell(next_chain);
    if (std::isnan(tmp_logpos)) tmp_logpos = R_NegInf;
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

void MigrateData_old_glm(arma::mat& usetheta,    // nchain x npar
                      arma::vec& uselp, arma::vec& usell,
                      std::vector<std::string> pnames,
                      arma::vec dists,
                      arma::mat usephi0, arma::mat usephi1, // nchain x npar
                      arma::vec lower, arma::vec upper, arma::uvec lg,
                      arma::vec allpar, std::vector<std::string> parnames,
                      arma::ucube model, std::string type,
                      std::vector<std::string> dim1,
                      std::vector<std::string> dim2,
                      std::vector<std::string> dim3,
                      arma::umat n1idx, arma::uvec ise, arma::umat cellidx,
                      arma::mat Y, arma::uvec matchcell, arma::uvec isr1,
                      arma::mat X, double rp, double gammamult, arma::uvec& rj)
{

  arma::mat phi0  = arma::trans(usephi0);
  arma::mat phi1  = arma::trans(usephi1);
  arma::mat theta = arma::trans(usetheta); // theta: npar x nchain;
  arma::vec uselp_ = uselp;
  arma::vec usell_ = usell;
  unsigned int npar    = theta.n_rows;
  unsigned int nchain  = theta.n_cols;
  arma::uvec subchains = GetSubchains(nchain, false); // eg. 0, 1, 3, 4, 8
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
           upper, lg);
    tmp_ll(i) = sumloglike_glm(tmp.col(i), type, X, Y);

  }

  tmp_logpos = tmp_ll(nsubchain - 1) + tmp_lp(nsubchain - 1);
  if (std::isnan(tmp_logpos)) tmp_logpos = R_NegInf;

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

arma::vec UpdatePriors_glm(arma::field<arma::mat>& usephi, arma::mat& theta,
                            arma::vec pdist, arma::vec p1, arma::vec p2,
                            arma::vec ldist, arma::vec sdist,
                             arma::vec lower,
                            arma::vec upper, arma::uvec lg) {

  arma::mat useloc = arma::trans(usephi(0)); // npar x nchain
  arma::mat usesca = arma::trans(usephi(1));
  arma::mat theta_ = arma::trans(theta); // npar x nchain

  unsigned int nchain = theta.n_rows;
  unsigned int npar   = theta.n_cols;
  arma::vec out(nchain);

  arma::uvec pidx_fin, lidx_fin, sidx_fin, pidx_non, lidx_non, sidx_non;
  pidx_fin = arma::find_finite(pdist);
  lidx_fin = arma::find_finite(ldist);
  sidx_fin = arma::find_finite(sdist);
  pidx_non = arma::find_nonfinite(p1);
  lidx_non = arma::find_nonfinite(ldist);
  sidx_non = arma::find_nonfinite(sdist);

  for (size_t i = 0; i < nchain; i++) {
      // theta_.col(i).row(2) = useloc.col(i).row(2);
      arma::vec tmp1 = theta_.col(i);
      arma::vec tmp2 = useloc.col(i);
      arma::vec tmp3 = usesca.col(i);

    out(i) = sumlogprior(tmp1(pidx_fin), pdist(pidx_fin), tmp2(pidx_fin),
        tmp3(pidx_fin), lower(pidx_fin), upper(pidx_fin), lg(pidx_fin));
  }
//  theta = arma::trans(theta_);

  return out ;
}


// [[Rcpp::export]]
List run_hglm(List samples, unsigned int report, double pm, double pm0,
                   double hpm, double hpm0, double gammamult) {

  List samples_in(clone(samples));
  CheckHyperPnames(samples_in);

  List hyper = samples_in.attr("hyper");
  unsigned int nchain = hyper["n.chains"];
  unsigned int nmc    = hyper["nmc"];
  unsigned int thin   = hyper["thin"];
  unsigned int start  = hyper["start"]; // start_R == 1;
  unsigned int nsub = samples.size();
  unsigned int start_C = start - 1;   // start_C == 0;
  unsigned int store_i = start_C;    // store_i == 0;
  unsigned int nsamp = 1 + (nmc - start) * thin;
  double rp = hyper["rp"];            // rp is defined in initialise

  /* data_hyper/thetas (nsub x npar x nchain) == cps (nchain x nsub x npar) */
  std::vector<std::string> pnames = hyper("p.names");
  List phi      = hyper["phi"];
  arma::mat hlp = hyper["h_summed_log_prior"]; // nmc x nchain
  arma::mat hll = hyper["h_log_likelihoods"];
  arma::cube location = phi[0]; // nchain x npar x nmc
  arma::cube scale    = phi[1];


  List subject0 = samples_in[0];
  unsigned int nhpar  = hyper["n.pars"];
  unsigned int npar = subject0["n.pars"];
  List pprior  = subject0["p.prior"];
  List ppprior = hyper["pp.prior"];     /* Extract pprior & ppprior */
  List lprior  = ppprior[0];
  List sprior  = ppprior[1];
  std::vector<std::string> types(nsub);

  arma::vec p1, p2, lower, upper,
            loc_p1, loc_p2, loc_l, loc_u,
            sca_p1, sca_p2, sca_l, sca_u,
            dist_pp, dist_lp, dist_sp;
  arma::uvec lg, loc_lg, sca_lg;
  GetPrior(pprior, dist_pp, p1, p2, lower, upper, lg);
  GetPrior(lprior, dist_lp, loc_p1, loc_p2, loc_l, loc_u, loc_lg);
  GetPrior(sprior, dist_sp, sca_p1, sca_p2, sca_l, sca_u, sca_lg);

  // Extract subject level data before entering the loop.
  arma::field<arma::cube> subtheta(nsub); // nchains x npar x nmc; nmc x nchain
  arma::field<arma::mat> usetheta(nsub), lp(nsub), ll(nsub), Xs(nsub), Ys(nsub);
  arma::field<arma::vec> uselp(nsub), usell(nsub), allpars(nsub);
  arma::field<arma::umat> n1idx(nsub), cellidx(nsub);
  arma::field<arma::uvec> matchcells(nsub), emptycells(nsub), isr1(nsub);

  arma::field<std::vector<std::string>> parnames(nsub), dim1s(nsub),
                                        dim2s(nsub), dim3s(nsub);
  arma::field<arma::ucube> models(nsub);

  arma::uvec posdrift(nsub), substore(nsub), npdas(nsub), gpuids(nsub);
  arma::vec bws(nsub);

  TransformSubjects_glm(samples_in, subtheta, usetheta, lp, uselp, ll, usell,
    substore, types, allpars, n1idx, matchcells, emptycells, cellidx, parnames,
    dim1s, dim2s, dim3s, isr1, posdrift, models, npdas, bws, gpuids, Ys, Xs);

  // nsub x npar x nchain; extract first nmc thetas from each participants
  arma::cube theta0 = GetTheta0(samples_in);

  arma::umat hyper_rejection_rate(nchain, nmc, arma::fill::zeros);
  arma::field<arma::umat> rejection_rate(nsub);
  InitializeSubjectRJ(samples_in, rejection_rate);
  arma::uvec hrj, rj; // nchain x nsamp,

  arma::uvec pidx_fin = arma::find_finite(dist_pp);
  arma::uvec lidx_fin = arma::find_finite(dist_lp);
  arma::uvec sidx_fin = arma::find_finite(dist_sp);
  arma::uvec pidx_non = arma::find_nonfinite(dist_pp);
  arma::uvec lidx_non = arma::find_nonfinite(dist_lp);
  arma::uvec sidx_non = arma::find_nonfinite(dist_sp);

  arma::field<arma::mat> usephi(2);
  arma::mat usephi0_ = location.slice(start_C); // nchain x nhpar
  arma::mat usephi1_ = scale.slice(start_C);
  usephi(0) = usephi0_; // nchain x nhpar
  usephi(1) = usephi1_;
  arma::vec usehlp = arma::trans(hlp.row(start_C)); // nchain
  arma::vec usehll = arma::trans(hll.row(start_C));


  for (size_t i = 1; i < nsamp; i++) {
          hrj = arma::zeros<arma::uvec>(nchain);

       if (R::runif(0, 1) < hpm) {

           MigrateHyper_glm(usephi, usehlp, usehll, theta0,
             dist_pp, lower, upper, lg,
             dist_lp, loc_p1, loc_p2, loc_l, loc_u, loc_lg,
             dist_sp, sca_p1, sca_p2, sca_l, sca_u, sca_lg,
             rp, hrj);

       } else if (R::runif(0, 1) < hpm0) {

         MigrateHyper_old_glm(usephi, usehlp, usehll, theta0,
             dist_pp, lower, upper, lg,
             dist_lp, loc_p1, loc_p2, loc_l, loc_u, loc_lg,
             dist_sp, sca_p1, sca_p2, sca_l, sca_u, sca_lg,
             rp, hrj);

       } else {
                 for (size_t l = 0; l < nhpar; l++) {

                    CrossoverHyper_glm(usephi, usehlp, usehll, theta0,
                        dist_pp, p1, p2, lower, upper, lg,
                        dist_lp, loc_p1, loc_p2, loc_l, loc_u, loc_lg,
                        dist_sp, sca_p1, sca_p2, sca_l, sca_u, sca_lg,
                        rp, gammamult, hrj, l);
                 }
        }
  //      /* ----------------------------------------------------------
  //       *  Update data level and hyper data
  //       * ----------------------------------------------------------*/
           for (size_t j = 0; j < nsub; j++) {
  //            // usephi: nchain x npar; usethetas(j): nchain x npar
                uselp(j) = UpdatePriors_glm(usephi, usetheta(j), dist_pp, p1, p2,
                      dist_lp, dist_sp, lower, upper, lg);
  //
                rj = arma::zeros<arma::uvec>(nchain);
  //
          if (R::runif(0, 1) < pm) {
              MigrateData_glm(usetheta(j), uselp(j), usell(j), pnames,
                dist_pp, usephi(0), usephi(1), lower, upper, lg, allpars(j),
                parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
                n1idx(j), emptycells(j), cellidx(j), Ys(j), matchcells(j),
                isr1(j), Xs(j), rp, gammamult, rj);

          } else if (R::runif(0, 1) < pm) {
            MigrateData_old_glm(usetheta(j), uselp(j), usell(j), pnames,
                dist_pp, usephi(0), usephi(1), lower, upper, lg, allpars(j),
                parnames(j), models(j), types[j], dim1s(j), dim2s(j), dim3s(j),
                n1idx(j), emptycells(j), cellidx(j), Ys(j), matchcells(j),
                isr1(j), Xs(j), rp, gammamult, rj);

          } else {
            CrossoverData_glm(usetheta(j), uselp(j), usell(j), pnames,
                              dist_pp, p1, p2, dist_lp, dist_sp,
                 usephi(0), usephi(1), lower, upper, lg, allpars(j),
                 parnames(j), models(j), types[j], dim1s(j), dim2s(j),
                 dim3s(j), n1idx(j), emptycells(j), cellidx(j), Ys(j),
                 matchcells(j), isr1(j), Xs(j), rp, gammamult, rj);

         }
  //
  //          // theta0s: nsub x npar x nchain == ps: nchain x nsub x npar
           for (size_t k = 0; k < nchain; k++) { // usetheta is nchain x npar
             theta0.slice(k).row(j) = usetheta(j).row(k);
           }
  //
  //
          if ( i % thin == 0 ) {
            substore(j)++;
            // nmc x nchain; nchain x 1
            lp(j).row(substore(j)) = uselp(j).t();
            ll(j).row(substore(j)) = usell(j).t();  // usetheta: nchain x npar
            subtheta(j).slice(substore(j)) = usetheta(j);
            rejection_rate(j).col(substore(j)) = rj;
          }
        } // end of subject loop
  //
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
  //
  //  /* ----------------------------------------------------------
  //   *  Reconstruct data-level and then hyper-level
  //   * ----------------------------------------------------------*/
   for (size_t ii = 0; ii < nsub; ii++) {
     List subject = samples_in[ii];
     subject["summed_log_prior"] = lp(ii);
     subject["log_likelihoods"]  = ll(ii);
     subject["theta"]  = subtheta(ii);
     subject["rejection_rate"]   = rejection_rate(ii);
     samples_in[ii] = subject;
   }
  //
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

