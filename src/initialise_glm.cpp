#include <ggdmc.hpp>
using namespace Rcpp;


// [[Rcpp::export]]
List init_new_glm(unsigned int nmc, List data, List start, List prior,
                  double rp, unsigned int thin, unsigned int nchain) {

  unsigned int npar = prior.size();
  std::vector<std::string> pnames = prior.attr("names");

  // extract data-model options
  NumericVector modelAttr = data.attr("model"); // To extract attributes
  std::string type = modelAttr.attr("type");    // model type: DDM or LBA
  arma::mat Y      = data.attr("Y");
  arma::mat X_     = data.attr("X"); // glm independent / predictor variable
  arma::mat X = arma::join_horiz(arma::ones(X_.n_rows), X_);

  /*---------------------------------*/
  // should past p1 and p2 value from hierarchical level, if called by
  // init_newhier
  arma::vec pp1, pp2, plower, pupper, dist_pstart, p1, p2, lower, upper, dist_pp;
  arma::uvec plg, lg;
  GetPrior(start, dist_pstart, pp1, pp2, plower, pupper, plg);
  GetPrior(prior, dist_pp, p1, p2, lower, upper, lg);

  /*---------------------------------*/
  arma::mat lp(nmc, nchain); lp.fill(-arma::datum::inf); // sum log prior
  arma::mat ll(nmc, nchain); ll.fill(-arma::datum::inf); // sum log likelihood
  arma::cube theta(nchain, npar, nmc); theta.fill(NA_REAL);

  double tmp0, tmp1;
  for (size_t i = 0; i < nchain; i++) {
    unsigned int j = 1;

    while (lp.row(0).col(i).has_inf() || ll.row(0).col(i).has_inf()) {
      arma::vec pvec = rprior_vec_(dist_pstart, pp1, pp2, plower, pupper);

      tmp0 = sumlogprior(pvec, dist_pp, p1, p2, lower, upper, lg);
      tmp1 = sumloglike_glm(pvec, type, X, Y);

      theta.slice(0).row(i) = arma::trans(pvec);
      lp.row(0).col(i) = tmp0;
      ll.row(0).col(i)  = tmp1;
      j++;
      if (j > 1e4) { stop("Fail to set up new samples."); }
    }
  }


  List out = List::create(
    Named("theta")            = theta,
    Named("summed_log_prior") = lp,
    Named("log_likelihoods")  = ll,
    Named("data")             = data,
    Named("p.prior")          = prior,
    Named("start")            = 1,
    Named("n.pars")           = npar,
    Named("p.names")          = pnames,
    Named("rp")               = rp,
    Named("nmc")              = nmc,
    Named("thin")             = thin,
    Named("n.chains")         = nchain);
  return out;
}


arma::mat init_logit_test(unsigned int nmc, List data, List start, List prior,
                double rp, unsigned int thin, unsigned int nchain) {
  // std::vector<std::string> a = data.attr("names");
  // std::vector<std::string> b = {"R", "N", "Y"};
  // Xnames   <- datnames[!datnames %in% c("R", "N", "Y")]
  StringVector a = data.attr("names");
  StringVector b = {"R", "N", "Y"};
  LogicalVector c = Rf_match(b, a, 0); // F F T T T
  StringVector d = a[!c];
  unsigned int xsize = d.size();

  arma::vec Y_ = data["Y"];
  arma::vec N_ = data["N"];
  arma::mat Y = arma::join_horiz(Y_, N_);
  arma::mat X(Y_.n_elem, xsize + 2);

  arma::vec tmp1 = data["S"];
  arma::vec tmp2 = data["R"];
  X.col(0).zeros();
  X.col(1) = tmp1;
  X.col(2) = tmp2;
  X.col(3) = tmp1 % tmp2;

  // return a[!c];
  return X;
}

// [[Rcpp::export]]
List init_logit(unsigned int nmc, List data, List start, List prior,
                    double rp, unsigned int thin,
                    unsigned int nchain) {


  NumericVector modelAttr = data.attr("model"); // To extract attributes
  std::string type = modelAttr.attr("type");    // model type: DDM or LBA

  arma::mat Y      = data.attr("Y");
  arma::mat X_     = data.attr("X"); // glm independent / predictor variable
  arma::mat X = arma::join_horiz(arma::ones(X_.n_rows), X_);


  unsigned int npar = prior.size();
  std::vector<std::string> pnames = prior.attr("names");
  arma::vec pp1, pp2, plower, pupper, dist_pstart, p1, p2, lower, upper, dist_pp;
  arma::uvec plg, lg;
  GetPrior(start, dist_pstart, pp1, pp2, plower, pupper, plg);
  GetPrior(prior, dist_pp, p1, p2, lower, upper, lg);

  arma::mat lp(nmc, nchain); lp.fill(R_NegInf);   // sum log prior
  arma::mat ll(nmc, nchain); ll.fill(R_NegInf);   // sum log likelihood
  arma::cube theta(nchain, npar, nmc); theta.fill(NA_REAL);
  arma::uvec pidx_non = arma::find_nonfinite(dist_pp);

  double tmp0, tmp1;
  for (size_t i = 0; i < nchain; i++) {
    unsigned int j = 1;
    while (lp.row(0).col(i).has_inf() || ll.row(0).col(i).has_inf()) {
      arma::vec pvec = rprior_vec_(dist_pstart, pp1, pp2, plower, pupper);

      tmp0 = sumlogprior(pvec, dist_pp, pp1, pp2, plower, pupper, plg);
      tmp1 = sumloglike_glm(pvec, type, X, Y);

      theta.slice(0).row(i) = arma::trans(pvec);
      lp.row(0).col(i) = tmp0;
      ll.row(0).col(i) = tmp1;
      j++;
      if (j > 1e4) { stop("Fail to set up new samples."); }
    }
  }

  List out = List::create(
    Named("theta")            = theta,
    Named("summed_log_prior") = lp,
    Named("log_likelihoods")  = ll,
    Named("data")         = data,
    Named("p.prior")      = prior,
    Named("start")        = 1,
    Named("n.pars")       = npar,
    Named("p.names")      = pnames,
    Named("rp")           = rp,
    Named("nmc")          = nmc,
    Named("thin")         = thin,
    Named("n.chains")     = nchain);
  return out;
}

// [[Rcpp::export]]
List init_logits(unsigned int nmc, List data, List start, List prior,
                 double rp, unsigned int thin,
                 unsigned int nchain) {

  NumericVector modelAttr = data.attr("model"); // To extract attributes
  std::string type = modelAttr.attr("type");    // model type: DDM or LBA

  arma::mat Y      = data.attr("Y");
  arma::mat X_     = data.attr("X"); // glm independent / predictor variable
  arma::mat X = arma::join_horiz(arma::ones(X_.n_rows), X_);


  unsigned int npar = prior.size();
  std::vector<std::string> pnames = prior.attr("names");
  arma::vec pp1, pp2, plower, pupper, dist_pstart, p1, p2, lower, upper, dist_pp;
  arma::uvec plg, lg;
  GetPrior(start, dist_pstart, pp1, pp2, plower, pupper, plg);
  GetPrior(prior, dist_pp, p1, p2, lower, upper, lg);

  arma::mat lp(nmc, nchain); lp.fill(R_NegInf);   // sum log prior
  arma::mat ll(nmc, nchain); ll.fill(R_NegInf);   // sum log likelihood
  arma::cube theta(nchain, npar, nmc); theta.fill(NA_REAL);
  arma::uvec pidx_non = arma::find_nonfinite(dist_pp);

  double tmp0, tmp1;
  for (size_t i = 0; i < nchain; i++) {
    unsigned int j = 1;
    while (lp.row(0).col(i).has_inf() || ll.row(0).col(i).has_inf()) {
      arma::vec pvec = rprior_vec_(dist_pstart, pp1, pp2, plower, pupper);

      tmp0 = sumlogprior(pvec, dist_pp, pp1, pp2, plower, pupper, plg);
      tmp1 = sumloglike_glm(pvec, type, X, Y);

      theta.slice(0).row(i) = arma::trans(pvec);
      lp.row(0).col(i) = tmp0;
      ll.row(0).col(i) = tmp1;
      j++;
      if (j > 1e4) { stop("Fail to set up new samples."); }
    }
  }

  List out = List::create(
    Named("theta")            = theta,
    Named("summed_log_prior") = lp,
    Named("log_likelihoods")  = ll,
    Named("data")         = data,
    Named("p.prior")      = prior,
    Named("start")        = 1,
    Named("n.pars")       = npar,
    Named("p.names")      = pnames,
    Named("rp")           = rp,
    Named("nmc")          = nmc,
    Named("thin")         = thin,
    Named("n.chains")     = nchain);
  return out;
}

// [[Rcpp::export]]
void init_hlogit(unsigned int nmc, List data, List start, List prior,
                 double rp, unsigned int thin, unsigned int nchain) {
  List pstart = start[0];  /* Extract pprior & ppprior */
  List lstart = start[1];
  List sstart = start[2];
  List pprior = prior[0];
  List lprior = prior[1];
  List sprior = prior[2];
  unsigned int npar = pprior.size();
  unsigned int nsub = data.size();
  std::vector<std::string> pnames = pprior.attr("names");

  /* In hierarchical model, p.prior may carry NA, indicating the parameters are
   drawn from an upper level distribution.  Therefore in the case of
   multi-level modeling, we must supply 'pstart' to init_newnonhier_start */
  // List out = init_logits(nmc, data, pstart, pprior, rp, thin, nchain);


}


List init_new_glm_bk(unsigned int nmc, List prior, List data, double rp,
                  unsigned int thin, unsigned int nchain) {
  arma::ucube model = data.attr("model"); // do not have attributes
  arma::vec Y       = data["Y"];
  arma::vec X_vec   = data["X"]; // glm independent / predictor variable
  arma::mat X = arma::join_horiz(arma::ones(X_vec.n_elem), X_vec);

  unsigned int npar = prior.size();
  std::vector<std::string> pnames = prior.attr("names");

  // extract data-model options
  NumericVector modelAttr = data.attr("model"); // To extract attributes
  std::string type = modelAttr.attr("type");    // model type: DDM or LBA
  arma::vec allpar = modelAttr.attr("all.par");
  List modelDim    = modelAttr.attr("dimnames");
  arma::umat n1idx = modelAttr.attr("n1.order");
  bool posdrift    = modelAttr.attr("posdrift"); // do not have attributes

  arma::uvec matchCell   = modelAttr.attr("match.cell");
  arma::uvec isCellEmpty = data.attr("cell.empty") ;
  arma::umat cellidx = cellIdx2Mat(data);
  std::vector<std::string> parnames = modelAttr.attr("par.names");
  std::vector<std::string> dim1 = modelDim[0] ; // row
  std::vector<std::string> dim2 = modelDim[1] ; // col; para
  std::vector<std::string> dim3 = modelDim[2] ; // slice; r1, r2
  arma::uvec isr1 = GetIsR1(modelAttr, type);

  /*---------------------------------*/
  // should past p1 and p2 value from hierarchical level, if called by init_newhier
  arma::vec p1(npar), p2(npar), lower(npar), upper(npar), dists(npar);
  arma::uvec islog(npar);
  GetPrior(prior, dists, p1, p2, lower, upper, islog);

  /*---------------------------------*/
  arma::mat slp(nmc, nchain); slp.fill(-arma::datum::inf); // sum log prior
  arma::mat ll(nmc, nchain); ll.fill(-arma::datum::inf); // log likelihood
  arma::cube theta(nchain, npar, nmc); theta.fill(NA_REAL);

  double tmp0, tmp1;

  for (size_t i = 0; i < nchain; i++) {
    unsigned int j = 1;

    while (slp.row(0).col(i).has_inf() || ll.row(0).col(i).has_inf()) {
      arma::vec pvec = rprior_vec_(dists, p1, p2, lower, upper);

      tmp0 = sumlogprior(pvec, dists, p1, p2, lower, upper, islog);
      tmp1 = sumloglike_glm(pvec, type, X, Y);

      theta.slice(0).row(i) = arma::trans(pvec);
      slp.row(0).col(i) = tmp0;
      ll.row(0).col(i)  = tmp1;
      j++;
      if (j > 1e4) {
        stop("Fail to set up new samples.");
      }
    }
  }


  List out = List::create(
    Named("theta")            = theta,
    Named("summed_log_prior") = slp,
    Named("log_likelihoods")  = ll,
    Named("data")             = data,
    Named("p.prior")          = prior,
    Named("start")            = 1,
    Named("n.pars")           = npar,
    Named("p.names")          = pnames,
    Named("rp")               = rp,
    Named("nmc")              = nmc,
    Named("thin")             = thin,
    Named("n.chains")         = nchain);
  return out;
}

// [[Rcpp::export]]
List init_add_glm(unsigned int nmc, List samples, double rp,
              unsigned int thin) {
  List samples_in(clone(samples));
  arma::cube theta = samples_in["theta"];
  arma::mat slp    = samples_in["summed_log_prior"];
  arma::mat ll     = samples_in["log_likelihoods"];
  unsigned int pnmc = samples_in["nmc"]; // previous nmc
  unsigned int npar = samples_in["n.pars"];
  unsigned int nchain    = samples_in["n.chains"];
  unsigned int newnmc    = nmc + pnmc;

  Rcout << "Add " << nmc << " new samples.\n";
  arma::cube newtheta = arma::resize(theta, nchain, npar, newnmc);
  arma::mat newlp     = arma::resize(slp, newnmc, nchain);
  arma::mat newll     = arma::resize(ll,  newnmc, nchain);

  newtheta.slices(pnmc, newnmc - 1).fill(NA_REAL);
  newlp.rows(pnmc, newnmc - 1).fill(-arma::datum::inf);
  newll.rows(pnmc, newnmc - 1).fill(-arma::datum::inf);

  List samples_out = List::create(
    Rcpp::Named("theta")            = newtheta,
    Rcpp::Named("summed_log_prior") = newlp,
    Rcpp::Named("log_likelihoods")  = newll,
    Rcpp::Named("data")             = samples_in["data"],
    Rcpp::Named("p.prior")          = samples_in["p.prior"],
    Rcpp::Named("start")            = pnmc,
    Rcpp::Named("n.pars")           = npar,
    Rcpp::Named("p.names")          = samples_in["p.names"],
    Rcpp::Named("rp")               = rp,
    Rcpp::Named("nmc")              = newnmc,
    Rcpp::Named("thin")             = thin,
    Rcpp::Named("n.chains")         = nchain);
  return samples_out;
}

// [[Rcpp::export]]
List init_old_glm(unsigned int nmc, List samples, double rp,
              unsigned int thin) {
  List samples_in(clone(samples));
  arma::cube theta = samples_in["theta"];
  arma::mat slp    = samples_in["summed_log_prior"];
  arma::mat ll     = samples_in["log_likelihoods"];
  unsigned int pnmc      = samples_in["nmc"];
  unsigned int npar      = samples_in["n.pars"];
  unsigned int nchain    = samples_in["n.chains"];
  unsigned int startfrom = pnmc - 1;

  // arma::cube newtheta = arma::resize(theta, nchain, npar, nmc);
  // arma::mat newlp     = arma::resize(slp, nmc, nchain);
  // arma::mat newll     = arma::resize(ll,  nmc, nchain);

  arma::cube newtheta(nchain, npar, nmc);
  arma::mat newlp(nmc, nchain);
  arma::mat newll(nmc, nchain);

  newtheta.fill(NA_REAL);
  newlp.fill(R_NegInf);
  newll.fill(R_NegInf);

  newtheta.slice(0) = theta.slice(startfrom);
  newlp.row(0)      = slp.row(startfrom);
  newll.row(0)      = ll.row(startfrom);

  List samples_out = List::create(
    Rcpp::Named("theta")            = newtheta,
    Rcpp::Named("summed_log_prior") = newlp,
    Rcpp::Named("log_likelihoods")  = newll,
    Rcpp::Named("data")             = samples_in["data"],
    Rcpp::Named("p.prior")          = samples_in["p.prior"],
    Rcpp::Named("start")            = 1,
    Rcpp::Named("n.pars")           = npar,
    Rcpp::Named("p.names")          = samples_in["p.names"],
    Rcpp::Named("rp")               = rp,
    Rcpp::Named("nmc")              = nmc,
    Rcpp::Named("thin")             = thin,
    Rcpp::Named("n.chains")         = nchain);
  return samples_out;
}


// [[Rcpp::export]]
List init_newnonhier_glm(unsigned int nmc, List data, List pprior,
                     double rp, unsigned int thin, unsigned int nchain) {

  unsigned int nsub = data.size();
  List out(nsub);

  for (size_t i = 0; i < nsub; i++) {
    out[i] = init_new_glm_bk(nmc, pprior, data[i], rp, thin, nchain);
  }

  out.attr("names") = data.attr("names");
  return out;
}

// [[Rcpp::export]]
List init_addnonhier_glm(unsigned int nmc, List samples, double rp,
                         unsigned int thin) {

  unsigned int nsub = samples.size();
  List out(nsub);

  for (size_t i = 0; i < nsub; i++) {
    out[i] = init_add_glm(nmc, samples[i], rp, thin);
  }
  out.attr("names") = samples.attr("names");
  return out;
}

// [[Rcpp::export]]
List init_oldnonhier_glm(unsigned int nmc, List samples, double rp,
                     unsigned int thin) {

  unsigned int nsub = samples.size();
  List out(nsub);

  for (size_t i = 0; i < nsub; i++) {
    out[i] = init_old_glm(nmc, samples[i], rp, thin);
  }

  out.attr("names") = samples.attr("names");
  return out;
}
// [[Rcpp::export]]
List init_oldhier_glm(unsigned int nmc, List samples, double rp,
                  unsigned int thin) {

  List samples_in(clone(samples));
  List hyper    = samples_in.attr("hyper");
  List subject0 = samples_in[0];
  List pprior   = subject0["p.prior"];

  List phi = hyper["phi"];
  arma::mat hslp = hyper["h_summed_log_prior"];
  arma::mat hll  = hyper["h_log_likelihoods"];
  unsigned int pnmc      = hyper["nmc"];
  unsigned int npar      = subject0["n.pars"];
  unsigned int nhpar     = hyper["n.pars"];
  unsigned int nchain    = hyper["n.chains"];
  unsigned int startfrom = pnmc - 1;

  arma::cube location = phi[0];
  arma::cube scale    = phi[1];

  List data0 = subject0["data"];
  NumericVector modelAttr = data0.attr("model"); // To extract attributes
  NumericVector pvec = modelAttr.attr("p.vector");
  std::vector<std::string> pnames = pvec.attr("names");

  arma::cube newlocation(nchain, npar, nmc);
  arma::cube newscale(nchain, npar, nmc);
  arma::mat newhlp(nmc, nchain);
  arma::mat newhll(nmc, nchain);

  newlocation.fill(NA_REAL); newscale.fill(NA_REAL);
  newhlp.fill(R_NegInf); newhll.fill(R_NegInf);
  newlocation.slice(0) = location.slice(startfrom - 1);
  newscale.slice(0)    = scale.slice(startfrom - 1);
  newhlp.row(0)        = hslp.row(startfrom - 1);
  newhll.row(0)        = hll.row(startfrom - 1);

  List out = init_oldnonhier_glm(nmc, samples_in, rp, thin);

  // Finish up
  phi(0) = newlocation;
  phi(1) = newscale;
  List newhyper = List::create(   // 16 elements
    Rcpp::Named("phi")                = phi,
    Rcpp::Named("h_summed_log_prior") = newhlp,
    Rcpp::Named("h_log_likelihoods")  = newhll,
    Rcpp::Named("pp.prior") = hyper["pp.prior"],
    Rcpp::Named("start")    = 1,
    Rcpp::Named("n.pars")   = nhpar,
    Rcpp::Named("p.names")  = pnames, // pnames are from p.prior
    Rcpp::Named("rp")       = rp,
    Rcpp::Named("nmc")      = nmc,
    Rcpp::Named("thin")     = thin,
    Rcpp::Named("n.chains") = nchain);
  out.attr("hyper") = newhyper;
  return out;
}

// [[Rcpp::export]]
List init_addhier_glm(unsigned int nmc, List samples, double rp,
                  unsigned int thin) {

  List samples_in(clone(samples));
  List hyper = samples_in.attr("hyper");
  List subject0 = samples_in[0];
  List pprior   = subject0["p.prior"];
  std::vector<std::string> pnames = pprior.attr("names");

  List phi = hyper["phi"];
  arma::mat hslp = hyper["h_summed_log_prior"];
  arma::mat hll  = hyper["h_log_likelihoods"];
  unsigned int pnmc      = hyper["nmc"]; // previous nmc
  unsigned int npar      = subject0["n.pars"];
  unsigned int nhpar     = hyper["n.pars"];
  unsigned int nchain    = hyper["n.chains"];
  unsigned int newnmc    = nmc + pnmc;

  arma::cube location = phi[0];
  arma::cube scale    = phi[1];

  Rcout << "Add " << nmc << " new samples onto the existing one.\n";
  arma::cube newlocation(nchain, npar, newnmc);
  arma::cube newscale(nchain, npar, newnmc);
  arma::mat newhlp(newnmc, nchain);
  arma::mat newhll(newnmc, nchain);

  newlocation.slices(pnmc, newnmc - 1).fill(NA_REAL);
  newscale.slices(pnmc, newnmc - 1).fill(NA_REAL);
  newhlp.rows(pnmc, newnmc - 1).fill(R_NegInf);
  newhll.rows(pnmc, newnmc - 1).fill(R_NegInf);

  List out = init_addnonhier_glm(nmc, samples_in, rp, thin);
  // Finish up
  phi(0) = newlocation;
  phi(1) = newscale;
  List newhyper = List::create(   // 16 elements
    Rcpp::Named("phi") = phi,
    Rcpp::Named("h_summed_log_prior") = newhlp,
    Rcpp::Named("h_log_likelihoods")  = newhll,
    Rcpp::Named("pp.prior") = hyper["pp.prior"],
    Rcpp::Named("start")    = pnmc,
    Rcpp::Named("n.pars")   = nhpar,
    Rcpp::Named("p.names")  = pnames, // pnames are from p.prior
    Rcpp::Named("rp")       = rp,
    Rcpp::Named("nmc")      = newnmc,
    Rcpp::Named("thin")     = thin,
    Rcpp::Named("n.chains") = nchain);
  out.attr("hyper") = newhyper;
  return out;
}

// [[Rcpp::export]]
List init_new_start(unsigned int nmc, List data, List start, List prior,
                    double rp, unsigned int thin, unsigned int nchain) {

  arma::ucube model = data.attr("model"); // do not have attributes
  arma::vec Y_      = data["Y"];
  arma::mat Y       = Y_;
  arma::vec X_vec   = data["X"]; // glm independent / predictor variable
  arma::mat X = arma::join_horiz(arma::ones(X_vec.n_elem), X_vec);
  unsigned int npar = prior.size();
  std::vector<std::string> pnames = prior.attr("names");

  // extract data-model options
  NumericVector modelAttr = data.attr("model"); // To extract attributes
  std::string type = modelAttr.attr("type");    // model type: DDM or LBA
  arma::vec allpar = modelAttr.attr("all.par");
  List modelDim    = modelAttr.attr("dimnames");
  arma::umat n1idx = modelAttr.attr("n1.order");
  bool posdrift    = modelAttr.attr("posdrift"); // do not have attributes
  arma::uvec matchCell   = modelAttr.attr("match.cell");
  arma::uvec isCellEmpty = data.attr("cell.empty") ;
  arma::umat cellidx = cellIdx2Mat(data);
  std::vector<std::string> parnames = modelAttr.attr("par.names");
  std::vector<std::string> dim1 = modelDim[0] ; // row
  std::vector<std::string> dim2 = modelDim[1] ; // col; para
  std::vector<std::string> dim3 = modelDim[2] ; // slice; r1, r2
  arma::uvec isr1 = GetIsR1(modelAttr, type);

  /*---------------------------------*/
  arma::vec pp1, pp2, plower, pupper, dist_pstart, p1, p2, lower, upper, dist_pp;
  arma::uvec plg, lg;
  GetPrior(start, dist_pstart, pp1, pp2, plower, pupper, plg);
  GetPrior(prior, dist_pp, p1, p2, lower, upper, lg);

  arma::mat lp(nmc, nchain); lp.fill(R_NegInf); // sum log prior
  arma::mat ll(nmc, nchain); ll.fill(R_NegInf);   // log likelihood
  arma::cube theta(nchain, npar, nmc); theta.fill(NA_REAL);
  /*---------------------------------*/
  arma::uvec pidx_non = arma::find_nonfinite(dist_pp);

  double tmp0, tmp1;

  for (size_t i = 0; i < nchain; i++) {
    unsigned int j = 1;

    while (lp.row(0).col(i).has_inf() || ll.row(0).col(i).has_inf()) {
      arma::vec pvec = rprior_vec_(dist_pstart, pp1, pp2, plower, pupper);
      // pp1(pidx_non) = p1(pidx_non);
      // pp2(pidx_non) = p2(pidx_non);
      // plower(pidx_non) = lower(pidx_non);
      // pupper(pidx_non) = upper(pidx_non);

      tmp0 = sumlogprior(pvec, dist_pp, pp1, pp2, plower, pupper, plg);
      tmp1 = sumloglike_glm(pvec, type, X, Y);

      theta.slice(0).row(i) = arma::trans(pvec);
      lp.row(0).col(i) = tmp0;
      ll.row(0).col(i) = tmp1;
      j++;
      if (j > 1e4) { stop("Fail to set up new samples."); }
    }
  }

  List out = List::create(
    Named("theta")            = theta,
    Named("summed_log_prior") = lp,
    Named("log_likelihoods")  = ll,
    Named("data")             = data,
    Named("p.prior")          = prior,
    Named("start")            = 1,
    Named("n.pars")           = npar,
    Named("p.names")          = pnames,
    Named("rp")               = rp,
    Named("nmc")              = nmc,
    Named("thin")             = thin,
    Named("n.chains")         = nchain);
  return out;
}



// [[Rcpp::export]]
List init_newnonhier_start(unsigned int nmc, List data, List start, List prior,
                            double rp, unsigned int thin, unsigned int nchain) {

  unsigned int nsub = data.size();
  List out(nsub);

  for (size_t i = 0; i < nsub; i++) {
    out[i] = init_new_start(nmc, data[i], start, prior, rp, thin, nchain);
  }

  out.attr("names") = data.attr("names");
  return out;
}

// [[Rcpp::export]]
List init_newhier_start(unsigned int nmc, List data, List start, List prior,
                         double rp, unsigned int thin, unsigned int nchain) {

  List pstart = start[0];  /* Extract pprior & ppprior */
  List lstart = start[1];
  List sstart = start[2];
  List pprior = prior[0];
  List lprior = prior[1];
  List sprior = prior[2];
  unsigned int npar = pprior.size();
  unsigned int nsub = data.size();
  std::vector<std::string> pnames = pprior.attr("names");

  /* In hierarchical model, p.prior may carry NA, indicating the parameters are
     drawn from an upper level distribution.  Therefore in the case of
     multi-level modeling, we must supply 'pstart' to init_newnonhier_start */
  List out = init_newnonhier_start(nmc, data, pstart, pprior, rp, thin, nchain);

  arma::cube theta = GetTheta0(out); // thetas: nsub x npar x nchain

  arma::vec dist_lstart, dist_sstart, dist_pstart;
  arma::uvec llg, slg, plg;
  arma::vec pp1, pp2, lp1, lp2, sp1, sp2;
  arma::vec lupper, llower, supper, slower, pupper, plower;
  GetPrior(lstart, dist_lstart, lp1, lp2, llower, lupper, llg);
  GetPrior(sstart, dist_sstart, sp1, sp2, slower, supper, slg);
  GetPrior(pstart, dist_pstart, pp1, pp2, plower, pupper, plg);

  arma::vec p1, p2, lower, upper,
        loc_p1, loc_p2, loc_l, loc_u,
        sca_p1, sca_p2, sca_l, sca_u;
  arma::vec dist_pp, dist_lp, dist_sp;
  arma::uvec lg, loc_lg, sca_lg;
  GetPrior(pprior, dist_pp, p1, p2, lower, upper, lg);
  GetPrior(lprior, dist_lp, loc_p1, loc_p2, loc_l, loc_u, loc_lg);
  GetPrior(sprior, dist_sp, sca_p1, sca_p2, sca_l, sca_u, sca_lg);

  arma::uvec pidx_fin, lidx_fin, sidx_fin, pidx_non, lidx_non, sidx_non;

  pidx_fin = arma::find_finite(dist_pp);
  lidx_fin = arma::find_finite(dist_lp);
  sidx_fin = arma::find_finite(dist_sp);
  pidx_non = arma::find_nonfinite(dist_pp);
  lidx_non = arma::find_nonfinite(dist_lp);
  sidx_non = arma::find_nonfinite(dist_sp);

  unsigned int nppar = pidx_fin.n_elem;
  unsigned int nlpar = lidx_fin.n_elem;
  unsigned int nspar = sidx_fin.n_elem;

  // if (nhpar < npar) {
  //   idx0 = arma::find_nonfinite(loc_p1);
  //   for(size_t i = 0; i < idx0.n_elem; i++) {
  //     Rcout << pnames[idx0(i)] << " is " << loc_p1(idx0(i)) << std::endl;
  //   }
  // }

  arma::cube location(nchain, npar, nmc); location.fill(NA_REAL);
  arma::cube scale(nchain, npar, nmc); scale.fill(NA_REAL);
  arma::mat hlp(nmc, nchain); hlp.fill(R_NegInf); // hyper sum log prior
  arma::mat hll(nmc, nchain); hll.fill(R_NegInf);   // hyper log likelihood

  arma::vec lvec, svec;
  double loclik, scalik;

  for (size_t i = 0; i < nchain; i++) {
    lvec = rprior_vec_(dist_lstart(lidx_fin), lp1(lidx_fin), lp2(lidx_fin),
                       llower(lidx_fin), lupper(lidx_fin));
    svec = rprior_vec_(dist_sstart(sidx_fin), sp1(sidx_fin), sp2(sidx_fin),
                       slower(sidx_fin), supper(sidx_fin));
    // sum log h prior
    double loclik = sumlogprior(lvec, dist_lp(lidx_fin), loc_p1(lidx_fin),
                                loc_p2(lidx_fin), loc_l(lidx_fin),
                                loc_u(lidx_fin), loc_lg(lidx_fin));
    double scalik = sumlogprior(svec, dist_sp(sidx_fin), sca_p1(sidx_fin),
                                sca_p2(sidx_fin), sca_l(sidx_fin),
                                sca_u(sidx_fin), sca_lg(sidx_fin));

    // first name and i chain
    hlp.row(0).col(i) = loclik + scalik;

    if (npar != nlpar && dist_lp.has_nan()) {
      lvec.resize(npar);
      lvec.elem( lidx_non ).fill(NA_REAL);
    }

    if (npar != nspar && dist_sp.has_nan()) {
      svec.resize(npar);
      svec.elem( sidx_non ).fill(NA_REAL);
    }

    location.slice(0).row(i) = lvec.t();
    scale.slice(0).row(i)    = svec.t();


    hll.row(0).col(i) = sumloghlike(theta.slice(i).cols(pidx_fin),
            dist_pp(pidx_fin), lvec(pidx_fin), svec(pidx_fin),
            lower(pidx_fin), upper(pidx_fin), lg(pidx_fin));
  } // end of looping chains

  if (hll.row(0).has_inf()) {
    Rcout << "\nSampling start points were valid at data level but not at "
          << "hyper level. Narrower (hyper) prior distributions may help";
  }

  // Finish up
  Rcout << std::endl;
  arma::field<arma::cube> phi(2);
  phi(0) = location;
  phi(1) = scale;

  List ppprior = List ::create(
    _["location"] = lprior,
    _["scale"] = sprior);

  List hyper = List::create(   // 16 elements
    Rcpp::Named("phi") = phi,
    Rcpp::Named("h_summed_log_prior") = hlp,
    Rcpp::Named("h_log_likelihoods")  = hll,
    Rcpp::Named("pp.prior") = ppprior,
    Rcpp::Named("start")    = 1,
    Rcpp::Named("n.pars")   = npar,
    Rcpp::Named("p.names")  = pnames,
    Rcpp::Named("rp")       = rp,
    Rcpp::Named("nmc")      = nmc,
    Rcpp::Named("thin")     = thin,
    Rcpp::Named("n.chains") = nchain);
  out.attr("hyper") = hyper;
  return out;
}
