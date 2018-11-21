#include <ggdmc.hpp>
using namespace Rcpp;

// [[Rcpp::export]]
List init_new(unsigned int nmc, List pprior, List data, double rp,
  unsigned int thin, unsigned int nchain, unsigned int ncore = 1,
  bool debug = false) {
  arma::ucube model  = data.attr("model"); // do not have attributes
  unsigned int npda  = data.attr("n.pda");
  double bw          = data.attr("bw");
  unsigned int gpuid = data.attr("gpuid");
  arma::vec RT       = data["RT"];
  unsigned int npar = pprior.size();
  std::vector<std::string> pnames = pprior.attr("names");

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
  std::vector<std::string> dists(npar);
  arma::vec p1(npar), p2(npar), lower(npar), upper(npar);
  arma::uvec islog(npar);
  GetPrior(pprior, dists, p1, p2, lower, upper, islog);

  /*---------------------------------*/
  arma::mat slp(nmc, nchain); slp.fill(-INFINITY); // sum log prior
  arma::mat ll(nmc, nchain); ll.fill(-INFINITY); // log likelihood
  arma::cube theta(nchain, npar, nmc); theta.fill(NA_REAL);
  for (size_t i = 0; i < nchain; i++) {
    unsigned int j = 1;
    while (slp.row(0).col(i).has_inf() || ll.row(0).col(i).has_inf()) {
      arma::vec pvec = rprior_vec(dists, p1, p2, lower, upper);

      theta.slice(0).row(i) = arma::trans(pvec);
      slp.row(0).col(i) = sumlogprior(pvec, dists, p1, p2, lower, upper, islog);
      ll.row(0).col(i)  = sumloglike(pvec, pnames, allpar, parnames, model,
        type, dim1, dim2, dim3, n1idx, isCellEmpty, cellidx, RT, matchCell,
        isr1, posdrift, npda, bw, ncore, gpuid, debug);
      j++;
      if (j > 1e4) { stop("Fail to set up new samples.");}
    }
  }


  List out = List::create(
    Named("theta")            = theta,
    Named("summed_log_prior") = slp,
    Named("log_likelihoods")  = ll,
    Named("data")             = data,
    Named("p.prior")          = pprior,
    Named("start")            = 1,
    Named("n.pars")           = npar,
    Named("p.names")          = pnames,
    Named("rp")               = rp,
    Named("nmc")              = nmc,
    Named("thin")             = thin,
    Named("n.chains")         = nchain,
    Named("npda")             = npda,
    Named("gpuid")            = gpuid,
    Named("bw")               = bw,
    Named("ncore")            = ncore);
  return out;
}

// [[Rcpp::export]]
List init_old(unsigned int nmc, List samples, double rp,
  unsigned int thin) {
  List samples_in(clone(samples));
  arma::cube theta = samples_in["theta"];
  arma::mat slp    = samples_in["summed_log_prior"];
  arma::mat ll     = samples_in["log_likelihoods"];
  unsigned int pnmc      = samples_in["nmc"];
  unsigned int npar      = samples_in["n.pars"];
  unsigned int nchain    = samples_in["n.chains"];
  unsigned int startfrom = pnmc - 1;

  arma::cube newtheta = arma::resize(theta, nchain, npar, nmc);
  arma::mat newlp     = arma::resize(slp, nmc, nchain);
  arma::mat newll     = arma::resize(ll,  nmc, nchain);
  newtheta.fill(NA_REAL); newlp.fill(-INFINITY); newll.fill(-INFINITY);

  newtheta.slice(0)   = theta.slice(startfrom);
  newlp.row(0)        = slp.row(startfrom);
  newll.row(0)        = ll.row(startfrom);

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
List init_add(unsigned int nmc, List samples, double rp,
  unsigned int thin) {
    List samples_in(clone(samples));
    arma::cube theta = samples_in["theta"];
    arma::mat slp    = samples_in["summed_log_prior"];
    arma::mat ll     = samples_in["log_likelihoods"];
    unsigned int pnmc = samples_in["nmc"]; // previous nmc
    unsigned int npar = samples_in["n.pars"];
    unsigned int nchain    = samples_in["n.chains"];
    unsigned int newnmc    = nmc + pnmc;

    Rcout << "Add " << nmc << " new samples onto the existing one.\n";
    arma::cube newtheta = arma::resize(theta, nchain, npar, newnmc);
    arma::mat newlp     = arma::resize(slp, newnmc, nchain);
    arma::mat newll     = arma::resize(ll,  newnmc, nchain);

    newtheta.slices(pnmc, newnmc - 1).fill(NA_REAL);
    newlp.rows(pnmc, newnmc - 1).fill(-INFINITY);
    newll.rows(pnmc, newnmc - 1).fill(-INFINITY);

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
List init_newnonhier(unsigned int nmc, List data, List pprior,
  double rp, unsigned int thin, unsigned int nchain) {

  unsigned int nsub = data.size();
  List out(nsub);

  for (size_t i = 0; i < nsub; i++) {
    out[i] = init_new(nmc, pprior, data[i], rp, thin, nchain);
  }

  out.attr("names") = data.attr("names");
  return out;
}

// [[Rcpp::export]]
List init_oldnonhier(unsigned int nmc, List samples, double rp,
  unsigned int thin) {

  unsigned int nsub = samples.size();
  List out(nsub);

  for (size_t i = 0; i < nsub; i++) {
    out[i] = init_old(nmc, samples[i], rp, thin);
  }

  out.attr("names") = samples.attr("names");
  return out;
}

// [[Rcpp::export]]
List init_addnonhier(unsigned int nmc, List samples, double rp,
  unsigned int thin) {

  unsigned int nsub = samples.size();
  List out(nsub);

  for (size_t i = 0; i < nsub; i++) {
    out[i] = init_add(nmc, samples[i], rp, thin);
  }
  out.attr("names") = samples.attr("names");
  return out;
}

// [[Rcpp::export]]
List init_newhier(unsigned int nmc, List data, List pprior, List ppprior,
  double rp,  unsigned int thin, unsigned int nchain) {

  List data0 = data[0];  // Temporary measure to save PDA options
  unsigned int npda  = data0.attr("n.pda");
  unsigned int gpuid = data0.attr("gpuid");
  double bw          = data0.attr("bw");
  unsigned int npar = pprior.size();
  unsigned int nsub = data.size();

  NumericVector modelAttr = data0.attr("model"); // To extract attributes
  // Rcout << modelAttr;

  NumericVector pvec = modelAttr.attr("p.vector");
  std::vector<std::string> pnames = pvec.attr("names");

  // for (size_t i = 0; i < pnames.size(); i++) {
  //   Rcout << pnames[i] << std::endl;
  // }

  arma::cube location(nchain, npar, nmc); location.fill(NA_REAL);
  arma::cube scale(nchain, npar, nmc); scale.fill(NA_REAL);

  List lprior = ppprior[0];   /* Extract pprior & ppprior */
  List sprior = ppprior[1];
  std::vector<std::string> pdists(npar), ldists(npar), sdists(npar);
  arma::vec pp1(npar), pp2(npar), plower(npar), pupper(npar),
            lp1(npar), lp2(npar), llower(npar), lupper(npar),
            sp1(npar), sp2(npar), slower(npar), supper(npar);
  arma::uvec plog(npar), llog(npar), slog(npar);

  GetPrior(pprior, pdists, pp1, pp2, plower, pupper, plog);
  GetPrior(lprior, ldists, lp1, lp2, llower, lupper, llog);
  GetPrior(sprior, sdists, sp1, sp2, slower, supper, slog);


  // Rcout <<" before init_newnon\n";

  // AH's lapply-style loop; cps subject list: nchain x npar; cps == thetas
  List out = init_newnonhier(nmc, data, pprior, rp, thin, nchain);
  arma::cube thetas = GetTheta0(out); // thetas: nsub x npar x nchain
  //
  // // Rcout <<" after init_newnon\n";
  arma::mat slp(nmc, nchain); slp.fill(-INFINITY);
  arma::mat hslp(nmc, nchain); hslp.fill(-INFINITY); // hyper sum log prior
  arma::mat hll(nmc, nchain); hll.fill(-INFINITY); // hyper log likelihood
  Rcout << "Generating hyper-start points for each chain: ";
  for (size_t i = 0; i < nchain; i++) {
    Rcout << ".";

     // Step 1: Generate random p1 == location; p2 == scale
    arma::vec lstart = rprior_vec(ldists, lp1, lp2, llower, lupper);
    arma::vec sstart = rprior_vec(sdists, sp1, sp2, slower, supper);
    location.slice(0).row(i) = arma::trans(lstart);
    scale.slice(0).row(i) = arma::trans(sstart);

     // Step 2: Calculate hyper-sum-log-prior and hyper-sum-log-likelihood
    hslp.row(0).col(i) = sumloghprior(lstart, sstart, ldists, sdists, lp1, sp1,
      lp2, sp2, llower, slower, lupper, supper, llog, slog);
    hll.row(0).col(i) = sumloghlike(thetas.slice(i), pdists, lstart, sstart,
      plower, pupper, plog);

     // Step 3: Reset sum-log_prior calculated in init_nonhier
     // sum-log-likelihood is based on data, so no reset
    for (size_t j = 0; j < nsub; j ++) {
      List subjecti = out[j];
      arma::mat tmp_slp = subjecti["summed_log_prior"];
      tmp_slp.row(0).col(i) = sumlogprior(arma::trans(thetas.slice(i).row(j)),
        pdists, lstart, sstart, plower, pupper, plog);
      subjecti["summed_log_prior"] = tmp_slp;
      out[j] = subjecti;
    }
  }

  if (hll.row(0).has_inf()) {
    Rcout << "\nSampling start points were valid at data level but not at hyper "
          << "level. Try a tighter pp.prior or hstart.prior";
  }
  Rcout << std::endl;   // Finish up
  arma::field<arma::cube> phi(2); phi(0) = location; phi(1) = scale;
  List hyper = List::create(   // 16 elements
    Rcpp::Named("phi") = phi,
    Rcpp::Named("h_summed_log_prior") = hslp,
    Rcpp::Named("h_log_likelihoods")  = hll,
    Rcpp::Named("pp.prior") = ppprior,
    Rcpp::Named("start")    = 1,
    Rcpp::Named("n.pars")   = npar,
    Rcpp::Named("p.names")  = pnames,
    Rcpp::Named("rp")       = rp,
    Rcpp::Named("nmc")      = nmc,
    Rcpp::Named("thin")     = thin,
    Rcpp::Named("n.chains") = nchain,
    Rcpp::Named("npda")     = npda,
    Rcpp::Named("gpuid")    = gpuid,
    Rcpp::Named("bw")       = bw,
    Rcpp::Named("ncore")    = 1);
  out.attr("hyper") = hyper;
  return out;
}

// [[Rcpp::export]]
List init_oldhier(unsigned int nmc, List samples, double rp,
  unsigned int thin) {

  List samples_in(clone(samples));
  List hyper = samples_in.attr("hyper");
  List subject0 = samples_in[0];
  List pprior   = subject0["p.prior"];

  List phi = hyper["phi"];
  arma::mat hslp = hyper["h_summed_log_prior"];
  arma::mat hll  = hyper["h_log_likelihoods"];
  unsigned int pnmc      = hyper["nmc"];
  unsigned int npar      = hyper["n.pars"];
  unsigned int nchain    = hyper["n.chains"];
  unsigned int startfrom = pnmc - 1;
  arma::cube location = phi[0];
  arma::cube scale    = phi[1];

  List data0 = subject0["data"];
  NumericVector modelAttr = data0.attr("model"); // To extract attributes
  NumericVector pvec = modelAttr.attr("p.vector");
  std::vector<std::string> pnames = pvec.attr("names");

  arma::cube newlocation = arma::resize(location, nchain, npar, nmc);
  arma::cube newscale    = arma::resize(scale, nchain, npar, nmc);
  arma::mat newhlp = arma::resize(hslp, nmc, nchain);
  arma::mat newhll = arma::resize(hll, nmc, nchain);

  newlocation.fill(NA_REAL); newscale.fill(NA_REAL);
  newhlp.fill(-INFINITY); newhll.fill(-INFINITY);
  newlocation.slice(0) = location.slice(startfrom - 1);
  newscale.slice(0)    = scale.slice(startfrom - 1);
  newhlp.row(0)        = hslp.row(startfrom - 1);
  newhll.row(0)        = hll.row(startfrom - 1);

  List out = init_oldnonhier(nmc, samples_in, rp, thin);

  // Finish up
  phi(0) = newlocation; phi(1) = newscale;
  List newhyper = List::create(   // 16 elements
    Rcpp::Named("phi") = phi,
    Rcpp::Named("h_summed_log_prior") = newhlp,
    Rcpp::Named("h_log_likelihoods")  = newhll,
    Rcpp::Named("pp.prior") = hyper["pp.prior"],
    Rcpp::Named("start")    = 1,
    Rcpp::Named("n.pars")   = npar,
    Rcpp::Named("p.names")  = pnames, // pnames are from p.prior
    Rcpp::Named("rp")       = rp,
    Rcpp::Named("nmc")      = nmc,
    Rcpp::Named("thin")     = thin,
    Rcpp::Named("n.chains") = nchain,
    Rcpp::Named("npda")     = hyper["npda"],
    Rcpp::Named("gpuid")    = hyper["gpuid"],
    Rcpp::Named("bw")       = hyper["bw"],
    Rcpp::Named("ncore")    = 1);
  out.attr("hyper") = newhyper;
  return out;
}

// [[Rcpp::export]]
List init_addhier(unsigned int nmc, List samples, double rp,
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
  unsigned int npar      = hyper["n.pars"];
  unsigned int nchain    = hyper["n.chains"];
  unsigned int newnmc    = nmc + pnmc;

  arma::cube location = phi[0];
  arma::cube scale    = phi[1];

  Rcout << "Add " << nmc << " new samples onto the existing one.\n";
  arma::cube newlocation = arma::resize(location, nchain, npar, newnmc);
  arma::cube newscale    = arma::resize(scale, nchain, npar, newnmc);
  arma::mat newhlp = arma::resize(hslp, newnmc, nchain);
  arma::mat newhll = arma::resize(hll, newnmc, nchain);

  newlocation.slices(pnmc, newnmc - 1).fill(NA_REAL);
  newscale.slices(pnmc, newnmc - 1).fill(NA_REAL);
  newhlp.rows(pnmc, newnmc - 1).fill(-INFINITY);
  newhll.rows(pnmc, newnmc - 1).fill(-INFINITY);

  List out = init_addnonhier(nmc, samples_in, rp, thin);
  // Finish up
  phi(0) = newlocation; phi(1) = newscale;
  List newhyper = List::create(   // 16 elements
    Rcpp::Named("phi") = phi,
    Rcpp::Named("h_summed_log_prior") = newhlp,
    Rcpp::Named("h_log_likelihoods")  = newhll,
    Rcpp::Named("pp.prior") = hyper["pp.prior"],
    Rcpp::Named("start")    = pnmc,
    Rcpp::Named("n.pars")   = npar,
    Rcpp::Named("p.names")  = pnames, // pnames are from p.prior
    Rcpp::Named("rp")       = rp,
    Rcpp::Named("nmc")      = newnmc,
    Rcpp::Named("thin")     = thin,
    Rcpp::Named("n.chains") = nchain,
    Rcpp::Named("npda")     = hyper["npda"],
    Rcpp::Named("gpuid")    = hyper["gpuid"],
    Rcpp::Named("bw")       = hyper["bw"],
    Rcpp::Named("ncore")    = 1);
  out.attr("hyper") = newhyper;
  return out;
}
