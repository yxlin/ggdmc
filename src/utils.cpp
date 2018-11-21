//#include <RcppArmadillo.h>
#include <ggdmc.hpp>
using namespace Rcpp;

//' Whether a hyper-prior distribution is set constant
//'
//' Check if a prior distribution for a location and scale parameter is set
//' constant. \code{testHyper} checks if the parameter names in location and
//' in scale match.
//'
//' @param ppprior hyper-parameter prior distributions. First element is a
//' location prior and second is a scale prior.
//' @return \code{isConstant} gives a npar x 2 matrix;
//' @examples
//' model <- BuildModel(p.map = list(A = "1", B = "R", t0 = "1",
//'   mean_v = c("F", "M"), sd_v = "M", st0 = "1"),
//'   match.map = list(M = list(s1=1, s2=2)),
//'   factors   = list(S = c("s1", "s2"),F = c("f1", "f2")),
//'     constants = c(sd_v.false = 1, st0 = 0),
//'     responses = c("r1", "r2"),
//'     type      = "norm")
//'   npar <- length(GetPNames(model))
//'
//' ## Population distribution, rate effect on F
//'   pop.mean <- c(A=.4, B.r1=.6, B.r2=.8, t0=.3, mean_v.f1.true=1.5,
//'     mean_v.f2.true=1, mean_v.f1.false=0, mean_v.f2.false=0,
//'     sd_v.true = .25)
//'     pop.scale <-c(A=.1, B.r1=.1, B.r2=.1, t0=.05, mean_v.f1.true=.2,
//'       mean_v.f2.true=.2, mean_v.f1.false=.2, mean_v.f2.false=.2,
//'       sd_v.true = .1)
//'     p.prior <- BuildPrior(
//'         dists = rep("tnorm", npar),
//'         p1    = pop.mean,
//'         p2    = pop.scale*5,
//'         lower = c(0,0,0,.1,NA,NA,NA,NA,0),
//'         upper = c(NA,NA,NA,NA,NA,NA,NA,NA,NA))
//'     mu.prior <- BuildPrior(
//'         dists = rep("tnorm", npar),
//'         p1    = pop.mean,
//'         p2    = c(1,1,1,1,2,2,2,2,1),
//'         lower = c(0,0,0,.1,NA,NA,NA,NA,0),
//'         upper = c(NA,NA,NA,NA,NA,NA,NA,NA,NA))
//'     sigma.prior <- BuildPrior(
//'         dists = rep("beta", npar),
//'         p1    = c(A=1, B.r1=1, B.r2=1, t0=1, mean_v.f1.true=1,
//'           mean_v.f2.true=1, mean_v.f1.false=1, mean_v.f2.false=1,
//'           sd_v.true = 1),
//'           p2    = rep(1, npar))
//'     pp.prior <- list(mu.prior, sigma.prior)
//'
//'     test1 <- GetConstIdx(pp.prior)
//'     test2 <- lapply(pp.prior,function(x){
//'       lapply(x,function(y){attr(y,"dist") == "constant"})})
//' @export
// [[Rcpp::export]]
arma::umat GetConstIdx(List ppprior) {

  if (ppprior.size() != 2) stop("pp.prior must be a two-element list.");
  List muprior, sigmaprior;
  muprior    = ppprior[0];
  sigmaprior = ppprior[1];
  if (muprior.size() != sigmaprior.size()) {
    stop("mu.prior and sigma.prior have different numbers of elements.");
  }

  // first column is location and second is scale
  unsigned int npar = muprior.size();
  arma::umat out(npar, 2);
  for (size_t i = 0; i < npar; i++) {
    List location = muprior[i];
    List scale    = sigmaprior[i];
    std::string string1 = location.attr("dist");
    std::string string2 = scale.attr("dist");
    out(i, 0) = (string1 == "constant");
    out(i, 1) = (string2 == "constant");
  }
  return out;
}


// [[Rcpp::export]]
LogicalVector MatchPnames(List samples) {
  List data          = samples["data"];
  NumericVector modelAttr = data.attr("model");
  NumericVector pvec = modelAttr.attr("p.vector");
  List pprior = samples["p.prior"];

  CharacterVector pnames1 = pvec.attr("names");
  CharacterVector pnames2 = pprior.attr("names");
  unsigned int npar1 = pnames1.size();

  for(size_t i = 0; i < npar1; i++) {
    Rcout << "p.vector " << pnames1[i] << " p.prior" << pnames2[i] << std::endl;
  }

  // SEXP Rf_match(SEXP itable, SEXP ix, int nmatch)
  //   where nmatch is the value to be returned if there is no match
  LogicalVector res = Rf_match(pnames1, pnames2, 0);
  return res;
}


// [[Rcpp::export]]
LogicalVector MatchPPPriorName(List ppprior) {
  if (ppprior.size() != 2) stop("pp.prior must be a 2-element list.");
  List muprior    = ppprior[0];
  List sigmaprior = ppprior[1];
  if (muprior.size() != sigmaprior.size()) stop("mu.prior & sigma.prior are inconsistent.");
  CharacterVector munames = muprior.names();
  CharacterVector sigmanames = sigmaprior.names();
  return Rf_match(munames, sigmanames, 0);
}

void InitializeSubjectRJ(List samples, arma::field<arma::umat>& rj) {

  List hyper = samples.attr("hyper");
  unsigned int nmc   = hyper["nmc"];
  // unsigned int startR= hyper["start"]; // start_R == 1;
  unsigned int nchain= hyper["n.chains"];
  unsigned int nsub = samples.size();
  // unsigned int nsamp = 1 + (nmc - startR) * thin;

  for (size_t i = 0; i < nsub; i++) {
    rj(i) = arma::umat(nchain, nmc, arma::fill::zeros);
  }
}


void InitializeOneSubject(List samples, arma::umat& rj) {
  unsigned int nchain  = samples["n.chains"];
  unsigned int nmc     = samples["nmc"];
  rj = arma::umat(nchain, nmc, arma::fill::zeros);
}

void TransformSubjects(List samples,
  arma::field<arma::cube>& thetas,
  arma::field<arma::mat>& usethetas,
  arma::field<arma::mat>& logpriors,
  arma::field<arma::vec>& uselogpriors,
  arma::field<arma::mat>& loglikes,
  arma::field<arma::vec>& useloglikes,
  arma::uvec& store_i,
  std::vector<std::string>& types,
  arma::field<arma::vec>& allpars, arma::field<arma::umat>& n1idxes,
  arma::field<arma::uvec>& matchcells, arma::field<arma::uvec>& emptycells,
  arma::field<arma::umat>& cellidxes, arma::field<std::vector<std::string>>& parnames,
  arma::field<std::vector<std::string>>& dim1s,
  arma::field<std::vector<std::string>>& dim2s,
  arma::field<std::vector<std::string>>& dim3s,
  arma::field<arma::uvec>& isr1s, arma::uvec& posdrift,
  arma::field<arma::ucube>& models,
  arma::uvec& npdas, arma::vec& bws, arma::uvec& gpuids,
  arma::field<arma::vec>& RTs) {

  unsigned int nsub = samples.size();

  for (size_t i = 0; i < nsub; i++) {
    List subjecti      = samples[i];
    List data          = subjecti["data"];  // extract data-model options
    arma::vec RT       = data["RT"];
    arma::ucube model  = data.attr("model");
    unsigned int npda  = data.attr("n.pda");
    double bw          = data.attr("bw");
    unsigned int gpuid = data.attr("gpuid");
    NumericVector modelAttr = data.attr("model");
    std::string type = modelAttr.attr("type");      // model type
    arma::vec allpar = modelAttr.attr("all.par");
    List modelDim    = modelAttr.attr("dimnames");
    arma::umat n1idx = modelAttr.attr("n1.order");
    arma::uvec mc    = modelAttr.attr("match.cell");
    bool posd        = modelAttr.attr("posdrift");
    arma::uvec ise   = data.attr("cell.empty");
    cellidxes(i)     = cellIdx2Mat(data);
    std::vector<std::string> parname = modelAttr.attr("par.names");
    std::vector<std::string> dim1 = modelDim[0]; // row;
    std::vector<std::string> dim2 = modelDim[1]; // col; parameters
    std::vector<std::string> dim3 = modelDim[2]; // slice; response types; r1, r2
    arma::uvec isr1 = GetIsR1(modelAttr, type);

    arma::cube thetai = subjecti["theta"]; // nchain x npar x nmc
    arma::mat lpi     = subjecti["summed_log_prior"]; // nmc x nchain
    arma::mat lli     = subjecti["log_likelihoods"];  // nmc x nchain
    unsigned int start   = subjecti["start"];        // R index
    unsigned int start_C = start - 1;
    store_i(i) = start - 1;
    // Rcout << "store_i: " << store_i(i) << std::endl;

    types[i] = type;
    allpars(i) = allpar;
    n1idxes(i) = n1idx;
    matchcells(i) = mc;
    emptycells(i) = ise;
    parnames(i) = parname;
    dim1s(i) = dim1;
    dim2s(i) = dim2;
    dim3s(i) = dim3;
    isr1s(i) = isr1;
    models(i) = model;
    npdas(i) = npda;
    bws(i) = bw;
    gpuids(i) = gpuid;
    RTs(i) = RT;
    posdrift(i) = posd;

    thetas(i)    = thetai; // nchain x npar x nmc
    usethetas(i) = thetai.slice(start_C); // nchain x npar
    logpriors(i) = lpi; // nmc x nchain
    loglikes(i)  = lli;  // nmc x nchain
    uselogpriors(i) = arma::trans(lpi.row(start_C)); // nchain x 1
    useloglikes(i)  = arma::trans(lli.row(start_C));
  }
}

void GetPrior(List pprior, std::vector<std::string>& dists, arma::vec& p1,
  arma::vec& p2, arma::vec& lower, arma::vec& upper, arma::uvec& islog) {

  std::vector<std::string> pnames = pprior.attr("names");
  List L1;
  unsigned int npar = pprior.size();
  for (size_t i = 0; i < npar; i++) {
    L1 = pprior[pnames[i]];
    std::string disti = L1.attr("dist");
    dists[i] = disti;
    p1(i) = L1[0];
    p2(i) = L1[1];
    lower(i) = L1[2];
    upper(i) = L1[3];
    islog(i) = L1[4];
  }
}


arma::uvec GetIsR1(NumericVector modelAttr, std::string type) {
  arma::uvec out;
  if (type == "rd") {
    arma::uvec tmp = modelAttr.attr("is.r1");
    out = tmp;
  } else {
    out.zeros(1);
  }
  return out;
}


// arma::uvec Get_r1(std::string type, NumericVector model_attr) {
//   arma::uvec out;
//   if (type == "rd") {
//     arma::uvec tmp = model_attr.attr("is.r1");
//     out = tmp;
//   } else {
//     List modelDim = model_attr.attr("dimnames");
//     std::vector<std::string> dim1 = modelDim[0] ; // row; s1.r1, etc
//     unsigned int ncondition = dim1.size();
//     out.zeros(ncondition);
//   }
//   return out;
// }


arma::uvec getArmaIdx(NumericVector x, unsigned int trueOrFalse) {
  arma::vec xArma(x.begin(), x.size(), false) ;
  arma::uvec index = find(xArma == trueOrFalse) ;
  return index ;
}


// [[Rcpp::export]]
void CheckPnames(List samples) {
  List data          = samples["data"];
  NumericVector modelAttr = data.attr("model");
  NumericVector pvec = modelAttr.attr("p.vector");
  List pprior = samples["p.prior"];
  CharacterVector pnames1 = pvec.attr("names");
  CharacterVector pnames2 = pprior.attr("names");

  unsigned int npar1 = pvec.size();
  unsigned int npar2 = pprior.size();
  if (npar1 != npar2) stop("Does p.vector has more elements than p.prior?");

  for (size_t i = 0; i < npar1; i++) {
    if (pnames1[i] != pnames2[i]) {
      Rcout << pnames1[i] << " in p.vector does not match " << pnames2[i] <<
        " in p.prior. Wrong sequence? " << std::endl;
      stop("Check p.vector and p.prior.");
    }
  }
}


// [[Rcpp::export]]
void CheckHyperPnames(List samples) {

  List hyper = samples.attr("hyper");
  CharacterVector pnames0 = hyper["p.names"];

  List subject0 = samples[0];
  List pprior   = subject0["p.prior"];
  List ppprior  = hyper["pp.prior"];
  List lprior   = ppprior[0];
  List sprior   = ppprior[1];

  List data     = subject0["data"];
  NumericVector modelAttr = data.attr("model");
  NumericVector pvec = modelAttr.attr("p.vector");
  CharacterVector pnames1 = pvec.attr("names");
  CharacterVector pnames2 = pprior.attr("names");
  CharacterVector pnames3 = lprior.attr("names");
  CharacterVector pnames4 = sprior.attr("names");

  unsigned int npar0 = pnames0.size();
  unsigned int npar1 = pvec.size();
  unsigned int npar2 = pprior.size();
  unsigned int npar3 = lprior.size();
  unsigned int npar4 = sprior.size();

  if (npar0 != npar1) stop("Does hyper attributes set correctly?");
  if (npar1 != npar2) stop("Does p.vector has more elements than p.prior?");
  if (npar1 != npar3) stop("Does p.vector has more elements than mu.prior?");
  if (npar1 != npar4) stop("Does p.vector has more elements than sigma.prior?");

  for (size_t i = 0; i < npar1; i++) {
    if (pnames1[i] != pnames0[i]) {
      Rcout << pnames1[i] << " in p.vector does not match " << pnames0[i] <<
        " in hyper p.prior. Wrong sequence? " << std::endl;
      stop("Check p.vector and p.prior in hyyper.");
    }

    if (pnames1[i] != pnames2[i]) {
      Rcout << pnames1[i] << " in p.vector does not match " << pnames2[i] <<
        " in p.prior. Wrong sequence? " << std::endl;
      stop("Check p.vector and p.prior.");
    }

    if (pnames1[i] != pnames3[i]) {
      Rcout << pnames1[i] << " in p.vector does not match " << pnames3[i] <<
        " in mu.prior. Wrong sequence? " << std::endl;
      stop("Check p.vector and mu.prior.");
    }

    if (pnames1[i] != pnames4[i]) {
      Rcout << pnames1[i] << " in p.vector does not match " << pnames4[i] <<
        " in sigma.prior. Wrong sequence? " << std::endl;
      stop("Check p.vector and sigma.prior.");
    }
  }
}


//' Convert cell index list to a matrix
//'
//' Convert cell index list to a matrix. This is a convenient function.
//'
//' @param data a model data instance
//' @return a matrix
//'
//' @examples
//' m1 <- BuildModel(
//'     p.map     = list(a="1",v="1",z="1",d="1",sz="1",sv="1",t0="1",st0="1"),
//'     match.map = list(M=list(s1="r1", s2="r2")),
//'     factors   = list(S=c("s1", "s2")),
//'     constants = c(st0=0, d=0),
//'     responses = c("r1","r2"),
//'     type      = "rd")
//' p.prior <- BuildPrior(
//'     dists = rep("tnorm", 6),
//'     p1    = c(a=2,  v=2.5, z=.5, sz=.3, sv=1,  t0=.3),
//'     p2    = c(a=.5, v=.5,  z=.1, sz=.1, sv=.3, t0=.05),
//'     lower = c(0,-5, 0, 0, 0, 0),
//'     upper = c(5, 7, 2, 2, 2, 2))
//' dat <- simulate(m1, nsim = 4, prior = p.prior)
//' dmi <- BuildDMI(dat, m1)
//' cellIdx2Mat(dmi)
//' @export
// [[Rcpp::export]]
arma::umat cellIdx2Mat(List data) {

  List cellidxlist = data.attr("cell.index"); // s1.r1,
  arma::uvec cell1 = cellidxlist[0];
  unsigned int ncondition  = cellidxlist.size();  // cell.index names are dim1
  unsigned int ntrial = cell1.n_elem;
  arma::umat out(ntrial, ncondition);

  for (size_t i = 0; i < ncondition; i++) {
    arma::uvec tmp = cellidxlist[i];
    out.col(i) =  tmp;
  }
  return(out);
}


NumericMatrix na_matrix(unsigned int nr, unsigned int nc) {
  NumericMatrix out(nr, nc);
  std::fill( out.begin(), out.end(), NumericVector::get_na() ) ;
  return out;
}


// [[Rcpp::export]]
arma::mat ac_(arma::vec x, unsigned int nlag) {
  unsigned int n = x.n_elem;

  arma::mat out(n, nlag); out.fill(NA_REAL);
  out.col(0) = x;

  for (size_t i = 1; i < nlag; i++) {
    arma::vec tmp1 = arma::shift(x, (int)i);
    tmp1.rows(0, i-1).fill(NA_REAL);
    out.col(i) = tmp1;
  }
  return(out);
}

