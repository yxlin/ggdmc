List run_hyper_glm(List samples, unsigned int report, double pm, double pm0,
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

  arma::field<arma::mat> usephi(2);
  usephi(0) = location.slice(start_C); // nchain x nhpar
  usephi(1) = scale.slice(start_C);
  arma::vec usehlp = arma::trans(hlp.row(start_C)); // nchain
  arma::vec usehll = arma::trans(hll.row(start_C));

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
    /* ----------------------------------------------------------
    *  Update data level and hyper data
    * ----------------------------------------------------------*/
    for (size_t j = 0; j < nsub; j++) {
      // usephi: nchain x npar; usethetas(j): nchain x npar
      uselp(j) = UpdatePriors_glm(usetheta(j), dist_pp, p1, p2,
            dist_lp, dist_sp, usephi(0), usephi(1), lower, upper, lg);

      rj = arma::zeros<arma::uvec>(nchain);

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