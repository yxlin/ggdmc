#include <RcppArmadillo.h>

void MigrateDMCHyperchains(arma::field<arma::mat>& usephi,
                           arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
                           arma::uvec pdists, arma::vec plower, arma::vec pupper,
                           arma::uvec plog, arma::uvec ldists, arma::vec lp1,
                           arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
                           arma::uvec sdists, arma::vec sp1, arma::vec sp2,
                           arma::vec slower, arma::vec supper, arma::uvec slog, double rp,
                           arma::uvec& rj, bool debug);
void MigrateDMCHyperchains_old(arma::field<arma::mat>& usephi,
                               arma::vec& usehlp, arma::vec& usehll,  arma::cube theta,
                               arma::uvec pdists, arma::vec plower, arma::vec pupper,
                               arma::uvec plog, arma::uvec ldists, arma::vec lp1,
                               arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
                               arma::uvec sdists, arma::vec sp1, arma::vec sp2,
                               arma::vec slower, arma::vec supper, arma::uvec slog, double rp,
                               arma::uvec& rj, bool debug);
void CrossoverDMCHyperchains_blocked(arma::field<arma::mat>& usephi,
                                     arma::vec& usehlp, arma::vec& usehll, arma::cube theta,
                                     arma::uvec pdists, arma::vec plower, arma::vec pupper,
                                     arma::uvec plog, arma::uvec ldists, arma::vec lp1,
                                     arma::vec lp2, arma::vec llower, arma::vec lupper, arma::uvec llog,
                                     arma::uvec sdists, arma::vec sp1, arma::vec sp2,
                                     arma::vec slower, arma::vec supper, arma::uvec slog, double rp,
                                     double gammaMult, arma::uvec& rj, unsigned int j);

arma::vec GetGamma(unsigned int npar, double gammamult, bool hyper = false);
arma::uvec PickChains(unsigned int k, unsigned int nchain, arma::uvec chains);
arma::uvec GetSubchains(unsigned int nchain, bool debug);
arma::uvec SelectEmigrants(unsigned int ngroup, unsigned int k);