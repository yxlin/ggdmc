#include <ggdmc.hpp>

using namespace Rcpp;

lba::lba (double A, double b, double mean_v, double sd_v, double t0,
          bool posdrift) :
  m_A(A), m_b(b), m_mean_v(mean_v), m_sd_v(sd_v), m_t0(t0), is_posv(posdrift)
// pdf, cdf
{
  // double m_A, m_b, m_mean_v, m_sd_v, m_t0, m_st0; // LBA distribution.
  // bool is_posv;
  //
  // double *m_meanv_vec, *m_sdv_vec;
  // unsigned int m_nmean_v;
}

lba::lba (double A, double b, double * mean_v, double * sd_v, double t0,
          double st0, unsigned int & nmean_v, bool posdrift) :
  m_A(A), m_b(b), m_t0(t0), m_st0(st0), is_posv(posdrift),
  m_meanv_vec(mean_v), m_sdv_vec(sd_v), m_nmean_v(nmean_v)
// rlba_norm. NOTE double * mean_v and double * sd_v
{
    if (m_st0 < 0) Rcpp::stop("st0 must be greater than 0.");
}

lba::~lba() {}


void lba::d (std::vector<double> & x, double * output)
{
    for(size_t i=0; i<x.size(); i++)
    {
      output[i] = d(x[i]);
    }
  }
void lba::d (std::vector<double> & x, std::vector<double> & output)
{
    for(size_t i=0; i<x.size(); i++)
    {
      output[i] = d(x[i]);
    }
  }
void lba::d (arma::vec & x, arma::vec & output)
{
    for(size_t i=0; i<x.size(); i++)
    {
      output[i] = d(x[i]);
    }
  }

void lba::p (std::vector<double> & x, double * output)
{
    for(size_t i=0; i<x.size(); i++)
    {
      output[i] = p(x[i]);
    }
  }
void lba::p (std::vector<double> & x, std::vector<double> & output)
{
    for(size_t i=0; i<x.size(); i++)
    {
      output[i] = p(x[i]);
    }
  }
void lba::p (arma::vec & x, arma::vec & output)
{
    for(size_t i=0; i<x.size(); i++)
    {
      output[i] = p(x[i]);
    }
  }

void lba::node1_pdf (arma::vec & x, arma::vec & output)
{
    if (x.size() != output.size()) Rcpp::stop("unequal sizes");

    for(size_t i=0; i<x.size(); i++)
    {
      output[i] *= 1. - p(x[i]);
    }
  }

void lba::r (unsigned int & n, arma::mat & output)
{
    // output n x 2
    arma::vec tmp(m_nmean_v);

    for (size_t i=0; i<n; i++)
    {
      rt(tmp);
      output(i, 0) = tmp.min();
      output(i, 1) = 1 + tmp.index_min(); // plus 1 to fit R indexing
    }

  }

void lba::print(const std::string & x) const
{
    Rcpp::Rcout << x << "[A, b, mean_v, sd_v, t0]: " << m_A << ", " << m_b <<
      ", " << m_mean_v << ", " << m_sd_v << ", " << m_t0 << std::endl;
  }

arma::vec n1PDFfixedt0(arma::vec rt, arma::vec A, arma::vec b,
                       arma::vec mean_v, arma::vec sd_v, arma::vec t0,
                       bool posdrift) {
  // No check for a rectangle matrix of parameters
  unsigned int nmean_v = mean_v.n_elem; // Number of accumulators/responses.
  unsigned int n       = rt.n_elem;     // Number of RTs

  arma::vec out(n);

  lba * obj = new lba(A[0], b[0], mean_v[0], sd_v[0], t0[0], posdrift);
  obj->d(rt, out);

  if (nmean_v > 1)
  {
    for (size_t i = 1; i < nmean_v; i++)
    {
      obj->m_A      = A[i];
      obj->m_b      = b[i];
      obj->m_mean_v = mean_v[i];
      obj->m_sd_v   = sd_v[i];
      obj->m_t0     = t0[i];

      obj->node1_pdf(rt, out);
    }
  }

  delete obj;
  return out;
}

//' Generate Random Deviates of the LBA Distribution
//'
//' \code{rlba_norm}, only slightly faster than \code{maker}, calls C++
//' function directly.
//'
//' @param n is the numbers of observation.
//' @param A start point upper bound, a vector of a scalar.
//' @param b decision threshold, a vector or a scalar.
//' @param mean_v mean drift rate vector
//' @param sd_v standard deviation of drift rate vector
//' @param t0 nondecision time, a vector.
//' @param st0 nondecision time variation, a vector.
//' @param posdrift if exclude negative drift rates
//'
//' @return a n x 2 matrix of RTs (first column) and responses (second column).
//' @export
// [[Rcpp::export]]
arma::mat rlba_norm(unsigned int n, arma::vec A, arma::vec b,
                    arma::vec mean_v, arma::vec sd_v, arma::vec t0,
                    arma::vec st0, bool posdrift) {
  unsigned int nmean_v = mean_v.size();
  unsigned int nA  = A.n_elem;
  unsigned int nb  = b.n_elem;
  unsigned int nt0 = t0.n_elem;
  unsigned int nst0= st0.n_elem;

  if (nA == 1) A = arma::repmat(A, nmean_v, 1);
  if (nb == 1) b = arma::repmat(b, nmean_v, 1);
  if (nt0 == 1) t0 = arma::repmat(t0, nmean_v, 1);
  if (sd_v.n_elem == 1) sd_v = arma::repmat(sd_v, nmean_v, 1);
  if (nst0 == 1) st0 = arma::repmat(st0, nmean_v, 1);

  double * mv  = new double[nmean_v];
  double * sdv = new double[nmean_v];
  for(size_t i=0; i<nmean_v; i++)
  {
    mv[i]  = mean_v[i];
    sdv[i] = sd_v[i];
  }

  lba * obj = new lba(A[0], b[0], mv, sdv, t0[0], st0[0], nmean_v, posdrift);

  arma::mat out(n, 2);
  obj->r(n, out);

  delete obj;
  delete [] sdv;
  delete [] mv;

  return out;
}


