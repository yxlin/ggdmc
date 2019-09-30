#ifndef DESIGN_HPP
#define DESIGN_HPP

#include <RcppArmadillo.h>

using namespace Rcpp;

class Design
// parameters = parnames; nParameter
// pvector of samples consists of part of the names of "parameters"; npar
// allpar parameters x conditions; np
{
public:
  unsigned int m_nc, m_np, m_nr, m_nRT, m_nParameter, m_npar, m_nallpar;
  double *m_allpar;
  std::string *m_pnames, *m_parameters, *m_dim0, *m_dim1, *m_dim2;
  bool *m_is_empty_cell, *m_is_matched_cell;

  arma::vec m_RT;
  arma::umat m_is_this_cell;
  arma::ucube m_model;

  Design(List & dmi)
  // Run & likelihood constructor
  {
    using namespace arma;
    using namespace std;

    // NOTE: all things from R must be casted to a valid R-recognizable type
    // first, so the messy syntax.
    // List dmi = samples["data"];  // data model instance

    NumericVector modelAttr  = dmi.attr("model");
    vector<bool> ise         = dmi.attr("cell.empty"); // nc
    List cidx                = dmi.attr("cell.index"); // nc elements
    ucube tmp_model          = dmi.attr("model");
    arma::vec RT             = dmi["RT"];

    List modelDim            = modelAttr.attr("dimnames");
    NumericVector pvector    = modelAttr.attr("p.vector"); // Carry NA values
    vector<string> parnames  = modelAttr.attr("par.names");
    vector<double> tmp_allpar= modelAttr.attr("all.par"); // with conditoin complications
    vector<bool> mc          = modelAttr.attr("match.cell");
    vector<string> pnames    = pvector.attr("names");
    // nc = dim0.size() ;  // number of condition, eg s1.r1 etc.
    // np = dim1.size() ;  // number of parameter x condition, eg v.f1, v.f2
    // nr = dim2.size() ;  // number of accumualtor/response, eg r1, r2
    vector<string> dim0 = modelDim[0];
    vector<string> dim1 = modelDim[1];
    vector<string> dim2 = modelDim[2];

    // 1. Get RTs (& response angles etc.)
    m_RT = RT;

    // 2. Get the sizes of all design-related variables
    m_nc          = dim0.size();
    m_np          = dim1.size();
    m_nr          = dim2.size();
    m_nRT         = RT.size();
    m_nParameter  = parnames.size(); // pure                model parameters (no condition)
    m_npar        = pvector.size();  // modeled             model-condition parameters
    m_nallpar     = tmp_allpar.size(); // modeled and fixed model-condition parameters

    // 3. Copy from stack to heap
    m_allpar          = new double[m_nallpar];
    m_pnames          = new string[m_npar];
    m_parameters      = new string[m_nParameter];
    m_dim0            = new string[m_nc];
    m_dim1            = new string[m_np];
    m_dim2            = new string[m_nr];
    m_is_matched_cell = new bool[m_nc];
    m_is_empty_cell   = new bool[m_nc];

    std::copy(tmp_allpar.begin(), tmp_allpar.end(), m_allpar);
    std::copy(pnames.begin(), pnames.end(), m_pnames);
    std::copy(parnames.begin(), parnames.end(), m_parameters);

    std::copy(dim0.begin(), dim0.end(), m_dim0);
    std::copy(dim1.begin(), dim1.end(), m_dim1);
    std::copy(dim2.begin(), dim2.end(), m_dim2);

    std::copy(ise.begin(), ise.end(), m_is_empty_cell);
    std::copy(mc.begin(), mc.end(), m_is_matched_cell);

    // 4. We rely on Armadillo methods to handle cellidx and model cube
    umat tmp_cidx(m_nRT, m_nc);
    for (size_t i = 0; i <m_nc; i++)
    {
      uvec tmp_vec = cidx[i];
      tmp_cidx.col(i) =  tmp_vec;
    }

    m_is_this_cell = tmp_cidx;
    m_model = tmp_model;

  }

  Design(std::vector<std::string> & pnames,
         std::vector<std::string> & parnames,
         std::vector<std::string> & dim0,
         std::vector<std::string> & dim1,
         std::vector<std::string> & dim2,
         std::vector<double> & allpar,
         arma::ucube & model) : m_model(model)
  // p_df constructor
  {
    m_nc  = dim0.size();  // number of condition, eg s1.r1 etc.
    m_np  = dim1.size();  // number of parameter x condition, eg v.f1, v.f2
    m_nr  = dim2.size();  // number of accumualtor/response, eg r1, r2
    m_nParameter = parnames.size();
    m_npar       = pnames.size();

    m_allpar     = new double[allpar.size()];
    m_pnames     = new std::string[m_npar];
    m_parameters = new std::string[m_nParameter];
    m_dim0       = new std::string[m_nc];
    m_dim1       = new std::string[m_np];
    m_dim2       = new std::string[m_nr];

    std::copy(pnames.begin(), pnames.end(), m_pnames);
    std::copy(parnames.begin(), parnames.end(), m_parameters);
    std::copy(dim0.begin(), dim0.end(), m_dim0);
    std::copy(dim1.begin(), dim1.end(), m_dim1);
    std::copy(dim2.begin(), dim2.end(), m_dim2);
    std::copy(allpar.begin(), allpar.end(), m_allpar);

  }

  ~Design()
  {
    // Rcout << "Design destructor\n";
  }
};

#endif
