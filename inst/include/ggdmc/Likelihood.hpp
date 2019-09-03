#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

#include <RcppArmadillo.h>

using namespace Rcpp;

class Likelihood
{
private:
  enum ModelType { DEFAULT, DDM, LBA };

  ModelType resolve_option(std::string type)
  {
    if(type == "rd")   return DDM; // 1
    if(type == "norm") return LBA; // 2
    return DEFAULT;  // 0
  }

public:
  Design * m_d;
  std::string m_mtype;

  // All additional arguments
  arma::uvec m_is_r1; // DDM
  double m_precision; // DDM

  arma::umat m_n1idx;  // LBA
  bool m_posdrift, m_n1order;

  Likelihood(List & dmi, Design * d) : m_d(d)
  // Run constructor
  {
    NumericVector modelAttr = dmi.attr("model");

    arma::umat tmp_n1idx = modelAttr.attr("n1.order");
    std::string type     = modelAttr.attr("type");

    m_n1idx = tmp_n1idx;

    // TODO: make it compatible to old version
    arma::uvec tmp_isr1 = modelAttr.attr("is.r1");
    m_mtype = type;

    if (m_mtype == "rd") m_is_r1 = tmp_isr1;

    m_precision = 3.0;
    m_posdrift  = modelAttr.attr("posdrift");
    m_n1order   = true;
  }

  Likelihood(std::string model_type,
             arma::uvec isr1,
             arma::umat n1idx,
             bool n1order,
             Design * d) :
    m_d(d), m_mtype(model_type), m_is_r1(isr1), m_n1idx(n1idx)
  // p_df constructor, random functions set n1order as false
  {
    m_precision  = 3.0;
    m_posdrift   = true;
    m_n1order    = n1order;
  }

  ~Likelihood()
  {
    delete m_d;
    // Rcout << "Likelihood destructor\n";
  }

  /* ---------------Model-specific methods--------------- */

  void transform (arma::mat & output, std::string & cell)
  // DDM transform
  {
    // parmat 8 x 2

    for (size_t i = 0; i < m_d->m_nc; i++)
    {
      // dim0: // "s1.r1" "s2.r1" "s1.r2" "s2.r2"
      // ir1 is an index vector, indicating:
      //
      // If this is the selected condition (ie cell), plus it is a
      // correct response, plus it is the first response, then flip zr.
      // TODO: NEED FURTHER INVESTIGATION. This may be a AH's bug.
      if (m_d->m_dim0[i] == cell && m_is_r1[i])
      {
        output.row(2) = 1 - output.row(2);
      }
    }

  }

  arma::mat transform (arma::mat & parmat, std::string & cell,
                       bool n1order)
  // LBA transform
  {
    // A b  t0 mean_v sd_v st0
    parmat.row(1) += parmat.row(0); // calculate b = A + B
    arma::mat out = parmat;

    // n1idx: nc x nr

    if (n1order)
    {
      for (size_t i=0; i<m_d->m_nc; i++)
      {
        if (m_d->m_dim0[i] == cell)
        {
          for (size_t j=0; j < m_d->m_nr; j++)
          {
            out.col(j) = parmat.col( m_n1idx(i, j) - 1 );
          }
        }
      }
    }

    return out;
  }


  arma::vec ddm (arma::vec & pvector)
  {
    double * para = new double[m_d->m_nParameter];

    arma::mat pmat(m_d->m_nParameter, m_d->m_nr); // 8 x 2
    arma::vec out(m_d->m_nRT);

    // Parameters * params;
    arma::uvec RTIdx, tmp;
    arma::vec selectedRT;

    // [a   v   zr   d  szr  sv  t0 st0]
    for (size_t i=0; i<m_d->m_nc; i++)
    {
      if (!m_d->m_is_empty_cell[i])
      {
        parameter_matrix(pvector, m_d->m_dim0[i], pmat);
        transform(pmat, m_d->m_dim0[i]);

        for (size_t j = 0; j < m_d->m_nParameter; j++) para[j] = pmat(j, 0);

        Parameters * params = new Parameters(para, m_precision);

        tmp = m_d->m_is_this_cell.col(i);
        RTIdx = arma::find(m_d->m_is_this_cell.col(i));
        selectedRT = m_d->m_RT(RTIdx);

        if (!params->ValidateParams(false)) // do not print invalid parameters
        {
          out(RTIdx).fill(1e-10);
        }
        else
        {
          if (m_d->m_is_matched_cell[i])  // choose g_plus
          {
            for (size_t k = 0; k < selectedRT.n_elem; k++)
            {
              out(RTIdx(k)) = R::fmax2(std::abs(g_plus(selectedRT(k), params)), 1e-10);
            }
          }
          else
          {
            for (size_t k = 0; k < selectedRT.n_elem; k++)
            {
              out(RTIdx(k)) = R::fmax2(std::abs(g_minus(selectedRT(k), params)), 1e-10);
            }
          }
        }

        delete params;


      }
    }

    delete [] para;
    return out;
  }

  arma::vec lba_ (arma::vec & pvector)
  {
    arma::mat pmat0(m_d->m_nParameter, m_d->m_nr);
    arma::vec out(m_d->m_nRT);
    arma::uvec RTIdx;

    // pmat: A   b  t0 mean_v sd_v st0
    for (size_t i=0; i<m_d->m_nc; i++)
    {
      if (!m_d->m_is_empty_cell[i])
      {
        parameter_matrix(pvector, m_d->m_dim0[i], pmat0);

        arma::mat pmat = transform(pmat0, m_d->m_dim0[i], m_n1order);
        pmat = arma::trans(pmat);

        RTIdx = arma::find(m_d->m_is_this_cell.col(i) == 1);
          out(RTIdx) = n1PDFfixedt0(m_d->m_RT(RTIdx), pmat.col(0), pmat.col(1),
              pmat.col(3), pmat.col(4), pmat.col(2), pmat.col(5), m_posdrift);
      }

    }

    // for (size_t i = 0; i < m_d->m_nRT; i++) out[i] = R::fmax2(out[i], 1e-10);
    return out;
  }

  /* ---------------Generic methods--------------- */
  void parameter_matrix(arma::vec & pvector, std::string & cell,
                        arma::mat & output)
  {
    arma::vec tmp(m_d->m_nParameter);

    // Iterate through accumulators, eg r1, r2, (if have any), r3, etc.
    // parnames is the parameter the user wishes to fit, which may/may not
    // differ from the parameter likelihood function expects
    // eg A, B, mean_v, sd_v, t0, st0 vs A, b, mean_v, sd_v, t0, st0
    // likelihood function expects to see b, but the user is asked to enter B.

    for (size_t i=0; i<m_d->m_nr; i ++) // r1, r2, r3, etc.
    {
      for (size_t j=0; j<m_d->m_nc; j++) // eg s1.r1, s2.r1, s1.r2, s2.r2 etc.
      {
        if (m_d->m_dim0[j] == cell)
        {
          size_t idx = 0;

          for (size_t k=0; k<m_d->m_np; k++) // eg allpars a, v.f1, v.f2, z, d, sz) ...
          {

            if (m_d->m_model(j, k, i))
            {
              // The values in allpar vector is either constant or NA.
              // When a value is NA, it indicates that it is not a constant.
              // In this case, its value is stored in pvector
              // allpar names == dim2 names

              tmp[idx] = m_d->m_allpar[k]; // constant values or NA

              for(size_t l = 0; l < m_d->m_npar; l++)
              {
                // replace NA with values in p.vector.
                if (m_d->m_pnames[l] == m_d->m_dim1[k] && ISNAN(tmp[idx]))
                {
                  tmp[idx] = pvector[l];
                }
              }

              idx++;

            }
          }
        }
      }

      output.col(i) = tmp;
    }
  }

  double sumloglike(arma::vec pvector) // mustn't pass memory location
  {
    arma::vec den;

    switch(resolve_option(m_mtype))
    {
    case DDM:
      den = ddm(pvector);
      break;
    case LBA:
      den = lba_(pvector);
      break;
    case DEFAULT:
      Rcout << "Undefined model type\n"; den.fill(1e-10);
      break;
    default:
      Rcout << "Unexpected situation\n"; den.fill(1e-10);
      break;
    }

    double out = arma::accu(arma::log(den));
    if (ISNAN(out)) out = R_NegInf;

    return out;
  }

  arma::vec likelihood(arma::vec & pvector)
  {

    arma::vec den;

    switch(resolve_option(m_mtype))
    {
    case DDM:
      den = ddm(pvector);
      break;
    case LBA:
      den = lba_(pvector);
      break;
    case DEFAULT:
      Rcout << "Undefined model type\n"; den.fill(1e-10);
      break;
    default:
      Rcout << "Unexpected situation\n"; den.fill(1e-10);
      break;
    }

    return den;
  }

  arma::mat get_pmat(arma::vec & pvector, std::string & cell)
  {
    arma::mat pmat(m_d->m_nParameter, m_d->m_nr);

    switch(resolve_option(m_mtype))
    {
    case DDM:   // 1
      parameter_matrix(pvector, cell, pmat);
      transform(pmat, cell);
      break;
    case LBA:  // 2
      parameter_matrix(pvector, cell, pmat);
      pmat = transform(pmat, cell, m_n1order);
      break;
    case DEFAULT:  // 0
      Rcout << "Undefined model.\n"; pmat.fill(NA_REAL);
      break;
    default:
      Rcout << "Unexpected situation.\n"; pmat.fill(NA_REAL);
      break;
    }

    return pmat.t();
  }


};

#endif
