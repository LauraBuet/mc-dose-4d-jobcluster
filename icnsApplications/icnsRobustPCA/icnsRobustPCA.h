/** \file imiRPCA.h
 *
 *  \b Initial \b Author: Matthias Wilms, Jonas Ortmueller \n\n
 *  \b Copyright (C) 2015 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/
#ifndef __imiRPCA_h
#define __imiRPCA_h

#include <armadillo>
// ITK includes

// VNL includes:

// Project includes:
//#include "imiObject.h"
//#include "imiITKImageTypeDefinitions.h"

namespace imi
{
/** \class imi::imiRPCA
 *  \brief This class implements the RPCA decomposition of matrix D into a low-rank
 *  matrix L and a sparse matrix S using the inexact augmented lagrangian multiplier method:
 *
 *  \min_{L,S} ||L||_* + \lambda ||S||_1, s.t. ||D-L-S||_{fro} < eps
 *
 * See
 * [1] Candes et al., "Robust Principal Component Analysis?",
 *   In: Journal of the ACM, Vol. 58, No. 3, 2011
 *
 * [2] Lin et al., "The Augmented Lagrangian Multiplier Method
 *  for Exact Recovery of Corrupted Low-Rank Matrices", 2011
 *
 *  for futher information.
 *
 *  <small> <!--Copyright Information: -->
 *  \b Initial \b Author: Matthias Wilms \n
 *  \b Copyright (C) 2015 Institute of Medical Informatics,
 *     University of Luebeck\n
 *  </small>
 ******************************************************************************/
class imiRPCA//: public imiObject
{
public:
  /* Class macro needed because class inherits from imiObject. */
  //imiClassMacro(imiRPCA,imiObject)
  //;

  /** Image type of the temporal images. */
  typedef double ValueType;
  typedef arma::vec ArmaVectorType;
  typedef arma::uvec ArmaUVectorType;
  typedef arma::mat ArmaMatrixType;

  /** \brief Create an instance of the object, use Delete() to destroy. */
  static imiRPCA* New();

  virtual double GetLambda(){ return m_Lambda; }
  virtual void SetLambda( double lambda ){ m_Lambda = lambda; }

  virtual double GetTolerance(){ return m_Tolerance; }
  virtual void SetTolerance( double tolerance ){ m_Tolerance = tolerance; }

  virtual int GetMaxIterations(){ return m_MaxIterations; }
  virtual void SetMaxIterations( int maxIterations ){ m_MaxIterations = maxIterations; };

  virtual ArmaMatrixType GetD(){ return m_D; }
  void SetD(ArmaMatrixType D){m_D=D;};

  virtual ArmaMatrixType GetL(){ return m_L; }
  virtual ArmaMatrixType GetS(){ return m_S; }

  bool Execute();

protected:
  imiRPCA();
  virtual ~imiRPCA();

  double m_Lambda;
  double m_Tolerance;
  int m_MaxIterations;
  ArmaMatrixType m_D;
  ArmaMatrixType m_L;
  ArmaMatrixType m_S;
};

} /* namespace imi */
#endif /* __imiRPCA_h */

