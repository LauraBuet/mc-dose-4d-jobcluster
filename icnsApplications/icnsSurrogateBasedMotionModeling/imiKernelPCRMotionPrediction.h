/** \file imiKernelPCRMotionPrediction.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiKernelPCRMotionPrediction_h
#define __imiKernelPCRMotionPrediction_h


// System includes:
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

// vnl includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_rank.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_cholesky.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_qr.h"


// Project includes:
#include "imiObject.h"
#include "imiITKImageTypeDefinitions.h"
#include "imiKernelPCA.h"
//#include "imiImageThreads.h"


namespace imi
{

class imiKernelPCRMotionPrediction : public imiObject
{

public:

  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiKernelPCRMotionPrediction,imiObject);

  /** Create an instance of the object, use Delete() to destroy. */
  static imiKernelPCRMotionPrediction* New();

  // ----------------------------
  //   Typedefs:
  // ----------------------------

  typedef imiRealType MatrixValueType;
  typedef vnl_vector<MatrixValueType> VnlVectorType;
  typedef vnl_matrix<MatrixValueType> VnlMatrixType;
  typedef vnl_diag_matrix<MatrixValueType> VnlDiagMatrixType;

  // ----------------------------
  //   Public Methods:
  // ----------------------------

  // --- Setter / Getter --------

//  void SetLSEstimator( VnlMatrixType& lsEstimator )
//  {
//    m_KernelPCREstimator = lsEstimator;
//  }
  void SetMeanRegressand( VnlMatrixType& meanRegressand )
  {
    m_meanRegressandVector = meanRegressand;
  }

  void SetRegressandTrainingMatrix( VnlMatrixType& regressandObservations )
  {
    m_regressandMatrix = regressandObservations;
  }
  void SetRegressorTrainingMatrix( VnlMatrixType& regressorObservations )
  {
    m_regressorMatrix = regressorObservations;
  }


//  /** \brief Get trained KernelPCR estimator. */
//  bool GetKernelPCREstimator( VnlMatrixType& lsEstimator )
//  {
//    if( m_KernelPCREstimator.empty() ) return false;
//    lsEstimator = m_KernelPCREstimator;
//    return true;
//  }

  /** \brief Get mean regressor matrix. */
  bool GetMeanRegressand( VnlMatrixType& meanRegressand )
  {
    if( m_meanRegressandVector.empty() ) return false;
    meanRegressand = m_meanRegressandVector;
    return true;
  }

  void SetNumberOfComponents( int numOfComponents )
  {
    m_numOfComponentsUsed=numOfComponents;
  }

  void SetVariabilityThreshold( MatrixValueType threshold )
  {
    m_variabilityThreshold=threshold;
  }

  void SetKernel(imiKernel* kernel)
  {
    m_kernel=kernel;
    m_kernelPCA->SetKernel(m_kernel);
  }

  // ----------------------------


  bool TrainKernelPCREstimator( void );

  /** \brief Predicts an output based on a regressor measurement.
   *
   * Predicts on output based on a regressor measurement.
   * It is assumed that the KernelPCR estimator is trained.
   *
   * @param regressorMeasurement
   * @param predictedOutput
   * @return success
   */
  bool PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput );

protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  imiKernelPCRMotionPrediction();

  /** \brief Destructor */
  virtual ~imiKernelPCRMotionPrediction();


  // ----------------------------
  //   Protected attributes.
  // ----------------------------



  VnlMatrixType m_regressorMatrix;
  VnlMatrixType m_regressandMatrix;

  VnlMatrixType m_meanRegressandVector;
  VnlMatrixType m_centeredRegressandMatrix;

  VnlMatrixType m_predictedOutput;

  VnlMatrixType m_principalComponents;
  VnlDiagMatrixType m_invEigenvalues;

  imiKernel* m_kernel;
  imiKernelPCA* m_kernelPCA;

  imiRealType m_variabilityThreshold;
  int m_numOfComponentsUsed;

}; // end class

}

#endif /* __imiKernelPCRMotionPrediction_h */
