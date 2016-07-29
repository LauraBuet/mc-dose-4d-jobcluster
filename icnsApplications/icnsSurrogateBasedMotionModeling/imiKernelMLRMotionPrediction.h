/** \file imiKernelMLRMotionPrediction.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiKernelMLRMotionPrediction_h
#define __imiKernelMLRMotionPrediction_h


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

// Project includes:
#include "imiObject.h"
#include "imiITKImageTypeDefinitions.h"
#include "imiKernel.h"
#include "imiGaussianKernel.h"
#include "imiLinearKernel.h"

//#include "imiImageThreads.h"


namespace imi
{

class imiKernelMLRMotionPrediction : public imiObject
{

public:

  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiKernelMLRMotionPrediction,imiObject);

  /** Create an instance of the object, use Delete() to destroy. */
  static imiKernelMLRMotionPrediction* New();

  // ----------------------------
  //   Typedefs:
  // ----------------------------

  typedef imiRealType MatrixValueType;
  typedef vnl_vector<MatrixValueType> VnlVectorType;
  typedef vnl_matrix<MatrixValueType> VnlMatrixType;

  // ----------------------------
  //   Public Methods:
  // ----------------------------

  // --- Setter / Getter --------

  void SetMeanRegressand( VnlMatrixType& meanRegressand )
  {
    m_meanObservationVector = meanRegressand;
  }

  void SetRegressandTrainingMatrix( VnlMatrixType& regressandObservations )
  {
    m_observationMatrix = regressandObservations;
  }
  void SetRegressorTrainingMatrix( VnlMatrixType& regressorObservations )
  {
    m_regressorMatrix = regressorObservations;
  }

  /** \brief Get predicted output based on trained LS estimator.
   *
   * Requires the LS estimator to be trained and the regressor measurement
   * to be set.
   */
  bool GetPredictedOutput( VnlMatrixType& predictedOutput )
  {
    bool ret = true;
    if( m_predictedOutput.empty() )
    {
      ret = this->PredictOutput();
      if( ret ) predictedOutput = m_predictedOutput;
    }
    return ret;
  }

  /** \brief Get mean regressor matrix. */
  bool GetMeanRegressand( VnlMatrixType& meanRegressand )
  {
    if( m_meanObservationVector.empty() ) return false;
    meanRegressand = m_meanObservationVector;
    return true;
  }


  void SetMulticollinearityCheck( bool flag )
  {
    m_multicollinearityCheckFlag = flag;
  }
  void MulticollinearityCheckOn( void )
  {
    m_multicollinearityCheckFlag = true;
  }
  void MulticollinearityCheckOff( void )
  {
    m_multicollinearityCheckFlag = false;
  }

  // ----------------------------

  bool TrainLSEstimator( void );

  bool MulticollinearityCheck( VnlMatrixType& regressorCovarianceMatrix );

  /** \brief Predicts an output based on a regressor measurement.
   *
   * Predicts on output based on a regressor measurement.
   * It is assumed that the LS estimator is trained.
   *
   * @param regressorMeasurement
   * @param predictedOutput
   * @return success
   */
  bool PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput );
  bool PredictOutput( void )
  {
    return PredictOutput( m_regressorMeasurement, m_predictedOutput );
  }

  void SetKernel( imiKernel* kernel )
  {
    m_kernel = kernel;
  }

  void SetTikhonovRegularizationParameter( MatrixValueType tikhonovRegularizationParameter )
  {
    m_TikhonovRegularizationParameter = tikhonovRegularizationParameter;
    m_invCenteredK.clear();
  }

protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  imiKernelMLRMotionPrediction();

  /** \brief Destructor */
  virtual ~imiKernelMLRMotionPrediction();


  // ----------------------------
  //   Protected attributes.
  // ----------------------------


  bool m_multicollinearityCheckFlag;
  MatrixValueType m_TikhonovRegularizationParameter;

  VnlMatrixType m_regressorMatrix;
  VnlMatrixType m_observationMatrix;

  VnlMatrixType m_meanObservationVector;
  VnlMatrixType m_centeredObservationMatrix;

  VnlMatrixType m_meanRegressorVector;
  VnlMatrixType m_centeredRegressorMatrix;

  VnlMatrixType m_K; //Kernel matrix
  VnlMatrixType m_centeredK; //Centered kernel matrix
  VnlMatrixType m_H; //Centering matrix
  VnlMatrixType m_invCenteredK; //inverse centered kernel matrix


  imiKernel* m_kernel;


  VnlMatrixType m_regressorMeasurement;
  VnlMatrixType m_predictedOutput;

}; // end class

}

#endif /* __imiKernelMLRMotionPrediction_h */
