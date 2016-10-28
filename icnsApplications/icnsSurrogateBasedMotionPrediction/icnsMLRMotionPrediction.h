/** \file icnsMLRMotionPredictionMethods.h
 *
 *  Original implementation:
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 *  Modified by Rene Werner, ICNS, 2016.
 *
 ****************************************************************************/

#ifndef __imiMLRMotionPredictionMethods_h
#define __imiMLRMotionPredictionMethods_h

// System includes:
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <vcl_iostream.h>

// vnl includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_rank.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"
#include "vnl/algo/vnl_cholesky.h"

namespace imi
{

class icnsMLRMotionPrediction
{

public:

  /** Create an instance of the object, use Delete() to destroy. */
  static icnsMLRMotionPrediction* New();

  // ----------------------------
  //   Typedefs:
  // ----------------------------

  typedef float MatrixValueType;
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

  void SetMeanRegressor( VnlMatrixType& meanRegressor )
  {
    m_meanRegressorVector = meanRegressor;
  }

  void SetRegressandTrainingMatrix( VnlMatrixType& regressandObservations )
  {
    // Compute mean observation vector:
    m_meanObservationVector.set_size( regressandObservations.rows(), 1 );
    m_meanObservationVector.fill( 0.0 );

    std::cout << "  Computing mean regressand vector ... " << std::endl;
    for( unsigned int i = 0; i < regressandObservations.rows(); i++ )
    {
      m_meanObservationVector( i, 0 ) = (regressandObservations.get_row( i )).mean();
    }

    // Compute mean regressand matrix entries:

    std::cout << "  Computing regressand matrix with centered entries ... " << std::endl;
    m_centeredObservationMatrix.set_size( regressandObservations.rows(), regressandObservations.cols() );
    for( unsigned int i = 0; i < regressandObservations.rows(); i++ )
    {
      m_centeredObservationMatrix.set_row( i, (regressandObservations.get_row( i ) - m_meanObservationVector( i, 0 )) );
    }
  }

  void SetRegressorTrainingMatrix( VnlMatrixType& regressorObservations )
  {
    std::cout<<"observations: "<<regressorObservations<<std::endl;

    // Compute mean observation vector:
    m_meanRegressorVector.set_size( regressorObservations.rows(), 1 );
    m_centeredRegressorMatrix.set_size( regressorObservations.rows(), regressorObservations.cols() );
    m_meanRegressorVector.fill( 0.0 );

    std::cout << "  Computing mean regressor vector ... " << std::endl;
    for( unsigned int i = 0; i < regressorObservations.rows(); i++ )
    {
      m_meanRegressorVector( i, 0 ) = (regressorObservations.get_row( i )).mean();
    }

    // Compute mean regressand matrix entries:

    //std::cout<<"regressorMatrix: "<<m_regressorMatrix<<std::endl;

    std::cout << "  Computing regressor matrix with centered entries ... " << std::endl;
    //m_centeredObservationMatrix = m_observationMatrix;
    for( unsigned int i = 0; i < regressorObservations.rows(); i++ )
    {
      m_centeredRegressorMatrix.set_row( i, (regressorObservations.get_row( i ) - m_meanRegressorVector( i, 0 )) );
    }
    /*std::cout<<"regressorObservations: "<<regressorObservations<<std::endl;
     std::cout<<" m_meanRegressorVector: "<< m_meanRegressorVector<<std::endl;
     std::cout<<" m_centeredRegressorMatrix: "<< m_centeredRegressorMatrix<<std::endl;*/
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
      if( ret )
        predictedOutput = m_predictedOutput;
    }
    return ret;
  }

  /** \brief Get mean regressor matrix. */
  bool GetMeanRegressand( VnlMatrixType& meanRegressand )
  {
    if( m_meanObservationVector.empty() )
      return false;
    meanRegressand = m_meanObservationVector;
    return true;
  }

  bool GetMeanRegressor( VnlMatrixType& meanRegressor )
  {
    if( m_meanRegressorVector.empty() )
      return false;
    meanRegressor = m_meanRegressorVector;
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

  bool GetTrainedEstimator( VnlMatrixType& estimator )
  {
    if( m_meanRegressorVector.empty() ) return false;
    estimator = m_trainedEstimator;
    return true;
  }

  void SetTrainedEstimator( VnlMatrixType& estimator )
  {
    m_trainedEstimator = estimator;
  }

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

  void SetTikhonovRegularizationParameter( MatrixValueType tikhonovRegularizationParameter )
  {
    m_TikhonovRegularizationParameter = tikhonovRegularizationParameter;
    m_invSquareMatrixPInv.clear();
  }

protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  icnsMLRMotionPrediction();

  /** \brief Destructor */
  virtual ~icnsMLRMotionPrediction();

  // ----------------------------
  //   Protected attributes.
  // ----------------------------

  bool m_multicollinearityCheckFlag;
  MatrixValueType m_TikhonovRegularizationParameter;

  VnlMatrixType m_meanObservationVector;
  VnlMatrixType m_centeredObservationMatrix;

  VnlMatrixType m_meanRegressorVector;
  VnlMatrixType m_centeredRegressorMatrix;

  VnlMatrixType m_squareMatrixPInv;
  VnlMatrixType m_invSquareMatrixPInv;

  VnlMatrixType m_regressorMeasurement;
  VnlMatrixType m_predictedOutput;

  VnlMatrixType m_trainedEstimator;

};
// end class

}

#endif /* __imiMLRMotionPredictionMethods_h */
