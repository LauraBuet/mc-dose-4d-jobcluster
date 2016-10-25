/** \file icnsMLRMotionPrediction.cxx
 *
 *  Original authors:
 *
 *  \b Initial \b Authors: Matthias Wilms, Rene Werner \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 *  Modified by: Rene Werner (ICNS, 2016)
 *
 ****************************************************************************/

// Project includes
#include "icnsMLRMotionPrediction.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

icnsMLRMotionPrediction* icnsMLRMotionPrediction::New()
{
  return new icnsMLRMotionPrediction();
}

icnsMLRMotionPrediction::icnsMLRMotionPrediction()
{
  // Set default parameters:
  m_multicollinearityCheckFlag = false;
  m_TikhonovRegularizationParameter = 0.001;
  m_squareMatrixPInv.clear();
  m_invSquareMatrixPInv.clear();
}

// imiMLRMotionPredictionMethods::~imiMLRMotionPredictionMethods()
icnsMLRMotionPrediction::~icnsMLRMotionPrediction()
{
}

// ----------------------------
//   Public methods
// ----------------------------

// imiMLRMotionPredictionMethods::TrainLSEstimator
bool icnsMLRMotionPrediction::TrainLSEstimator( void )
{
  std::cout <<  "-----------------------------------" << std::endl;
  std::cout <<  "  TRAINING OF LS ESTIMATOR"          << std::endl;
  std::cout <<  "-----------------------------------" << std::endl;

  //-------------------------------------------
  // FIRST STEP:
  // Center regressor data
  //-------------------------------------------

  if( m_centeredRegressorMatrix.empty() )
  {
    std::cerr << "MLR: Regressor matrix empty! Aborting computation." << std::endl;
    return false;
  }

  //-------------------------------------------
  // First STEP:
  // Center observation data
  //-------------------------------------------

  if( m_centeredObservationMatrix.empty() )
  {
    std::cerr << "MLR: Observation matrix empty! Aborting computation." << std::endl;
    return false;
  }

  //-------------------------------------------
  // Second STEP:
  // Compute square matrix (covariance matrix).
  // Note: It is assumed that the regressor
  // data has already been centered. This is
  // automatically done when setting the
  // regressor data.
  //-------------------------------------------

  std::cout << "  Computing square matrix ZZ^T for pseudoinverse ... " << std::flush;
  if (m_centeredRegressorMatrix.rows() <= m_centeredRegressorMatrix.cols())
  {
    m_squareMatrixPInv = m_centeredRegressorMatrix * m_centeredRegressorMatrix.transpose();
  }
  else
  {
    m_squareMatrixPInv = m_centeredRegressorMatrix.transpose() * m_centeredRegressorMatrix;
  }
  std::cout << " Ok.\n ZZ^T is of dimension " << m_squareMatrixPInv.rows() << "x" << m_squareMatrixPInv.cols() << std::endl;
  
  //-------------------------------------------
  // Third STEP:
  // Invert covariance matrix -> (ZZ^T)^-1
  //-------------------------------------------
  
  std::cout << "  Inverting covariance matrix ... " << std::flush;
  
  if (m_multicollinearityCheckFlag)
  {
    MulticollinearityCheck( m_squareMatrixPInv );
  }
  else
  {
    std::cout << "Ridge parameter: " << m_TikhonovRegularizationParameter << std::endl;
    VnlMatrixType regularizingMatrix( m_squareMatrixPInv.rows(), m_squareMatrixPInv.cols() );
    regularizingMatrix.fill( 0.0 );
    regularizingMatrix.fill_diagonal( m_TikhonovRegularizationParameter );
    m_squareMatrixPInv = m_squareMatrixPInv + regularizingMatrix;
  }
  m_invSquareMatrixPInv = vnl_matrix_inverse<MatrixValueType>( m_squareMatrixPInv );
  
  std::cout << "Ok." << std::endl;

  //-------------------------------------------
  // Fourth STEP:
  // Finally compute estimator B = V(ZZ^T)^-1Z
  //-------------------------------------------
  
  if (m_centeredRegressorMatrix.rows() <= m_centeredRegressorMatrix.cols())
  {
    m_trainedEstimator = m_centeredObservationMatrix * (m_centeredRegressorMatrix.transpose() * m_squareMatrixPInv);
  }
  else
  {
    m_trainedEstimator = (m_centeredObservationMatrix * m_invSquareMatrixPInv) * m_centeredRegressorMatrix.transpose();
  }

  return true;
}

// imiMLRMotionPredictionMethods::CheckForMultiCollinearities
bool icnsMLRMotionPrediction::MulticollinearityCheck( VnlMatrixType& testMatrix )
{
  std::cout <<  "  Checking matrix properties ... " << std::endl;
  MatrixValueType conditionNumber = 0.0;

  // Preparation: Obviously, vnl rank works only for double matrices ...
  vnl_matrix<double> testMatrix_double( testMatrix.rows(), testMatrix.cols() );
  for( unsigned int i = 0; i < testMatrix.rows(); i++ )
  {
    for( unsigned int j = 0; j < testMatrix.cols(); j++ )
    {
      testMatrix_double[i][j] = static_cast<double>( testMatrix[i][j] );
    }
  }

  // FIRST TEST: Compute rank of covariance matrix.
  // If rank is not full, the matrix is (in principle) not invertible.
  // Thus, a workaround is required.

  std::cout <<  "Dimension of matrix: " << testMatrix.rows() << std::endl;
  unsigned int covMatrixRank = vnl_rank( testMatrix_double );
  std::cout <<  "Rank of tested matrix: " << covMatrixRank << std::endl;
  if( covMatrixRank < testMatrix.rows() )
  {
    conditionNumber = -1.0;
    std::cout << "Rank of tested matrix smaller than dimension." << std::endl;
  }

  // SECOND TEST: Compute eigensystem and condition number of covariance matrix.
  // The condition number is defined as sqrt of the ratio of the largest and the smallest
  // eigenvalue of the covariance matrix. A number larger than 30 is usually assumed to
  // indicate multicollinearities (high correlation between at least two regressor variables).

  if( conditionNumber != -1.0 )
  {
    std::cout <<  "Computing eigensystem: " << std::endl << std::endl;
    vnl_symmetric_eigensystem<MatrixValueType> eigSystem( testMatrix );
    int testValue = 0;
    for( unsigned int i = 0; i < testMatrix.rows(); i++ )
    {
      vcl_cerr << "  " << i << "-th eigenvalue: " << eigSystem.get_eigenvalue( i ) << " (" << eigSystem.get_eigenvector( i ) << ")" << std::endl;
      if( eigSystem.get_eigenvalue( i ) < 0 )
      {
        testValue++;
      }
    }

    conditionNumber = sqrt( fabs( eigSystem.get_eigenvalue( testMatrix.rows() - 1 ) ) / fabs( eigSystem.get_eigenvalue( 0 ) ) );
    vcl_cerr << "  --> Condition number: " << conditionNumber << std::endl;

    if( conditionNumber > 30 )
    {
      std::cout << "Condition number > 30." << std::endl;
    }
    else
    {
      std::cout <<  "No multicollinearities in tested matrix!" << std::endl;
      return true;
    }
  }

  // WORKAROUND in case of multicollinearities:
  // Tikhonov regularisation, i.e. AA^T --> AA^T + \lambda*E_n
  // Here we attempt to reduce the condition number to < 30 in an iterative manner.

  std::cout <<  "Performing Tikhonov regularization!" << std::endl;
  MatrixValueType finalTikhonovRegularizationParameter = 0.0;

  //conditionNumber=1;

  while( (conditionNumber > 30))
  {

    //--------------------------------------------------------------------------
    //code max: for more than 500 input points -> due to performance reasons
    //          the eigensystem of the covariance matrix no longer gets checked,
    //          simply adding of 0.2 as regularization parameter

    if( testMatrix.rows() > 500 )
    {
      finalTikhonovRegularizationParameter = 0.2;
      VnlMatrixType regularizingMatrix( testMatrix.rows(), testMatrix.cols() );
      regularizingMatrix.fill( 0.0 );
      regularizingMatrix.fill_diagonal( finalTikhonovRegularizationParameter );
      testMatrix = testMatrix + regularizingMatrix;

      std::cout << " Reg-param set to 0.2 ... " << std::endl;
      std::cout << " Aborting regularization!" << std::endl;
      break;
    }

    //--------------------------------------------------------------------------

    finalTikhonovRegularizationParameter += m_TikhonovRegularizationParameter;

    VnlMatrixType regularizingMatrix( testMatrix.rows(), testMatrix.cols() );
    regularizingMatrix.fill( 0.0 );
    regularizingMatrix.fill_diagonal( m_TikhonovRegularizationParameter );

    testMatrix = testMatrix + regularizingMatrix;

    // Again: computing eigensystem and condition number of updated covariance matrix:
    vnl_symmetric_eigensystem<MatrixValueType> eigSystem( testMatrix );
    conditionNumber = sqrt( fabs( eigSystem.get_eigenvalue( testMatrix.rows() - 1 ) ) / fabs( eigSystem.get_eigenvalue( 0 ) ) );

    if( finalTikhonovRegularizationParameter >= 1.0 )
    {
      std::cout << " Reg-param reached 1.0 ... " << std::endl;
      std::cout << " Condition number: " << conditionNumber << std::endl;
      std::cout << " Aborting regularization!" << std::endl;
      break;
    }
    std::cout << "  Tikhonov regularization factor: " << finalTikhonovRegularizationParameter << "; condition number: " << conditionNumber << std::endl;
  }
  vcl_cerr << "  --> Final Tikhonov regularization parameter: " << finalTikhonovRegularizationParameter << std::endl;
  return true;
}

// imiMLRMotionPredictionMethods::PredictOutput
bool icnsMLRMotionPrediction::PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput )
{
  if (m_trainedEstimator.empty())
  {

    if (m_centeredRegressorMatrix.rows() <= m_centeredRegressorMatrix.cols())
    {
        predictedOutput = m_meanObservationVector + (m_centeredObservationMatrix * (m_centeredRegressorMatrix.transpose() * (m_invSquareMatrixPInv * (regressorMeasurement-m_meanRegressorVector))));
    } else {
        predictedOutput = m_meanObservationVector + ((m_centeredObservationMatrix * m_invSquareMatrixPInv) * (m_centeredRegressorMatrix.transpose()* (regressorMeasurement-m_meanRegressorVector)));
    }

  } else
  {
    predictedOutput = m_meanObservationVector + (m_trainedEstimator * (regressorMeasurement-m_meanRegressorVector));
  }

  std::cout <<  "  finished." << std::endl;
  return true;
}

// ----------------------------
//   Protected / Private methods
// ----------------------------

}// namespace imi
