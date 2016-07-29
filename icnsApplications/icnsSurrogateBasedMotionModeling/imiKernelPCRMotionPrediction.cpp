/** \file imiKernelPCRMotionPrediction.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// Project includes
#include "imiKernelPCRMotionPrediction.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

imiKernelPCRMotionPrediction* imiKernelPCRMotionPrediction::New()
{
  return new imiKernelPCRMotionPrediction();
}

imiKernelPCRMotionPrediction::imiKernelPCRMotionPrediction()
{
  // Set default parameters:
  m_kernelPCA=imiKernelPCA::New();
}

// imiKernelPCRMotionPrediction::~imiKernelPCRMotionPrediction()
imiKernelPCRMotionPrediction::~imiKernelPCRMotionPrediction()
{

}

// ----------------------------
//   Public methods
// ----------------------------

// imiKernelPCRMotionPrediction::TrainKernelPCREstimator
bool imiKernelPCRMotionPrediction::TrainKernelPCREstimator( void )
{
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "  TRAINING OF KernelPCR ESTIMATOR" );
  imiDEBUGINFO( 5, "-----------------------------------" );

  if( m_regressorMatrix.empty() )
  {
    imiERROR( "KernelPCR: Regressor matrix empty! Aborting computation." );
    return false;
  }

  //-------------------------------------------
  // First STEP:
  // Center observation data
  //-------------------------------------------

  m_meanRegressandVector.set_size( m_regressandMatrix.rows(), 1 );
  m_meanRegressandVector.fill( 0.0 );

  if( m_regressandMatrix.empty() )
  {
    imiERROR( "KernelPCR: Observation matrix empty! Aborting computation." );
    return false;
  }

  // Compute mean observation vector:

  imiDEBUGINFO( 5, "  Computing mean regressand vector ... " );
  for( unsigned int i = 0; i < m_regressandMatrix.rows(); i++ )
  {
    m_meanRegressandVector( i, 0 ) = (m_regressandMatrix.get_row( i )).mean();
  }

  // Compute mean regressand matrix entries:

  imiDEBUGINFO( 5, "  Computing regressand matrix with centered entries ... " );
  //m_centeredObservationMatrix = m_observationMatrix;
  m_centeredRegressandMatrix=m_regressandMatrix;
  for( unsigned int i = 0; i < m_regressandMatrix.rows(); i++ )
  {
    m_centeredRegressandMatrix.set_row( i, (m_regressandMatrix.get_row( i ) - m_meanRegressandVector( i, 0 )) );
  }

  //-------------------------------------------
  // Second STEP:
  // Compute KernelPCR-estimator matrix (Due to memory limitations, only the parts of the estimator matrix B=V*Z^T*U_r*(SS^T)^-1*U_r^T are calculated.)
  //-------------------------------------------

  imiDEBUGINFO( 5, "  Computing eigenvectors of the kernel matrix..." << std::endl );

  int tempDebugLevel=this->GetGlobalDebugLevel();

  this->SetGlobalDebugLevel(0);
  m_kernelPCA->SetCentered(true);
  m_kernelPCA->SetInputMatrix(m_regressorMatrix);
  m_kernelPCA->ComputePrincipalComponents();
  this->SetGlobalDebugLevel(tempDebugLevel);

  VnlDiagMatrixType eigenvalues=m_kernelPCA->GetEigenvalues();

  // Calculate Sum of Eigenvalues and relative Eigenvalue
  double sumOfEigenvalues = 0;
  VnlVectorType RelativeEigenValues( eigenvalues.size() );
  for( unsigned int i = 0; i < eigenvalues.size(); i++ )
  {
    sumOfEigenvalues += eigenvalues[i];
  }
  for( unsigned int i = 0; i < eigenvalues.size(); i++ )
  {
    RelativeEigenValues[i] = eigenvalues[i] / sumOfEigenvalues;
  }

  imiDEBUGINFO( 5, "Eigenvalues lambda are : [ " );
  if( this->CheckGlobalDebugLevel( 5 ) )
  {
    for( unsigned int i = 0; i < eigenvalues.size() - 1; i++ )
    {
      std::cout << eigenvalues[i] << "(" << RelativeEigenValues[i] << "%) , ";
    }
    std::cout << eigenvalues[eigenvalues.size() - 1] << "(" << RelativeEigenValues[eigenvalues.size() - 1] << "%)]\n";
  }

  //determine the number of modes needed to describe x% of the variability
  //m_variabilityThreshold = 0.95;
  //m_numOfComponentsUsed = 0;
  double temp = 0;

  if (m_numOfComponentsUsed < 0)
  {
    m_numOfComponentsUsed=0;
    for( unsigned int i = 0; (i < m_invEigenvalues.size() - 1) && temp < m_variabilityThreshold; i++, m_numOfComponentsUsed++ )
    {
      temp += RelativeEigenValues[i];
    }
  }

  //////////////////////////////////////////////////////////////////

  imiDEBUGINFO( 5, "Get first "<<m_numOfComponentsUsed<<" Eigenvectors ... " );
  m_principalComponents = VnlMatrixType( m_numOfComponentsUsed, m_regressorMatrix.cols());

  for( unsigned int j = 0; j < m_principalComponents.cols(); j++ )
  {
    VnlVectorType temp = m_regressorMatrix.get_column(j);
    VnlVectorType projection=m_kernelPCA->ProjectDataPointOntoEigenvectors(temp);
    //std::cout<<"proj "<<j<<":"<<projection<<std::endl;
    for( unsigned int i = 0; i < m_principalComponents.rows(); i++ )
    {
      m_principalComponents[i][j] = projection[i];
    }
  }

  imiDEBUGINFO( 5, "Get first "<<m_numOfComponentsUsed<<" Eigenvalues (SS^T) and invert the resulting diagonal matrix... " );
  m_invEigenvalues.set_size( m_numOfComponentsUsed );

  for( int j = 0; j < m_numOfComponentsUsed; j++ )
  {
    m_invEigenvalues[j] = 1. / eigenvalues[j];
  }

  imiDEBUGINFO( 5, "  Training finished." );

  return true;
}

// imiKernelPCRMotionPrediction::PredictOutput
bool imiKernelPCRMotionPrediction::PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput )
{
  imiDEBUGINFO( 5, "  Predicting output for given regressor measurement ... " );

  if( m_meanRegressandVector.empty() || m_centeredRegressandMatrix.empty() || m_principalComponents.empty() )
  {
    imiERROR( "Prediction not possible: KernelPCR estimator not trained. Call TrainKernelPCREstimator() before trying to predict." );
    return false;
  }


    VnlVectorType temp=regressorMeasurement.get_column(0);
    VnlVectorType projection=m_kernelPCA->ProjectDataPointOntoEigenvectors(temp);
    VnlMatrixType projectionMatrix;
    projectionMatrix.set_size(m_numOfComponentsUsed,1);
    projectionMatrix.set_column(0,projection);
    for( int i = 0; i < m_numOfComponentsUsed; i++ )
    {
      projectionMatrix[i][0] = projection[i];
    }

/*  std::cout<<"m_meanRegressandVector: ("<<m_meanRegressandVector.rows()<<","<<m_meanRegressandVector.cols()<<")"<<std::endl;
  std::cout<<"m_centeredRegressandMatrix: ("<<m_centeredRegressandMatrix.rows()<<","<<m_centeredRegressandMatrix.cols()<<")"<<std::endl;
  std::cout<<"m_principalComponents: ("<<m_principalComponents.rows()<<","<<m_principalComponents.cols()<<")"<<std::endl;
  std::cout<<"m_invEigenvalues: ("<<m_invEigenvalues.rows()<<","<<m_invEigenvalues.cols()<<")"<<std::endl;
  std::cout<<"projectionMatrix: ("<<projectionMatrix.rows()<<","<<projectionMatrix.cols()<<")"<<std::endl;*/
  /*std::cout<<"m_principalComponents:"<<m_principalComponents<<std::endl;
  std::cout<<"m_principalComponents.T:"<<m_principalComponents.transpose()<<std::endl;
  std::cout<<"m_invEigenvalues:"<<m_invEigenvalues<<std::endl;
  std::cout<<"projectionMatrix:"<<projectionMatrix<<std::endl;*/


  predictedOutput = m_meanRegressandVector
      + (m_centeredRegressandMatrix * m_principalComponents.transpose() * m_invEigenvalues * projectionMatrix );

  imiDEBUGINFO( 5, "  finished." );
  return true;
}

// ----------------------------
//   Protected / Private methods
// ----------------------------

}// namespace imi
