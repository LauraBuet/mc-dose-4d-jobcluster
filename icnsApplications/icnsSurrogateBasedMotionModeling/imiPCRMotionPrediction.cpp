/** \file imiPCRMotionPrediction.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// Project includes
#include "imiPCRMotionPrediction.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

imiPCRMotionPrediction* imiPCRMotionPrediction::New()
{
  return new imiPCRMotionPrediction();
}

imiPCRMotionPrediction::imiPCRMotionPrediction()
{
  // Set default parameters:
}

// imiPCRMotionPrediction::~imiPCRMotionPrediction()
imiPCRMotionPrediction::~imiPCRMotionPrediction()
{
}

// ----------------------------
//   Public methods
// ----------------------------

// imiPCRMotionPrediction::TrainPCREstimator
bool imiPCRMotionPrediction::TrainEstimator( void )
{
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "  TRAINING OF PCR ESTIMATOR" );
  imiDEBUGINFO( 5, "-----------------------------------" );

  //-------------------------------------------
  // FIRST STEP:
  // Center regressor data
  //-------------------------------------------

  m_meanRegressorVector.set_size( m_regressorMatrix.rows(), 1 );
  m_meanRegressorVector.fill( 0.0 );

  if( m_regressorMatrix.empty() )
  {
    imiERROR( "PCR: Regressor matrix empty! Aborting computation." );
    return false;
  }

  // Compute mean regressor vector:

  imiDEBUGINFO( 5, "  Computing mean regressor vector ... " );
  for( unsigned int i = 0; i < m_regressorMatrix.rows(); i++ )
  {
    m_meanRegressorVector( i, 0 ) = (m_regressorMatrix.get_row( i )).mean();
  }

  // Compute mean regressor matrix entries:

  imiDEBUGINFO( 5, "  Computing regressor matrix with centered entries ... " );
  m_centeredRegressorMatrix = m_regressorMatrix;
  for( unsigned int i = 0; i < m_regressorMatrix.rows(); i++ )
  {
    m_centeredRegressorMatrix.set_row( i, (m_regressorMatrix.get_row( i ) - m_meanRegressorVector( i, 0 )) );
  }

  //-------------------------------------------
  // SECOND STEP:
  // Center observation data
  //-------------------------------------------

  m_meanRegressandVector.set_size( m_regressandMatrix.rows(), 1 );
  m_meanRegressandVector.fill( 0.0 );

  if( m_regressandMatrix.empty() )
  {
    imiERROR( "PCR: Observation matrix empty! Aborting computation." );
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
  // THIRD STEP:
  // Compute PCR-estimator matrix (Due to memory limitations, only the parts of the estimator matrix B=V*Z^T*U_r*(SS^T)^-1*U_r^T are calculated.)
  //-------------------------------------------

  imiDEBUGINFO( 5, "  Computing eigenvectors of the regressor covariance matrix ZZ^T (SVD of Z)..." << std::endl );

  vnl_svd<MatrixValueType> svdSolver( m_centeredRegressorMatrix );

  VnlMatrixType VMatrix = svdSolver.V();
  VnlMatrixType UMatrix = svdSolver.U();
  VnlDiagMatrixType SigmaMatrix = svdSolver.W();

  imiDEBUGINFO( 5, " Check if VV^T = Id ..." );
  // set tolerance dependend on the matrix-value type
  const imiRealType tolerance = 1e-5;
  VnlMatrixType Ident = VMatrix * VMatrix.transpose();
  for( unsigned int i = 0; i < Ident.rows(); i++ )
    for( unsigned int j = 0; j < Ident.cols(); j++ )
    {
      if( (i != j && fabs( Ident[i][j] ) > tolerance) || (i == j && fabs( Ident[i][j] - 1.0 ) > tolerance) )
      {
        imiWARNING( "VV^T["<<i<<","<<j<<"]="<<Ident[i][j]<<" is not identity! " );
      }
    }

  imiDEBUGINFO( 5, " Check if U^T*U = Id ..." );
  Ident = UMatrix.transpose() * UMatrix;
  for( unsigned int i = 0; i < Ident.rows(); i++ )
    for( unsigned int j = 0; j < Ident.cols(); j++ )
    {
      if( (i != j && fabs( Ident[i][j] ) > tolerance) || (i == j && fabs( Ident[i][j] - 1.0 ) > tolerance) )
      {
        imiWARNING( "U^T*U["<<i<<","<<j<<"]="<<Ident[i][j]<<" is not identity! " );
      }
    }

  SigmaMatrix = SigmaMatrix * SigmaMatrix;

  // Calculate Sum of Eigenvalues and relative Eigenvalue
  double sumOfEigenvalues = 0;
  VnlVectorType RelativeEigenValues( SigmaMatrix.size() );
  for( unsigned int i = 0; i < SigmaMatrix.size(); i++ )
  {
    sumOfEigenvalues += SigmaMatrix[i];
  }
  for( unsigned int i = 0; i < SigmaMatrix.size(); i++ )
  {
    RelativeEigenValues[i] = SigmaMatrix[i] / sumOfEigenvalues;
  }

  imiDEBUGINFO( 5, "Eigenvalues lambda are : [ " );
  if( this->CheckGlobalDebugLevel( 5 ) )
  {
    for( unsigned int i = 0; i < SigmaMatrix.size() - 1; i++ )
    {
      std::cout << SigmaMatrix[i] << "(" << RelativeEigenValues[i] << "%) , ";
    }
    std::cout << SigmaMatrix[SigmaMatrix.size() - 1] << "(" << RelativeEigenValues[SigmaMatrix.size() - 1] << "%)]\n";
  }

  //determine the number of modes needed to describe x% of the variability
  m_variabilityThreshold = 0.95;
  //m_numOfComponents = 0;
//  double temp = 0;

  /*for( unsigned int i = 0; (i < SigmaMatrix.size() - 1) && temp < m_variabilityThreshold; i++, m_numOfComponents++ )
  {
    temp += RelativeEigenValues[i];
  }*/

  //m_numOfModesUsed=m_centeredRegressandMatrix.cols()-1;
  //m_numOfComponents=2;
  //////////////////////////////////////////////////////////////////

  imiDEBUGINFO( 5, "Get first "<<m_numOfComponents<<" Eigenvektors (cols of U) ... " );
  m_regressorComponents = VnlMatrixType( m_centeredRegressorMatrix.rows(), m_numOfComponents );

  for( unsigned int j = 0; j < m_regressorComponents.cols(); j++ )
  {
    for( unsigned int i = 0; i < m_regressorComponents.rows(); i++ )
    {
      m_regressorComponents[i][j] = svdSolver.U( i, j );
    }
  }

  imiDEBUGINFO( 5, "Get first "<<m_numOfComponents<<" Eigenvalues (SS^T) and invert the resulting diagonal matrix... " );
  m_invEigenvalues.set_size( m_numOfComponents );

  for( unsigned int j = 0; j < m_numOfComponents; j++ )
  {
    m_invEigenvalues[j] = 1. / SigmaMatrix[j];
  }

  imiDEBUGINFO( 5, "  Training finished." );

  return true;
}

// imiPCRMotionPrediction::PredictOutput
bool imiPCRMotionPrediction::PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput )
{
  imiDEBUGINFO( 5, "  Predicting output for given regressor measurement ... " );

  if( m_meanRegressandVector.empty() || m_meanRegressorVector.empty() || m_centeredRegressandMatrix.empty() || m_centeredRegressorMatrix.empty() || m_regressorComponents.empty() )
  {
    imiERROR( "Prediction not possible: PCR estimator not trained. Call TrainPCREstimator() before trying to predict." );
    return false;
  }

  /*std::cout<<"componentes: "<<m_centeredRegressorMatrix.transpose() *m_regressorComponents<<std::endl;
  std::cout<<"projection: "<<(m_regressorComponents.transpose() * (regressorMeasurement - m_meanRegressorVector))<<std::endl;*/

  VnlMatrixType predRegressorComponents;
  VnlDiagMatrixType predInvEigenvalues;
  predRegressorComponents.set_size(m_regressorComponents.rows(),m_numOfComponents);
  predInvEigenvalues.set_size(m_numOfComponents);

  for( unsigned int j = 0; j < m_numOfComponents; j++ )
  {
   predRegressorComponents.set_column(j,m_regressorComponents.get_column(j));
   predInvEigenvalues[j] = m_invEigenvalues[j];
  }

  predictedOutput = m_meanRegressandVector
      + (m_centeredRegressandMatrix
          * (m_centeredRegressorMatrix.transpose() * (predRegressorComponents * (predInvEigenvalues * (predRegressorComponents.transpose() * (regressorMeasurement - m_meanRegressorVector))))));

  imiDEBUGINFO( 5, "  finished." );
  return true;
}

// ----------------------------
//   Protected / Private methods
// ----------------------------

}// namespace imi
