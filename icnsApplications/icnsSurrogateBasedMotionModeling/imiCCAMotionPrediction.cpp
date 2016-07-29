/** \file imiCCAMotionPrediction.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/
// Project includes
#include "imiCCAMotionPrediction.h"

#include "imiLinearKernel.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

imiCCAMotionPrediction* imiCCAMotionPrediction::New()
{
  return new imiCCAMotionPrediction();
}

imiCCAMotionPrediction::imiCCAMotionPrediction()
{
  // Set default parameters:
  m_numOfComponents = 2;
  m_useBothBases = false;

  m_maxIters = 500;
  m_tolerance = 1e-6;

  m_TikhonovRegularizationParameter = 0.001;

  m_regressorPCAvariabilityThreshold = 1.00;
  m_regressandPCAvariabilityThreshold = 1.00;

}

// imiCCAMotionPrediction::~imiCCAMotionPrediction()
imiCCAMotionPrediction::~imiCCAMotionPrediction()
{
}

// ----------------------------
//   Public methods
// ----------------------------

// imiCCAMotionPrediction::TrainCCAEstimator
bool imiCCAMotionPrediction::TrainEstimator( void )
{
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "  TRAINING OF CCA ESTIMATOR" );
  imiDEBUGINFO( 5, "-----------------------------------" );

  //-------------------------------------------
  // FIRST STEP:
  // Perform PCA on regressor data
  //-------------------------------------------

  imiDEBUGINFO( 5, " Computing SVD of mean centered regressor data matrix... " );

  if( m_regressorMatrix.empty() )
  {
    imiERROR( "CCA: Regressor matrix empty! Aborting computation." );
    return false;
  }

  //START OF REGRESSOR PCA

 /* m_regressorPCAmeanVector.set_size( m_regressorMatrix.rows(), 1 );

  for( unsigned int i = 0; i < m_regressorMatrix.rows(); i++ )
  {
    m_regressorPCAmeanVector( i, 0 ) = (m_regressorMatrix.get_row( i )).mean();
  }

  // Compute mean regressand matrix entries:

  imiDEBUGINFO( 5, "  Computing regressand matrix with centered entries ... " );
  VnlMatrixType regressorPCAcenteredMatrix = m_regressorMatrix;
  for( unsigned int i = 0; i < m_regressorMatrix.rows(); i++ )
  {
    regressorPCAcenteredMatrix.set_row( i, (m_regressorMatrix.get_row( i ) - m_regressorPCAmeanVector( i, 0 )) );
  }

  imiDEBUGINFO( 5, "  Computing eigenvectors of the regressor covariance matrix ZZ^T (SVD of Z)..." << std::endl );

  vnl_svd<MatrixValueType> svdSolver( regressorPCAcenteredMatrix );

  VnlMatrixType VMatrix = svdSolver.V();
  VnlMatrixType UMatrix = svdSolver.U();
  VnlDiagMatrixType SigmaMatrix = svdSolver.W();

  //std::cout<<"W"<<svdSolver.W()<<std::endl;

  SigmaMatrix = SigmaMatrix * SigmaMatrix;

  // Calculate sum of Eigenvalues and relative Eigenvalues
  double sumOfEigenvalues = 0;
  VnlVectorType RelativeEigenValues( SigmaMatrix.size() );
  for( unsigned int i = 0; i < SigmaMatrix.size(); i++ )
  {
    sumOfEigenvalues += SigmaMatrix[i];
    std::cout<<"eigenvals: "<<SigmaMatrix[i]<<std::endl;
  }
  for( unsigned int i = 0; i < SigmaMatrix.size(); i++ )
  {
    RelativeEigenValues[i] = SigmaMatrix[i] / sumOfEigenvalues;
  }

  //determine the number of modes needed to describe x% of the variability
  unsigned int regressorPCAnumOfModesUsed = 0;
  double temp = 0;

  for( unsigned int i = 0; (i < SigmaMatrix.size() - 1) && temp < m_regressorPCAvariabilityThreshold; i++, regressorPCAnumOfModesUsed++ )
  {
    temp += RelativeEigenValues[i];
  }

//  regressorPCAnumOfModesUsed=regressorPCAcenteredMatrix.cols()-1;
  //////////////////////////////////////////////////////////////////

  imiDEBUGINFO( 5, "Get first "<<regressorPCAnumOfModesUsed<<" Eigenvektors (cols of U) ... " );
  m_regressorPCAComponents = VnlMatrixType( regressorPCAcenteredMatrix.rows(), regressorPCAnumOfModesUsed );

  for( unsigned int j = 0; j < m_regressorPCAComponents.cols(); j++ )
  {
    for( unsigned int i = 0; i < m_regressorPCAComponents.rows(); i++ )
    {
      m_regressorPCAComponents[i][j] = svdSolver.U( i, j );
    }
  }

  //Projection of regressor data matrix on pca "subspace"
  m_regressorMatrix = m_regressorPCAComponents.transpose() * regressorPCAcenteredMatrix;

  //END OF REGRESSOR PCA

  //START OF REGRESSAND PCA

  m_regressandPCAmeanVector.set_size( m_regressandMatrix.rows(), 1 );

  for( unsigned int i = 0; i < m_regressandMatrix.rows(); i++ )
  {
    m_regressandPCAmeanVector( i, 0 ) = (m_regressandMatrix.get_row( i )).mean();
  }

  // Compute mean regressand matrix entries:

  imiDEBUGINFO( 5, "  Computing regressand matrix with centered entries ... " );
  VnlMatrixType regressandPCAcenteredMatrix = m_regressandMatrix;
  for( unsigned int i = 0; i < m_regressandMatrix.rows(); i++ )
  {
    regressandPCAcenteredMatrix.set_row( i, (m_regressandMatrix.get_row( i ) - m_regressandPCAmeanVector( i, 0 )) );
  }

  imiDEBUGINFO( 5, "  Computing eigenvectors of the implicit regressand covariance matrix V^T*V..." << std::endl );

  vnl_svd<MatrixValueType> svdSolver2( regressandPCAcenteredMatrix.transpose()*regressandPCAcenteredMatrix );

  VnlMatrixType VMatrix2 = svdSolver2.V();
  VnlMatrixType UMatrix2 = svdSolver2.U();
  VnlDiagMatrixType SigmaMatrix2 = svdSolver2.W();



  // Calculate sum of Eigenvalues and relative Eigenvalues
  double sumOfEigenvalues2 = 0;
  VnlVectorType RelativeEigenValues2( SigmaMatrix2.size() );
  for( unsigned int i = 0; i < SigmaMatrix2.size(); i++ )
  {
    sumOfEigenvalues2 += SigmaMatrix2[i];
    std::cout<<"eigenvals: "<<SigmaMatrix2[i]<<std::endl;
  }
  for( unsigned int i = 0; i < SigmaMatrix2.size(); i++ )
  {
    RelativeEigenValues2[i] = SigmaMatrix2[i] / sumOfEigenvalues2;
  }

  //determine the number of modes needed to describe x% of the variability
  unsigned int regressandPCAnumOfModesUsed = 0;
  temp = 0;

  for( unsigned int i = 0; (i < SigmaMatrix2.size() - 1) && temp < m_regressandPCAvariabilityThreshold; i++, regressandPCAnumOfModesUsed++ )
  {
    temp += RelativeEigenValues2[i];
  }

//  regressorPCAnumOfModesUsed=regressorPCAcenteredMatrix.cols()-1;
  //////////////////////////////////////////////////////////////////

  imiDEBUGINFO( 5, "Get first "<<regressandPCAnumOfModesUsed<<" Eigenvektors (cols of U) ... " );
  m_regressandPCAComponents = VnlMatrixType( regressandPCAcenteredMatrix.rows(), regressandPCAnumOfModesUsed );

  std::cout<<"TEST1"<<std::endl;

  VnlMatrixType newEigenvectors=regressandPCAcenteredMatrix*UMatrix2;

  std::cout<<"TEST"<<std::endl;

  for( unsigned int j = 0; j < m_regressandPCAComponents.cols(); j++ )
  {
    for( unsigned int i = 0; i < m_regressandPCAComponents.rows(); i++ )
    {
      m_regressandPCAComponents[i][j] = newEigenvectors[ i][j]/sqrt(SigmaMatrix2[j]);
    }
  }

  //Projection of regressand data matrix on pca "subspace"
  m_regressandMatrix = m_regressandPCAComponents.transpose() * regressandPCAcenteredMatrix;*/

 /* imiDEBUGINFO( 5, "  Computing eigenvectors of the regressand covariance matrix VV^T (SVD of V)..." << std::endl );

  m_regressandPCA->SetInputMatrix(m_regressandMatrix);
  m_regressandPCA->ComputePrincipalComponents();

  imiDEBUGINFO( 5, "Get first "<<m_regressandPCA->NumOfCompsNeededForVarThresh()<<" Eigenvektors (cols of U) ... " );

  VnlMatrixType regressandPCAProjections;
  regressandPCAProjections.set_size(m_regressandPCA->NumOfCompsNeededForVarThresh(),m_regressandMatrix.cols());

  for (unsigned int i=0; i<m_regressandMatrix.cols();i++)
  {
    VnlVectorType tempVector;
    tempVector=m_regressandMatrix.get_column(i);
    tempVector=m_regressandPCA->ProjectDataPointOntoEigenvectors(tempVector);
    regressandPCAProjections.set_column(i,tempVector);
  }

  m_regressandMatrix=regressandPCAProjections;

  //END OF REGRESSAND PCA*/



  /*m_regressorMatrix.set_size(3,10);
   imiRealType zarray[3*10]={0.7517,    0.0172,    0.5387,    0.0945,    0.2943,    0.0682,    0.6513,    0.8169,    0.2124,    0.9564,
   0.3684,    0.8291,    0.6505,    0.8776,    0.1799,    0.5811,    0.8646,    0.5289,    0.5433,    0.4445,
   0.9418,    0.6266,    0.7266,    0.0144,    0.9263,    0.6372,    0.0560,    0.6944,    0.7025,    0.0854};
   m_regressorMatrix=VnlMatrixType(zarray,3,10);


   m_regressandMatrix.set_size(5,10);
   imiRealType varray[5*10]={0.8900,    0.2284,    0.8801,    0.1139,    0.3257,    0.5999,    0.1080,    0.7009,    0.9585,    0.7409,
   0.3302,    0.6520,    0.4443,    0.9786,    0.6302,    0.4484,    0.4599,    0.8722,    0.7900,    0.5068,
   0.2297,    0.0662,    0.7559,    0.8486,    0.2303,    0.0354,    0.4509,    0.0522,    0.4519,    0.1999,
   0.1139,    0.2754,    0.6033,    0.0506,    0.5799,    0.5138,    0.5511,    0.2197,    0.3334,    0.4272,
   0.3109,    0.2818,    0.7833,    0.4662,    0.6032,    0.4077,    0.8054,    0.4596,    0.0591,    0.1687};
   m_regressandMatrix=VnlMatrixType(varray,5,10);*/

  m_meanRegressorVector.set_size( m_regressorMatrix.rows(), 1 );
  m_meanRegressorVector.fill( 0.0 );

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
    imiERROR( "CCA: Observation matrix empty! Aborting computation." );
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
  m_centeredRegressandMatrix = m_regressandMatrix;
  for( unsigned int i = 0; i < m_regressandMatrix.rows(); i++ )
  {
    m_centeredRegressandMatrix.set_row( i, (m_regressandMatrix.get_row( i ) - m_meanRegressandVector( i, 0 )) );
  }

  //-------------------------------------------
  // FOURTH STEP:
  // Compute CCA-estimator matrix (Due to memory limitations, we perform a PCA on the regressor and regressand data first.)
  //-------------------------------------------

  imiDEBUGINFO( 5, " Computing eigenvectors of the CCA eigenproblem... " );

  vnl_matrix<double> centeredRegressorMatrix_double( m_centeredRegressorMatrix.rows(), m_centeredRegressorMatrix.cols() );
  for( unsigned int i = 0; i < m_centeredRegressorMatrix.rows(); i++ )
  {
    for( unsigned int j = 0; j < m_centeredRegressorMatrix.cols(); j++ )
    {
      centeredRegressorMatrix_double[i][j] = static_cast<double>( m_centeredRegressorMatrix[i][j] );
    }
  }
  vnl_matrix<double> centeredRegressandMatrix_double( m_centeredRegressandMatrix.rows(), m_centeredRegressandMatrix.cols() );
  for( unsigned int i = 0; i < m_centeredRegressandMatrix.rows(); i++ )
  {
    for( unsigned int j = 0; j < m_centeredRegressandMatrix.cols(); j++ )
    {
      centeredRegressandMatrix_double[i][j] = static_cast<double>( m_centeredRegressandMatrix[i][j] );
    }
  }
  vnl_matrix<double> regularizingMatrixX( m_centeredRegressorMatrix.rows(), m_centeredRegressorMatrix.cols() );
  regularizingMatrixX.fill( 0.0 );
  regularizingMatrixX.fill_diagonal( m_tolerance );
  vnl_matrix<double> regularizingMatrixY( m_centeredRegressandMatrix.rows(), m_centeredRegressandMatrix.cols() );
  regularizingMatrixY.fill( 0.0 );
  regularizingMatrixY.fill_diagonal( m_tolerance );
  vnl_qr<double> qr_vv( (centeredRegressandMatrix_double * centeredRegressandMatrix_double.transpose()) + regularizingMatrixY );
  vnl_matrix<double> VVinv = qr_vv.inverse(); //* centeredRegressandMatrix_double.transpose();

  vnl_qr<double> qr_zz( (centeredRegressorMatrix_double * centeredRegressorMatrix_double.transpose()) + regularizingMatrixY );
  vnl_matrix<double> ZZinv = qr_zz.inverse();

  imiDEBUGINFO( 5, "Computing regressor eigensystem... " << std::endl );

  //computing new regressor basis
  vnl_matrix<double> A = centeredRegressorMatrix_double * centeredRegressandMatrix_double.transpose() * VVinv * centeredRegressandMatrix_double * centeredRegressorMatrix_double.transpose();
  vnl_matrix<double> B = ((centeredRegressorMatrix_double) * (centeredRegressorMatrix_double).transpose()) + regularizingMatrixX;

  /*std::cout<<"A: "<<A<<std::endl;
   std::cout<<"B: "<<B<<std::endl;*/

  vnl_generalized_eigensystem eigensystemRegressor( A, B );

  /*std::cout<<"V: "<<eigensystem.V<<std::endl;
   std::cout<<"D: "<<eigensystem.D<<std::endl;*/

  //std::cout<<"X: "<<m_regressorMatrix<<std::endl;
  //std::cout<<"Y: "<<m_regressandMatrix<<std::endl;
  /*std::cout<<"meanX: "<<m_meanRegressorVector<<std::endl;
   std::cout<<"meanY: "<<m_meanRegressandVector<<std::endl;*/

  //imiDEBUGINFO( 5, "\t ...done. Eigenvalues ("<<eigensystem.D.size()<<"): ");
  m_numOfComponents = std::min<unsigned int>( m_numOfComponents, std::min<unsigned int>( m_centeredRegressandMatrix.rows(), m_centeredRegressorMatrix.rows() ) );

  imiDEBUGINFO( 5, "\t Using "<<m_numOfComponents<<" components (max. possible "<<std::min<unsigned int>( m_centeredRegressandMatrix.rows(),m_centeredRegressorMatrix.rows())<<" )" );

  m_regressandComponents.set_size( m_centeredRegressandMatrix.rows(), m_numOfComponents );
  m_regressorComponents.set_size( m_centeredRegressorMatrix.rows(), m_numOfComponents );
  m_regressandEigenvalues.set_size( m_numOfComponents );
  m_regressorEigenvalues.set_size( m_numOfComponents );

  for( unsigned int i = 0; i < m_numOfComponents; i++ )
  {
    unsigned int evIndex = std::max<unsigned int>( eigensystemRegressor.D.size() - 1 - i, 0 );
    vnl_vector<double> evTemp;

    evTemp = eigensystemRegressor.V.get_column( evIndex ).normalize();
    for( unsigned int j = 0; j < evTemp.size(); j++ )
    {
      m_regressorComponents[j][i] = static_cast<float>( evTemp[j] );
    }

    m_regressorEigenvalues[i] = static_cast<float>( eigensystemRegressor.D[evIndex] );

    imiDEBUGINFO( 5, "\t\t "<<i+1<<": "<<m_regressorEigenvalues[i] <<"  correlation: "<<sqrt(m_regressorEigenvalues[i]) );

  }

  if( m_useBothBases )
  {
    imiDEBUGINFO( 5, "Computing regressand eigensystem... " << std::endl );

    //computing new regressor basis
    A = centeredRegressandMatrix_double * centeredRegressorMatrix_double.transpose() * ZZinv * centeredRegressorMatrix_double * centeredRegressandMatrix_double.transpose();
    B = ((centeredRegressandMatrix_double) * (centeredRegressandMatrix_double).transpose()) + regularizingMatrixX;

    vnl_generalized_eigensystem eigensystemRegressand( A, B );

    imiDEBUGINFO( 5, "\t Using "<<m_numOfComponents<<" components (max. possible "<<std::min<unsigned int>( m_centeredRegressandMatrix.rows(),m_centeredRegressorMatrix.rows())<<" )" );

    for( unsigned int i = 0; i < m_numOfComponents; i++ )
    {
      unsigned int evIndex = std::max<unsigned int>( eigensystemRegressand.D.size() - 1 - i, 0 );
      vnl_vector<double> evTemp;

      evTemp = eigensystemRegressand.V.get_column( evIndex ).normalize();
      for( unsigned int j = 0; j < evTemp.size(); j++ )
      {
        m_regressandComponents[j][i] = static_cast<float>( evTemp[j] );
      }

      m_regressandEigenvalues[i] = static_cast<float>( eigensystemRegressand.D[evIndex] );

      imiDEBUGINFO( 5, "\t\t "<<i+1<<": "<<m_regressandEigenvalues[i] <<"  correlation: "<<sqrt(m_regressandEigenvalues[i]) );

    }

  }

  imiDEBUGINFO( 5, "\t Checking if regressor weights are not pairwise orthogonal..." );
  bool check = true;
  for( unsigned int i = 0; i < m_regressorComponents.cols() - 1; i++ )
  {
    for( unsigned int j = i + 1; j < m_regressorComponents.cols(); j++ )
    {
      if( dot_product( m_regressorComponents.get_column( i ), m_regressorComponents.get_column( j ) ) > m_tolerance )
      {
        check = false;
        //std::cout << "i: " << i << " j: " << j << "dot: " << dot_product( m_regressorComponents.get_column( i ), m_regressorComponents.get_column( j ) ) << std::endl;
      }
    }
  }
  if( !check )
  {
    imiDEBUGINFO( 5, "\t Succesful." );
  }
  else
  {
    imiWARNING( "Regressor weights are not orthogonal!" )
  }

  if( m_useBothBases )
  {
    imiDEBUGINFO( 5, "\t Checking if regressand weights are pairwise orthogonal..." );
    check = true;
    for( unsigned int i = 0; i < m_regressandComponents.cols() - 1; i++ )
    {
      for( unsigned int j = i + 1; j < m_regressandComponents.cols(); j++ )
      {
        //std::cout<<"i: "<<i<<" j: "<<j<<"dot: "<<dot_product(m_regressandComponents.get_column(i),m_regressandComponents.get_column(j))<<std::endl;
        if( dot_product( m_regressandComponents.get_column( i ), m_regressandComponents.get_column( j ) ) > m_tolerance )
        {
          check = false;
          //std::cout << "i: " << i << " j: " << j << "dot: " << dot_product( m_regressandComponents.get_column( i ), m_regressandComponents.get_column( j ) ) << std::endl;
        }
      }
    }
    if( check )
    {
      imiDEBUGINFO( 5, "\t Succesful." );
    }
    else
    {
      imiWARNING( "Regressand weights are not pairwise orthogonal!" )
    }
  }

  return true;
}

// imiCCAMotionPrediction::PredictOutput
bool imiCCAMotionPrediction::PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput )
{
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "           PREDICTION              " );
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "  Predicting output for given regressor measurement ... " );

  if( m_meanRegressandVector.empty() || m_meanRegressorVector.empty() || m_centeredRegressandMatrix.empty() || m_centeredRegressorMatrix.empty() )
  {
    imiERROR( "Prediction not possible: CCA estimator not trained. Call TrainCCAEstimator() before trying to predict." );
    return false;
  }

  VnlMatrixType predRegressorComponents;
  VnlMatrixType predRegressandComponents;
  predRegressorComponents.set_size(m_regressorComponents.rows(),m_numOfComponents);
  predRegressandComponents.set_size(m_regressandComponents.rows(),m_numOfComponents);

  for( unsigned int j = 0; j < m_numOfComponents; j++ )
  {
   predRegressorComponents.set_column(j,m_regressorComponents.get_column(j));
   predRegressandComponents.set_column(j,m_regressandComponents.get_column(j));

  }

  VnlMatrixType scores = m_centeredRegressorMatrix.transpose() * predRegressorComponents;
  scores = scores.transpose() * scores;
//  this->MulticollinearityCheck( scores );

  vnl_svd<imiRealType> svd_scores( scores );

  VnlMatrixType scoresInv = svd_scores.inverse();

  if( m_useBothBases )
  {
 /*   predictedOutput = m_meanRegressandVector
        + (predRegressandComponents*(predRegressandComponents.transpose() * m_centeredRegressandMatrix
            * (m_centeredRegressorMatrix.transpose()
                * (predRegressorComponents
                    * (scoresInv * (predRegressorComponents.transpose() * ((m_regressorPCAComponents.transpose() * (regressorMeasurement - m_regressorPCAmeanVector)) - m_meanRegressorVector)))))));*/
    predictedOutput = m_meanRegressandVector
            + (predRegressandComponents
                * (predRegressandComponents.transpose() * m_centeredRegressandMatrix
                    * (m_centeredRegressorMatrix.transpose() * (predRegressorComponents * (scoresInv * (predRegressorComponents.transpose() * (regressorMeasurement - m_meanRegressorVector)))))));
  }
  else
  {
/*    predictedOutput = m_meanRegressandVector
        + (m_centeredRegressandMatrix
            * (m_centeredRegressorMatrix.transpose()
                * (predRegressorComponents
                    * (scoresInv * (predRegressorComponents.transpose() * ((m_regressorPCAComponents.transpose() * (regressorMeasurement - m_regressorPCAmeanVector)) - m_meanRegressorVector))))));*/
    predictedOutput = m_meanRegressandVector
            + (m_centeredRegressandMatrix
                * (m_centeredRegressorMatrix.transpose() * (predRegressorComponents * (scoresInv * (predRegressorComponents.transpose() * (regressorMeasurement - m_meanRegressorVector))))));
  }

//  predictedOutput=m_regressandPCAmeanVector+m_regressandPCAComponents*predictedOutput;
  //VnlVectorType tempVector=predictedOutput.get_column(0);
  //predictedOutput.set_column(0,m_regressandPCA->ComputePreImageOfProjection(tempVector));

  /*std::cout<<"Prediction"<<predictedOutput<<std::endl;
   std::cout<<"m_regressandMatrix"<<m_regressandMatrix<<std::endl;
   std::cout<<"m_regressorMatrix"<<m_regressorMatrix<<std::endl;*/
  imiDEBUGINFO( 5, "  finished." );
  return true;
}

// ----------------------------
//   Protected / Private methods
// ----------------------------


} // namespace imi

#ifdef FALSE
// imiCCAMotionPrediction::TrainCCAEstimator
bool imiCCAMotionPrediction::TrainEstimator( void )
{
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "  TRAINING OF CCA ESTIMATOR" );
  imiDEBUGINFO( 5, "-----------------------------------" );

  /*double testdataV[36] = {
   30.0000,   -3.4273,   13.9254,   13.7049,   -2.4446,   20.2380,
   -3.4273,   13.7049,   -2.4446,    1.3659,    3.6702,   -0.2282,
   13.9254,   -2.4446,   20.2380,    3.6702,   -0.2282,   28.6779,
   13.7049,    1.3659,    3.6702,   12.5273,   -1.6045,    3.9419,
   -2.4446,    3.6702,   -0.2282,   -1.6045,    3.9419,    2.5821,
   20.2380,   -0.2282,   28.6779,    3.9419,    2.5821,   44.0636,
   };

   vnl_matrix<double> testV(testdataV, 6,6);*/

  imiRealType testdataV[9] =
  {
    0.8147, 0.9134, 0.2785,
    0.9058, 0.6324, 0.5469,
    0.1270, 0.0975, 0.9575
  };

  vnl_matrix<imiRealType> testV(testdataV, 3,3);

  /*imiRealType testdataZ[36] = {
   30.0000,   -3.4273,   1.9254,   23.7049,   -2.4446,   0.2380,
   -13.4273,   13.7049,   -20.4446,    1.3659,    3.6702,   10.2282,
   12.9254,   -22.4446,   0.2380,    8.702,   -0.2282,   8.6779,
   13.7049,    19.3659,    13.6702,   1.5273,   -1.6045,    3.9419,
   -20.4446,    3.6702,   -0.2282,   1.6045,    13.19,    2.5821,
   0.2380,   -1.2282,   8.6779,    17.9419,    0.5821,   4.0636,
   };*/

  /*double testdataZ[36] = {
   0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  2,
   0,  0,  0,  0, -1,  0,
   0,  0,  0,  2,  0,  0,
   };



   vnl_matrix<double> testZ(testdataZ, 6,6);*/

  imiRealType testdataZ[9] =
  { 1.8979, 1.7768, 1.0922,
    1.7768, 1.7922, 1.0734,
    1.0922, 1.0734, 1.0366
  };

  vnl_matrix<imiRealType> testZ(testdataZ, 3,3);

  m_regressorMatrix=testZ;
  m_regressandMatrix=testV;

  //-------------------------------------------
  // FIRST STEP:
  // Center regressor data
  //-------------------------------------------

  m_meanRegressorVector.set_size( m_regressorMatrix.rows(), 1 );
  m_meanRegressorVector.fill( 0.0 );

  if( m_regressorMatrix.empty() )
  {
    imiERROR( "CCA: Regressor matrix empty! Aborting computation." );
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
    imiERROR( "CCA: Observation matrix empty! Aborting computation." );
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
  m_centeredRegressandMatrix=m_regressandMatrix;
  for( unsigned int i = 0; i < m_regressandMatrix.rows(); i++ )
  {
    m_centeredRegressandMatrix.set_row( i, (m_regressandMatrix.get_row( i ) - m_meanRegressandVector( i, 0 )) );
  }

  //-------------------------------------------
  // THIRD STEP:
  // Compute CCA-estimator matrix (Due to memory limitations, only the parts of the estimator matrix B=V*Z^T*U_r*(SS^T)^-1*U_r^T are calculated.)
  //-------------------------------------------

  imiDEBUGINFO( 5, " Computing eigenvectors... ");

  if (m_mode==CCAModeA && m_deflationMode==CCADeflationModeRegression)
  {
    imiDEBUGINFO( 5, "(classical CCA, Mode A and deflation by regression)");
  }
  else if (m_mode==CCAModeB && m_deflationMode==CCADeflationModeCanonical)
  {
    imiDEBUGINFO( 5, "(CCA, Mode B and canonical deflation)");
  }

  imiDEBUGINFO( 5, " Starting NIPALS outer loop...");

  m_regressandComponents.set_size(m_centeredRegressandMatrix.rows(),m_numOfComponents);
  m_regressorComponents.set_size(m_centeredRegressorMatrix.rows(),m_numOfComponents);
  m_regressandEigenvalues.set_size(m_numOfComponents);
  m_regressorEigenvalues.set_size(m_numOfComponents);

  VnlMatrixType centeredRegressandMatrixResiduals=m_centeredRegressandMatrix;
  VnlMatrixType centeredRegressorMatrixResiduals=m_centeredRegressorMatrix;

  /*centeredRegressandMatrixResiduals=testV;
   centeredRegressorMatrixResiduals=testZ;
   m_regressandComponents.set_size(testV.rows(),m_numOfComponents);
   m_regressorComponents.set_size(testZ.rows(),m_numOfComponents);
   m_regressandScores.set_size(testV.cols(),m_numOfComponents);
   m_regressorScores.set_size(testZ.cols(),m_numOfComponents);*/

  /*imiDEBUGINFO( 5, "Computing eigensystem using svd: " << std::endl );
   VnlMatrixType matrix=testZ*testV.transpose();
   std::cout<<"matrix: "<<matrix<<std::endl;
   vnl_svd<imiRealType> svdSystem( matrix );
   for( unsigned int k = 0; k < 1; k++ )
   {
   vcl_cerr << "  " << k << "-th eigenvalue: " << svdSystem.W()[k] << " (" << svdSystem.U().get_column( k ) << ")" << std::endl;
   }*/

  VnlVectorType tempRegressorLoadings;
  for (unsigned int i=0;i<m_numOfComponents;i++)
  {

    if (m_mode==CCAModeB)
    {

      vnl_matrix<double> centeredRegressorMatrixResiduals_double( centeredRegressorMatrixResiduals.rows(), centeredRegressorMatrixResiduals.cols() );
      for( unsigned int i = 0; i < centeredRegressorMatrixResiduals.rows(); i++ )
      {
        for( unsigned int j = 0; j < centeredRegressorMatrixResiduals.cols(); j++ )
        {
          centeredRegressorMatrixResiduals_double[i][j] = static_cast<double>( centeredRegressorMatrixResiduals[i][j] );
        }
      }
      vnl_matrix<double> centeredRegressandMatrixResiduals_double( centeredRegressandMatrixResiduals.rows(), centeredRegressandMatrixResiduals.cols() );
      for( unsigned int i = 0; i < centeredRegressandMatrixResiduals.rows(); i++ )
      {
        for( unsigned int j = 0; j < centeredRegressandMatrixResiduals.cols(); j++ )
        {
          centeredRegressandMatrixResiduals_double[i][j] = static_cast<double>( centeredRegressandMatrixResiduals[i][j] );
        }
      }
      std::cout<<"X: "<<centeredRegressorMatrixResiduals<<std::endl;
      std::cout<<"Y: "<<centeredRegressandMatrixResiduals<<std::endl;
      vnl_matrix<double> regularizingMatrixX( centeredRegressorMatrixResiduals.rows(), centeredRegressorMatrixResiduals.cols() );
      regularizingMatrixX.fill( 0.0 );
      regularizingMatrixX.fill_diagonal( m_tolerance );
      vnl_matrix<double> regularizingMatrixY( centeredRegressandMatrixResiduals.rows(), centeredRegressandMatrixResiduals.cols() );
      regularizingMatrixY.fill( 0.0 );
      regularizingMatrixY.fill_diagonal( m_tolerance );
      vnl_qr<double> qr_v((centeredRegressandMatrixResiduals_double.transpose()*centeredRegressandMatrixResiduals_double)+regularizingMatrixY);
      vnl_matrix<double> Vpinv=qr_v.inverse()*centeredRegressandMatrixResiduals_double.transpose();

      vnl_qr<double> qr_z((centeredRegressorMatrixResiduals_double.transpose()*centeredRegressorMatrixResiduals_double)+regularizingMatrixX);
      vnl_matrix<double> Zpinv=qr_z.inverse()*centeredRegressorMatrixResiduals_double.transpose();

      std::cout<<"Zpinv: "<<Zpinv<<std::endl;
      std::cout<<"Vpinv: "<<Vpinv<<std::endl;

      imiDEBUGINFO( 5, "Computing eigensystem X using svd: " << std::endl );
      vnl_matrix<double> A=centeredRegressorMatrixResiduals_double*centeredRegressandMatrixResiduals_double.transpose()*(Vpinv.transpose()*Vpinv)*centeredRegressandMatrixResiduals_double*centeredRegressorMatrixResiduals_double.transpose();
      vnl_matrix<double> B=((centeredRegressorMatrixResiduals_double)*(centeredRegressorMatrixResiduals_double).transpose())+regularizingMatrixX;

      std::cout<<"A: --------------->"<<A<<std::endl;
      std::cout<<"B: --------------->"<<B<<std::endl;

      vnl_matrix<double> A_double( A.rows(), A.cols() );
      for( unsigned int i = 0; i < A.rows(); i++ )
      {
        for( unsigned int j = 0; j < A.cols(); j++ )
        {
          A_double[i][j] = static_cast<double>( A[i][j] );
        }
      }

      vnl_matrix<double> B_double( B.rows(), B.cols() );
      for( unsigned int i = 0; i < B.rows(); i++ )
      {
        for( unsigned int j = 0; j < B.cols(); j++ )
        {
          B_double[i][j] = static_cast<double>( B[i][j] );
        }
      }

      vnl_matrix<double> Zpinv_double( Zpinv.rows(), Zpinv.cols() );
      for( unsigned int i = 0; i < Zpinv.rows(); i++ )
      {
        for( unsigned int j = 0; j < Zpinv.cols(); j++ )
        {
          Zpinv_double[i][j] = static_cast<double>( Zpinv[i][j] );
        }
      }

      gsl_matrix *A_gsl=gsl_matrix_calloc(A_double.rows(),A_double.cols());
      for( unsigned int i = 0; i < A.rows(); i++ )
      {
        for( unsigned int j = 0; j < A.cols(); j++ )
        {
          gsl_matrix_set(A_gsl,i,j,A_double[i][j]);
        }
      }

      gsl_matrix *B_gsl=gsl_matrix_calloc(B_double.rows(),B_double.cols());
      for( unsigned int i = 0; i < B.rows(); i++ )
      {
        for( unsigned int j = 0; j < B.cols(); j++ )
        {
          gsl_matrix_set(B_gsl,i,j,B_double[i][j]);
        }
      }

      gsl_eigen_gensymmv_workspace *system = gsl_eigen_gensymmv_alloc (B_double.rows());
      gsl_vector *eval = gsl_vector_alloc (3);
      gsl_matrix *evec = gsl_matrix_alloc (3, 3);
      gsl_eigen_gensymmv(A_gsl, B_gsl,eval,evec, system);

      {
        int i;

        for (i = 0; i < 3; i++)
        {
          double eval_i
          = gsl_vector_get (eval, i);
          gsl_vector_view evec_i
          = gsl_matrix_column (evec, i);

          printf ("eigenvalue = %g\n", eval_i);
          printf ("eigenvector = \n");
          gsl_vector_fprintf (stdout,
              &evec_i.vector, "%g");
        }
      }

      vnl_generalized_eigensystem systemX(A_double,B_double);
      //vnl_real_eigensystem systemX((Zpinv_double.transpose()*Zpinv_double)*A_double);
      //std::cout<<"matrix X: "<<matrix<<std::endl;
      //vnl_svd<imiRealType> svdSystem( matrix );

      vcl_cout << "V = " << systemX.V << vcl_endl
      << "D = " << systemX.D << vcl_endl;

      return false;
      /*for( unsigned int k = 0; k < 6; k++ )
       {
       vcl_cerr << "  " << k << "-th eigenvalue: " << systemX.D[k] << " (" << systemX.V.get_column( k ) << ")" << std::endl;
       }*/
      //return true;
      /*imiDEBUGINFO( 5, "Computing eigensystem Y using svd: " << std::endl );
       matrix=centeredRegressandMatrixResiduals*centeredRegressorMatrixResiduals.transpose();
       std::cout<<"matrix Y: "<<matrix<<std::endl;
       vnl_svd<imiRealType> svdSystem2( matrix );
       for( unsigned int k = 0; k < 6; k++ )
       {
       vcl_cerr << "  " << k << "-th eigenvalue: " << svdSystem.W()[k] << " (" << svdSystem.V().get_column( k ) << ")" << std::endl;
       }*/

      /*imiDEBUGINFO( 5, "Computing eigensystem: " << std::endl );
       vnl_symmetric_eigensystem<MatrixValueType> eigSystem( matrix );
       for( unsigned int i = 0; i < matrix.rows(); i++ )
       {
       vcl_cerr << "  " << i << "-th eigenvalue: " << eigSystem.get_eigenvalue( i ) << " (" << eigSystem.get_eigenvector( i ) << ")" << std::endl;
       }*/

      /*imiDEBUGINFO( 5, "Computing eigensystem Y using svd: " << std::endl );
       matrix=centeredRegressandMatrixResiduals*centeredRegressorMatrixResiduals.transpose();
       std::cout<<"matrix Y: "<<matrix<<std::endl;
       vnl_svd<imiRealType> svdSystem2( matrix );
       for( unsigned int k = 0; k < 6; k++ )
       {
       vcl_cerr << "  " << k << "-th eigenvalue: " << svdSystem.W()[k] << " (" << svdSystem.V().get_column( k ) << ")" << std::endl;
       }*/

    }
    else
    {
      /*imiDEBUGINFO( 5, "Computing eigensystem X using svd: " << std::endl );
       VnlMatrixType matrix=centeredRegressorMatrixResiduals*centeredRegressandMatrixResiduals.transpose();
       std::cout<<"matrix X: "<<matrix<<std::endl;
       vnl_svd<imiRealType> svdSystem( matrix );
       for( unsigned int k = 0; k < 6; k++ )
       {
       vcl_cerr << "  " << k << "-th eigenvalue: " << svdSystem.W()[k] << " (" << svdSystem.U().get_column( k ) << ")" << std::endl;
       }

       imiDEBUGINFO( 5, "Computing eigensystem Y using svd: " << std::endl );
       matrix=centeredRegressandMatrixResiduals*centeredRegressorMatrixResiduals.transpose();
       std::cout<<"matrix Y: "<<matrix<<std::endl;
       vnl_svd<imiRealType> svdSystem2( matrix );
       for( unsigned int k = 0; k < 6; k++ )
       {
       vcl_cerr << "  " << k << "-th eigenvalue: " << svdSystem.W()[k] << " (" << svdSystem.V().get_column( k ) << ")" << std::endl;
       }*/

    }

    imiDEBUGINFO( 5, "\t Computing component "<<i+1<<"...");
    VnlVectorType tempRegressandWeights,tempRegressorWeights;
    VnlVectorType tempRegressandScores,tempRegressorScores;
    imiRealType tempEigenvalueRegressor,tempEigenvalueRegressand;

    NIPALSInnerLoop(centeredRegressorMatrixResiduals,centeredRegressandMatrixResiduals,tempRegressorWeights,tempRegressandWeights,tempRegressorScores,tempRegressandScores);

    tempEigenvalueRegressand=dot_product(tempRegressandScores,tempRegressorScores);
    tempEigenvalueRegressor=dot_product(tempRegressandScores,tempRegressorScores);

    imiDEBUGINFO( 5, "\t Done. Eigenvalues: ");
    //imiDEBUGINFO( 5, "\t\t Regressor: "<<tempEigenvalueRegressor);
    //imiDEBUGINFO( 5, "\t\t Regressand: "<<tempEigenvalueRegressand);
    imiDEBUGINFO( 5, "\t Done. Regressor: "<<tempRegressorWeights<<"; "<<tempEigenvalueRegressor<<std::endl);
    imiDEBUGINFO( 5, "\t Done. Regressand: "<<tempRegressandWeights<<"; "<<tempEigenvalueRegressand<<std::endl);

    m_regressandEigenvalues[i]=tempEigenvalueRegressand;
    m_regressorEigenvalues[i]=tempEigenvalueRegressor;

    m_regressandComponents.set_column(i,tempRegressandWeights/tempRegressandWeights.magnitude());
    m_regressorComponents.set_column(i,tempRegressorWeights);

    m_regressorScores.set_column(i,tempRegressorScores);
    m_regressandScores.set_column(i,tempRegressandScores);

    //Deflation part
    imiDEBUGINFO( 5, "\t Deflation...");

    //regress regressor matrix on regressor scores
    /*tempRegressorLoadings = (centeredRegressorMatrixResiduals* tempRegressorScores) / dot_product(tempRegressorScores,tempRegressorScores);


     //subtract rank-one approximation to obtain residual matrix
     centeredRegressorMatrixResiduals-=outer_product(tempRegressorLoadings,tempRegressorScores);

     VnlVectorType tempRegressandLoadings;
     if (m_deflationMode==CCADeflationModeCanonical)
     {
     // regress regressand matrix on regressAND scores
     tempRegressandLoadings=(centeredRegressandMatrixResiduals* tempRegressandScores) / dot_product(tempRegressandScores,tempRegressandScores);
     //subtract rank-one approximation to obtain residual matrix
     centeredRegressandMatrixResiduals-=outer_product(tempRegressandLoadings,tempRegressandScores);

     } else if (m_deflationMode==CCADeflationModeRegression)
     {
     // regress regressand matrix on regressOR scores
     tempRegressandLoadings=(centeredRegressandMatrixResiduals* tempRegressorScores) / dot_product(tempRegressorScores,tempRegressorScores);
     //subtract rank-one approximation to obtain residual matrix
     //centeredRegressandMatrixResiduals-=outer_product(tempRegressandLoadings,tempRegressorScores);
     centeredRegressandMatrixResiduals-=outer_product(tempRegressandWeights/tempRegressandWeights.magnitude(),(dot_product(tempRegressandScores,tempRegressorScores)/dot_product(tempRegressorScores,tempRegressorScores))*tempRegressorScores);
     //tempRegressandWeights/tempRegressandWeights.magnitude()
     } else
     {
     return false;
     }*/

    centeredRegressorMatrixResiduals-=((outer_product(m_regressorComponents.get_column(i),m_regressorComponents.get_column(i)))*centeredRegressorMatrixResiduals);
    //centeredRegressandMatrixResiduals-=((outer_product(m_regressandComponents.get_column(i),m_regressandComponents.get_column(i)))*centeredRegressandMatrixResiduals);

  }

  imiDEBUGINFO( 5, "\t Checking if regressor weights are pairwise orthogonal...");
  bool check=true;
  for(unsigned int i=0;i<m_regressorComponents.cols()-1;i++)
  {
    for(unsigned int j=i+1;j<m_regressorComponents.cols();j++)
    {
      if (dot_product(m_regressorComponents.get_column(i),m_regressorComponents.get_column(j))>m_tolerance)
      {
        check=false;
        std::cout<<"i: "<<i<<" j: "<<j<<"dot: "<<dot_product(m_regressorComponents.get_column(i),m_regressorComponents.get_column(j))<<std::endl;
      }
    }
  }
  if (check)
  {
    imiDEBUGINFO( 5, "\t Succesful.");
  }
  else
  {
    imiWARNING("Regressor weights are not pairwise orthogonal!")
  }

  /*imiDEBUGINFO( 5, "\t Checking if regressor weights are solutions of the eigenproblem...");
   for(unsigned int i=0;i<m_regressorComponents.cols()-1;i++)
   {
   //VnlMatrixType problem=testZ*testV.transpose()*testV*testZ.transpose();
   //VnlVectorType solution=problem*m_regressorComponents.get_column(i);
   vnl_qr<imiRealType> qr_v(m_centeredRegressandMatrix.transpose()*m_centeredRegressandMatrix);
   VnlMatrixType Vpinv=qr_v.inverse()*m_centeredRegressandMatrix.transpose();

   vnl_qr<imiRealType> qr_z(m_centeredRegressorMatrix.transpose()*m_centeredRegressorMatrix);
   VnlMatrixType Zpinv=qr_z.inverse()*m_centeredRegressorMatrix.transpose();

   imiDEBUGINFO( 5, "Computing eigensystem X using svd: " << std::endl );
   VnlMatrixType problem=m_centeredRegressorMatrix*m_centeredRegressandMatrix.transpose()*(Vpinv.transpose()*Vpinv)*m_centeredRegressandMatrix*m_centeredRegressorMatrix.transpose()*Zpinv.transpose();
   VnlVectorType solution=problem*m_regressorComponents.get_column(i);

   vnl_symmetric_eigensystem<MatrixValueType> eigSystem( problem );

   std::cout<<"---------------------------"<<std::endl;
   for (unsigned int j=0;j<solution.size();j++)
   {
   //std::cout<<solution[j]/m_regressorComponents.get_column(i)[j]<<std::endl;
   std::cout<<solution[j]/eigSystem.get_eigenvector(i)[j]<<std::endl;
   }

   std::cout<<"---------------------------"<<std::endl;



   }*/

  imiDEBUGINFO( 5, "\t Checking if regressand weights are pairwise orthogonal...");
  check=true;
  for(unsigned int i=0;i<m_regressandComponents.cols()-1;i++)
  {
    for(unsigned int j=i+1;j<m_regressandComponents.cols();j++)
    {
      //std::cout<<"i: "<<i<<" j: "<<j<<"dot: "<<dot_product(m_regressandComponents.get_column(i),m_regressandComponents.get_column(j))<<std::endl;
      if (dot_product(m_regressandComponents.get_column(i),m_regressandComponents.get_column(j))>m_tolerance)
      {
        check=false;
        std::cout<<"i: "<<i<<" j: "<<j<<"dot: "<<dot_product(m_regressandComponents.get_column(i),m_regressandComponents.get_column(j))<<std::endl;
      }
    }
  }
  if (check)
  {
    imiDEBUGINFO( 5, "\t Succesful.");
  }
  else
  {
    imiWARNING("Regressand weights are not pairwise orthogonal!")
  }

  // Preparation: Obviously, vnl rank works only for double matrices ...
  vnl_matrix<double> matrix_double( m_regressorComponents.rows(), m_regressorComponents.cols() );
  for( unsigned int i = 0; i < matrix_double.rows(); i++ )
  {
    for( unsigned int j = 0; j < matrix_double.cols(); j++ )
    {
      matrix_double[i][j] = static_cast<double>( m_regressorComponents[i][j] );
    }
  }

  // FIRST TEST: Compute rank of covariance matrix.
  // If rank is not full, the matrix is (in principle) not invertible.
  // Thus, a workaround is required.

  /* imiDEBUGINFO(5, "Dimension of matrix: " << m_regressorComponents.rows() );
   unsigned int covMatrixRank = vnl_rank( matrix_double );
   imiDEBUGINFO(5, "Rank of covariance matrix: " << covMatrixRank );*/

  check=false;

  return true;
}
#endif
