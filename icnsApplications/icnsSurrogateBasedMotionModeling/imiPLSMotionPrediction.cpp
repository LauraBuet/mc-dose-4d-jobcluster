/** \file imiPLSMotionPrediction.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 *  This code is partially based on (Python) code from
 *
 *  http://scikit-learn.org,
 *
 *  especially https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/pls.py,
 *
 *  and implements Partial Least Squares (PLS).
 *
 *  @article{scikit-learn,
 *  title={{Scikit-learn: Machine Learning in Python }},
 *  author={Pedregosa, F. and Varoquaux, G. and Gramfort, A. and Michel, V.
 *         and Thirion, B. and Grisel, O. and Blondel, M. and Prettenhofer, P.
 *        and Weiss, R. and Dubourg, V. and Vanderplas, J. and Passos, A. and
 *        Cournapeau, D. and Brucher, M. and Perrot, M. and Duchesnay, E.},
 *  journal={Journal of Machine Learning Research},
 *  volume={12},
 *  pages={2825--2830},
 *  year={2011}
 *  }
 *
 ****************************************************************************/
// Project includes
#include "imiPLSMotionPrediction.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

imiPLSMotionPrediction* imiPLSMotionPrediction::New()
{
  return new imiPLSMotionPrediction();
}

imiPLSMotionPrediction::imiPLSMotionPrediction()
{
  // Set default parameters:
  m_mode = PLSModeA;
  m_deflationMode = PLSDeflationModeCanonical;
  m_numOfComponents = 2;
  m_useBothBases = false;

  m_maxIters = 500;
  m_tolerance = 1e-6;

  m_TikhonovRegularizationParameter = 0.001;

}

// imiPLSMotionPrediction::~imiPLSMotionPrediction()
imiPLSMotionPrediction::~imiPLSMotionPrediction()
{
}

// ----------------------------
//   Public methods
// ----------------------------

// imiPLSMotionPrediction::TrainPLSEstimator
bool imiPLSMotionPrediction::TrainEstimator( void )
{
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "  TRAINING OF PLS ESTIMATOR" );
  imiDEBUGINFO( 5, "-----------------------------------" );

  //-------------------------------------------
  // FIRST STEP:
  // Center regressor data
  //-------------------------------------------

  m_meanRegressorVector.set_size( m_regressorMatrix.rows(), 1 );
  m_meanRegressorVector.fill( 0.0 );

  if( m_regressorMatrix.empty() )
  {
    imiERROR( "PLS: Regressor matrix empty! Aborting computation." );
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
    imiERROR( "PLS: Observation matrix empty! Aborting computation." );
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
  // THIRD STEP:
  // Compute PLS-estimator matrix (Due to memory limitations, only the parts of the estimator matrix B=V*Z^T*U_r*(SS^T)^-1*U_r^T are calculated.)
  //-------------------------------------------

  imiDEBUGINFO( 5, " Computing eigenvectors... " );

  if( m_mode == PLSModeA && m_deflationMode == PLSDeflationModeRegression )
  {
    imiDEBUGINFO( 5, "(classical PLS, Mode A and deflation by regression)" );
  }
  else if( m_mode == PLSModeB && m_deflationMode == PLSDeflationModeCanonical )
  {
    imiDEBUGINFO( 5, "(CCA, Mode B and canonical deflation)" );
  }

  imiDEBUGINFO( 5, " Starting NIPALS outer loop..." );

  m_regressandComponents.set_size( m_centeredRegressandMatrix.rows(), m_numOfComponents );
  m_regressorComponents.set_size( m_centeredRegressorMatrix.rows(), m_numOfComponents );
  m_regressandEigenvalues.set_size( m_numOfComponents );
  m_regressorEigenvalues.set_size( m_numOfComponents );
  m_regressandScores.set_size( m_centeredRegressandMatrix.cols(),m_numOfComponents);
  m_regressorScores.set_size( m_centeredRegressorMatrix.cols(), m_numOfComponents );

  VnlMatrixType centeredRegressandMatrixResiduals = m_centeredRegressandMatrix;
  VnlMatrixType centeredRegressorMatrixResiduals = m_centeredRegressorMatrix;

  VnlVectorType tempRegressorLoadings;
  for( unsigned int i = 0; i < m_numOfComponents; i++ )
  {

    imiDEBUGINFO( 5, "\t Computing component "<<i+1<<"..." );
    VnlVectorType tempRegressandWeights, tempRegressorWeights;
    VnlVectorType tempRegressandScores, tempRegressorScores;
    imiRealType tempEigenvalueRegressor, tempEigenvalueRegressand;

    NIPALSInnerLoop( centeredRegressorMatrixResiduals, centeredRegressandMatrixResiduals, tempRegressorWeights, tempRegressandWeights, tempRegressorScores, tempRegressandScores );

    tempEigenvalueRegressand = dot_product( tempRegressandScores, tempRegressorScores );
    tempEigenvalueRegressor = dot_product( tempRegressandScores, tempRegressorScores );

    imiDEBUGINFO( 5, "\t Done. Eigenvalues: " );
    imiDEBUGINFO( 5, "\t\t Regressor: "<<tempEigenvalueRegressor );
    imiDEBUGINFO( 5, "\t\t Regressand: "<<tempEigenvalueRegressand );
    //imiDEBUGINFO( 5, "\t Done. Regressor: "<<tempRegressorWeights<<"; "<<tempEigenvalueRegressor<<std::endl);
    //imiDEBUGINFO( 5, "\t Done. Regressand: "<<tempRegressandWeights<<"; "<<tempEigenvalueRegressand<<std::endl);
    imiDEBUGINFO( 5, "\t Done. Regressor: "<<tempRegressorScores<<std::endl);

    m_regressandEigenvalues[i] = tempEigenvalueRegressand;
    m_regressorEigenvalues[i] = tempEigenvalueRegressor;

    m_regressandComponents.set_column( i, tempRegressandWeights / tempRegressandWeights.magnitude() );
    m_regressorComponents.set_column( i, tempRegressorWeights );

    m_regressorScores.set_column( i, tempRegressorScores );
    m_regressandScores.set_column( i, tempRegressandScores );

    imiDEBUGINFO( 5, "\t Done. Regressor: "<<m_regressorScores.get_column( i )<<std::endl);

    //Deflation part
    imiDEBUGINFO( 5, "\t Deflation..." );

    //regress regressor matrix on regressor scores
    /*tempRegressorLoadings = (centeredRegressorMatrixResiduals* tempRegressorScores) / dot_product(tempRegressorScores,tempRegressorScores);


     //subtract rank-one approximation to obtain residual matrix
     centeredRegressorMatrixResiduals-=outer_product(tempRegressorLoadings,tempRegressorScores);

     VnlVectorType tempRegressandLoadings;
     if (m_deflationMode==PLSDeflationModeCanonical)
     {
     // regress regressand matrix on regressAND scores
     tempRegressandLoadings=(centeredRegressandMatrixResiduals* tempRegressandScores) / dot_product(tempRegressandScores,tempRegressandScores);
     //subtract rank-one approximation to obtain residual matrix
     centeredRegressandMatrixResiduals-=outer_product(tempRegressandLoadings,tempRegressandScores);

     } else if (m_deflationMode==PLSDeflationModeRegression)
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

   VnlMatrixType regressorComponentMatrix( m_regressorComponents.get_column( i ).size(), 1 );
    regressorComponentMatrix.set_column( 0, m_regressorComponents.get_column( i ) );

    VnlMatrixType regressandComponentMatrix( m_regressandComponents.get_column( i ).size(), 1 );
    regressandComponentMatrix.set_column( 0, m_regressandComponents.get_column( i ) );


    centeredRegressorMatrixResiduals -= regressorComponentMatrix * (regressorComponentMatrix.transpose() * centeredRegressorMatrixResiduals);
    centeredRegressandMatrixResiduals -= regressandComponentMatrix * (regressandComponentMatrix.transpose() * centeredRegressandMatrixResiduals);
    //centeredRegressorMatrixResiduals-=((outer_product(m_regressorComponents.get_column(i),m_regressorComponents.get_column(i)))*centeredRegressorMatrixResiduals);
    //centeredRegressandMatrixResiduals-=((outer_product(m_regressandComponents.get_column(i),m_regressandComponents.get_column(i)))*centeredRegressandMatrixResiduals);

  }

  imiDEBUGINFO( 5, "\t Checking if regressor weights are pairwise orthogonal..." );
  bool check = true;
  for( unsigned int i = 0; i < m_regressorComponents.cols() - 1; i++ )
  {
    for( unsigned int j = i + 1; j < m_regressorComponents.cols(); j++ )
    {
      if( dot_product( m_regressorComponents.get_column( i ), m_regressorComponents.get_column( j ) ) > m_tolerance )
      {
        check = false;
        std::cout << "i: " << i << " j: " << j << "dot: " << dot_product( m_regressorComponents.get_column( i ), m_regressorComponents.get_column( j ) ) << std::endl;
      }
    }
  }
  if( check )
  {
    imiDEBUGINFO( 5, "\t Succesful." );
  }
  else
  {
    imiWARNING( "Regressor weights are not pairwise orthogonal!" )
  }

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
        std::cout << "i: " << i << " j: " << j << "dot: " << dot_product( m_regressandComponents.get_column( i ), m_regressandComponents.get_column( j ) ) << std::endl;
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

  return true;
}

// imiPLSMotionPrediction::PredictOutput
bool imiPLSMotionPrediction::PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput )
{
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "           PREDICTION              " );
  imiDEBUGINFO( 5, "-----------------------------------" );
  imiDEBUGINFO( 5, "  Predicting output for given regressor measurement ... " );

  if( m_meanRegressandVector.empty() || m_meanRegressorVector.empty() || m_centeredRegressandMatrix.empty() || m_centeredRegressorMatrix.empty() )
  {
    imiERROR( "Prediction not possible: PLS estimator not trained. Call TrainEstimator() before trying to predict." );
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
    predictedOutput = m_meanRegressandVector
        + (predRegressandComponents
            * (predRegressandComponents.transpose() * m_centeredRegressandMatrix
                * (m_centeredRegressorMatrix.transpose() * (predRegressorComponents * (scoresInv * (predRegressorComponents.transpose() * (regressorMeasurement - m_meanRegressorVector)))))));
  }
  else
  {
    predictedOutput = m_meanRegressandVector
        + (m_centeredRegressandMatrix
            * (m_centeredRegressorMatrix.transpose() * (predRegressorComponents * (scoresInv * (predRegressorComponents.transpose() * (regressorMeasurement - m_meanRegressorVector))))));
  }

  imiDEBUGINFO( 5, "  finished." );
  return true;
}

/*unsigned int imiPLSMotionPrediction::PerformLOOCV(unsigned int minNumberOfComponents, unsigned int maxNumberOfComponents)
{
  int oldDebugLevel=imiObject::GetGlobalDebugLevel();
  imiObject::SetGlobalDebugLevel( 0 );



  for (unsigned int i)


  imiObject::SetGlobalDebugLevel( oldDebugLevel );
return 0;
}*/

// ----------------------------
//   Protected / Private methods
// ----------------------------

/***
 * Inner loop of the iterative NIPALS algorithm. Provides an iterative way to compute eigenvectors
 * of the generalized eigenproblem associated with PLS and CCA.
 *
 * For details see:
 *
 * Overview and recent advances in partial least squares
 * R Rosipal, N Krämer
 * Subspace, Latent Structure and Feature Selection, 34-51, 2006
 *
 * A survey of Partial Least Squares (PLS) methods, with emphasis on the two-block case.
 * Jacob A. Wegelin
 * Technical Report 371, Department of Statistics, University of Washington, Seattle, 2000
 *
 *
 * @param regressorWeights
 * @param regressandWeights
 * @return
 */
bool imiPLSMotionPrediction::NIPALSInnerLoop( VnlMatrixType &regressorMatrix, VnlMatrixType &regressandMatrix, VnlVectorType &regressorWeights, VnlVectorType &regressandWeights,
    VnlVectorType &regressorScores, VnlVectorType &regressandScores )
{

  //Z: regressor, V: regressand

  VnlVectorType vWeights, zWeights;

  VnlVectorType vScores, zScores;
  //random initialisation
  vScores.set_size( regressandMatrix.cols() );
  vScores.fill( 0.0 );
  vScores = vScores + (regressandMatrix.get_row( 0 ));

  VnlMatrixType Vpinv, Zpinv;

  //For mode b (CCA), pseudo inverse matrices Z^+ and V^+ have to be computed
  if( m_mode == PLSModeB )
  {
    vnl_qr<imiRealType> qr_v( regressandMatrix.transpose() * regressandMatrix );
    Vpinv = qr_v.inverse() * regressandMatrix.transpose();

    vnl_qr<imiRealType> qr_z( regressorMatrix.transpose() * regressorMatrix );
    Zpinv = qr_z.inverse() * regressorMatrix.transpose();
  }

  unsigned int iters = 0;
  imiRealType diff = 1.0;
  VnlVectorType zWeightsOld( regressorMatrix.rows() );
  VnlVectorType zScoresOld( regressorMatrix.cols() );

  while( iters < m_maxIters && diff > m_tolerance )
  {
    //update regressor weights
    if( m_mode == PLSModeA )
    {
      zWeights = (regressorMatrix * vScores) / dot_product( vScores, vScores );
    }
    else if( m_mode == PLSModeB )
    {
      zWeights = (Zpinv.transpose() * vScores);
    }
    else
    {
      return false;
    }

    //normalize z weights
    zWeights.normalize();

    //Update the Z latent scores
    zScores = regressorMatrix.transpose() * zWeights;

    //Update regressand weights
    if( m_mode == PLSModeA )
    {
      vWeights = (regressandMatrix * zScores) / dot_product( zScores, zScores );
    }
    else if( m_mode == PLSModeB )
    {
      vWeights = (Vpinv.transpose() * zScores);
    }
    else
    {
      return false;
    }

    if( m_deflationMode == PLSDeflationModeCanonical )
    {
      vWeights.normalize();
    }

    //update the V latent scores
    vScores = (regressandMatrix.transpose() * vWeights) / dot_product( vWeights, vWeights );

    diff = dot_product( zWeights - zWeightsOld, zWeights - zWeightsOld );
    //diff=dot_product(zScores-zScoresOld,zScores-zScoresOld);
    iters++;
    zWeightsOld = zWeights;
    //zScoresOld=zScores;
    imiDEBUGINFO( 9, "\t\tIteration: "<<iters<<" diff: "<<diff );
  }

  regressorWeights = zWeights;
  regressandWeights = vWeights;
  regressorScores = zScores;
  regressandScores = vScores;

  return true;
}

} // namespace imi

/*
 * bool imiPLSMotionPrediction::NIPALSInnerLoop(VnlMatrixType &regressorMatrix, VnlMatrixType &regressandMatrix, VnlVectorType &regressorWeights,VnlVectorType &regressandWeights,VnlVectorType &regressorScores, VnlVectorType &regressandScores)
 {


 //Z: regressor, V: regressand

 VnlVectorType vWeights,zWeights;

 VnlVectorType vScores,zScores;
 //random initialisation
 vScores.set_size(regressandMatrix.cols());
 vScores.fill(0.0);
 vScores=vScores+(regressandMatrix.get_row(0));

 VnlMatrixType Vpinv,Zpinv;

 //For mode b (CCA), pseudo inverse matrices Z^+ and V^+ have to be computed
 if (m_mode==PLSModeB)
 {
 vnl_qr<imiRealType> qr_v(regressandMatrix.transpose()*regressandMatrix);
 Vpinv=qr_v.inverse()*regressandMatrix.transpose();

 vnl_qr<imiRealType> qr_z(regressorMatrix.transpose()*regressorMatrix);
 Zpinv=qr_z.inverse()*regressorMatrix.transpose();
 }

 unsigned int iters=0;
 imiRealType diff=1.0;
 VnlVectorType zWeightsOld(regressorMatrix.rows());
 VnlVectorType zScoresOld(regressorMatrix.cols());


 while (iters<m_maxIters && diff > m_tolerance )
 {
 //update regressor weights
 if (m_mode==PLSModeA)
 {
 zWeights=(regressorMatrix*vScores)/dot_product(vScores,vScores);
 } else if (m_mode==PLSModeB) {
 zWeights=(Zpinv*vScores);
 } else {
 return false;
 }

 //normalize z weights
 zWeights.normalize();


 //Update the Z latent scores
 zScores=regressorMatrix.transpose()*zWeights;


 //Update regressand weights
 if (m_mode==PLSModeA)
 {
 vWeights=(regressandMatrix*zScores)/dot_product(zScores,zScores);
 } else if (m_mode==PLSModeB) {
 vWeights=(Vpinv.transpose()*zScores);
 } else {
 return false;
 }

 if (m_deflationMode==PLSDeflationModeCanonical)
 {
 vWeights.normalize();
 }

 //update the V latent scores
 vScores=(regressandMatrix.transpose()*vWeights)/dot_product(vWeights,vWeights);

 //diff=dot_product(zWeights-zWeightsOld,zWeights-zWeightsOld);
 diff=dot_product(zScores-zScoresOld,zScores-zScoresOld);
 iters++;
 //zWeightsOld=zWeights;
 zScoresOld=zScores;
 imiDEBUGINFO( 9, "\t\tIteration: "<<iters<<" diff: "<<diff);
 }

 regressorWeights=zWeights;
 regressandWeights=vWeights;
 regressorScores=zScores;
 regressandScores=vScores;

 return true;
 }*/

/*NIPALS deBie
 * bool imiPLSMotionPrediction::NIPALSInnerLoop(VnlMatrixType &regressorMatrix, VnlMatrixType &regressandMatrix, VnlVectorType &regressorWeights,VnlVectorType &regressandWeights,VnlVectorType &regressorScores, VnlVectorType &regressandScores)
 {


 //Z: regressor, V: regressand

 VnlVectorType vWeights,zWeights;

 VnlVectorType vScores,zScores;
 //random initialisation
 vScores.set_size(regressandMatrix.cols());
 vScores.fill(0.0);
 vScores=vScores+(regressandMatrix.get_row(0));

 vWeights.set_size(regressandMatrix.rows());
 vWeights.fill(1.0);

 VnlMatrixType Vpinv,Zpinv;

 //For mode b (CCA), pseudo inverse matrices Z^+ and V^+ have to be computed
 if (m_mode==PLSModeB)
 {
 vnl_qr<imiRealType> qr_v(regressandMatrix.transpose()*regressandMatrix);
 Vpinv=qr_v.inverse()*regressandMatrix.transpose();

 vnl_qr<imiRealType> qr_z(regressorMatrix.transpose()*regressorMatrix);
 Zpinv=qr_z.inverse()*regressorMatrix.transpose();
 }

 unsigned int iters=0;
 imiRealType diff=1.0;
 VnlVectorType zWeightsOld(regressorMatrix.rows());
 VnlVectorType zScoresOld(regressorMatrix.cols());



 while (iters<m_maxIters && diff > m_tolerance )
 {

 zScores=regressandMatrix.transpose()*vWeights;


 //update regressor weights
 if (m_mode==PLSModeA)
 {
 zWeights=(regressorMatrix*zScores);
 } else if (m_mode==PLSModeB) {
 zWeights=(Zpinv*vScores);
 } else {
 return false;
 }

 std::cout<<"Eigenvalue x: "<<zWeights.magnitude()<<std::endl;

 //normalize z weights
 zWeights.normalize();


 vScores=regressorMatrix.transpose()*zWeights;


 //Update regressand weights
 if (m_mode==PLSModeA)
 {
 vWeights=(regressandMatrix*vScores);
 std::cout<<"Eigenvalue y: "<<vWeights.magnitude()<<std::endl;
 vWeights.normalize();
 } else if (m_mode==PLSModeB) {
 vWeights=(Vpinv.transpose()*zScores);
 } else {
 return false;
 }

 if (m_deflationMode==PLSDeflationModeCanonical)
 {
 vWeights.normalize();
 }

 //update the V latent scores
 vScores=(regressandMatrix.transpose()*vWeights)/dot_product(vWeights,vWeights);

 //diff=dot_product(zWeights-zWeightsOld,zWeights-zWeightsOld);
 diff=dot_product(zScores-zScoresOld,zScores-zScoresOld);
 iters++;
 //zWeightsOld=zWeights;
 zScoresOld=zScores;
 imiDEBUGINFO( 9, "\t\tIteration: "<<iters<<" diff: "<<diff);
 }

 regressorWeights=zWeights;
 regressandWeights=vWeights;
 regressorScores=zScores;
 regressandScores=vScores;

 return true;
 }
 */
