/** \file imiKernelPCA.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// Project includes
#include "imiKernelPCA.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

imiKernelPCA* imiKernelPCA::New()
{
  return new imiKernelPCA();
}

imiKernelPCA::imiKernelPCA()
{
  // Set default parameters:
  m_inputMatrix.set_size( 0, 0 );
  m_kernel = NULL;
  m_variabilityThreshold = 0.95;
  m_numOfComponentsUsed = -1;
  m_centered = true;
}

// imiKernelPCA::~imiKernelPCA()
imiKernelPCA::~imiKernelPCA()
{

}

// ----------------------------
//   Public methods
// ----------------------------

bool imiKernelPCA::ComputePrincipalComponents()
{
  imiDEBUGINFO( 5, "-------------------------------------------" );
  imiDEBUGINFO( 5, "  Kernel PCA (ComputePrincipalComponents)  " );
  imiDEBUGINFO( 5, "-------------------------------------------" );

  if( m_inputMatrix.empty() )
  {
    imiERROR( "KPCA: Data matrix empty! Aborting computation." );
    return false;
  }

  if( m_kernel == NULL )
  {
    imiERROR( "KPCA: No kernel specified! Aborting computation." );
    return false;
  }

  //-------------------------------------------
  // FIRST STEP:
  // Generate kernel matrix K and centering matrix H
  //-------------------------------------------

  //!!optimization potential: K is symmetric!! maybe...

  imiDEBUGINFO( 5, "  Computing kernel matrix ... " );

  m_K.set_size( m_inputMatrix.cols(), m_inputMatrix.cols() );

  for( unsigned int i = 0; i < m_K.rows(); i++ )
  {
    for( unsigned int j = 0; j < m_K.cols(); j++ )
    {
      m_K.put( i, j, m_kernel->Evaluate( m_inputMatrix.get_column( i ), m_inputMatrix.get_column( j ) ) );
    }
  }

  m_H.set_size( m_inputMatrix.cols(), m_inputMatrix.cols() );
  m_H.fill( -1.0 / m_inputMatrix.cols() );
  m_H.fill_diagonal( 1.0 - (1.0 / m_inputMatrix.cols()) );

  VnlMatrixType onesMat;
  onesMat.set_size( m_inputMatrix.cols(), m_inputMatrix.cols() );
  onesMat.fill( 1. / m_inputMatrix.cols() );

  if( m_centered )
  {
    m_centeredK = m_H*m_K*m_H;
//    m_centeredK = (m_K - (onesMat * m_K)) - (m_K * onesMat) + (onesMat * m_K * onesMat);
  }
  else
  {
    m_centeredK = m_K;
  }
  //VnlMatrixType m_centeredK = m_K-(onesMat*m_K)-(m_K*onesMat)+(onesMat*m_K*onesMat);

  //std::cout<<"kernelMatrix: "<<m_centeredK<<std::endl;

  //-------------------------------------------
  // SECOND STEP:
  // Eigenvalue decomposition of the centered kernel matrix
  //-------------------------------------------

  imiDEBUGINFO( 5, "  Eigenvalue decomposition of kernel matrix ... " );

  vnl_symmetric_eigensystem<MatrixValueType> eigSystem( m_centeredK );

  m_eigenvectors = eigSystem.V;
  m_eigenvectors.fliplr();

  //std::cout<<"eigenvectors: \n"<<m_eigenvectors<<std::endl;

  m_eigenvalues.set_size( eigSystem.D.size() );

  for( unsigned int i = 0, j = eigSystem.D.size() - 1; i < eigSystem.D.size(); i++, j-- )
  {
    m_eigenvalues[i] = eigSystem.D[j];
  }

  /*for (unsigned int i=0;i<m_eigenvectors.cols();i++)
   {
   m_eigenvectors.set_column(i,1./std::sqrt(m_eigenvalues[i])*m_eigenvectors.get_column(i));
   }*/

  imiDEBUGINFO( 5, "  Eigenvalues:" );
  for( unsigned int i = 0; i < m_eigenvalues.size(); i++ )
  {
    imiDEBUGINFO( 5, "    "<<i<<": "<<m_eigenvalues[i] );
  }

  return true;
}

VnlVectorType imiKernelPCA::ProjectDataPointOntoEigenvectors( VnlVectorType& data )
{
  /*imiDEBUGINFO( 5, "-------------------------------------------------" );
   imiDEBUGINFO( 5, "  Kernel PCA (ProjectDataPointOntoEigenvectors)  " );
   imiDEBUGINFO( 5, "-------------------------------------------------" );*/

  VnlVectorType projection;
  VnlMatrixType k_x;
  VnlMatrixType centeredk_x;
  VnlMatrixType ones;

  //std::cout<<"data: "<<data<<std::endl;

  ones.set_size( m_inputMatrix.cols(), 1 );
  ones.fill( -1. / m_inputMatrix.cols() );

  k_x.set_size( m_inputMatrix.cols(), 1 );

  for( unsigned int i = 0; i < k_x.size(); i++ )
  {
    k_x[i][0] = m_kernel->Evaluate( data, m_inputMatrix.get_column( i ) );
  }

  if( m_centered )
  {
    centeredk_x = m_H * (k_x + (m_K * ones));
  }
  else
  {
    centeredk_x = k_x;
  }

  unsigned int compsUsed = 0;

  if( m_numOfComponentsUsed < 0 )
  {
    compsUsed = NumOfCompsNeededForVarThresh();

  }
  else
  {
    compsUsed = m_numOfComponentsUsed;
  }


  projection.set_size( compsUsed );

  for( unsigned int i = 0; i < projection.size(); i++ )
  {

    //(1./std::sqrt(m_eigenvalues[i]))*
    projection[i] = (1. / std::sqrt( m_eigenvalues[i] )) * dot_product( m_eigenvectors.get_column( i ), centeredk_x.get_column( 0 ) );
  }

  //imiDEBUGINFO( 5, "  Projection: "<<projection<<" ("<<projection.magnitude()<<")" );

  /*VnlMatrixType onesVec;
   onesVec.set_size(m_inputMatrix.cols(),1);
   onesVec.fill(1./m_inputMatrix.cols());

   VnlMatrixType identMat;
   identMat.set_size(m_inputMatrix.cols(),m_inputMatrix.cols());
   identMat.fill(0.0);
   identMat.fill_diagonal(1.0);

   VnlMatrixType onesMat;
   onesMat.set_size(m_inputMatrix.cols(),m_inputMatrix.cols());
   onesMat.fill(1./m_inputMatrix.cols());

   k_x.set_size(m_inputMatrix.cols(),1);

   for (unsigned int i=0;i<k_x.size();i++)
   {
   k_x[i][0]=m_kernel->Evaluate(data,m_inputMatrix.get_column(i));
   }

   VnlMatrixType projMat = m_eigenvectors.transpose()*(identMat-onesMat)*(k_x-(m_K*onesVec));

   projection=projMat.get_column(0);

   std::cout<<"size: "<<projection.size()<<std::endl;*/

  return projection;
}

VnlVectorType imiKernelPCA::ComputePreImageOfProjection( VnlVectorType& projection )
{

  /*imiDEBUGINFO( 5, "-------------------------------------------------" );
   imiDEBUGINFO( 5, "  Kernel PCA (ComputePreImageOfProjection)       " );
   imiDEBUGINFO( 5, "-------------------------------------------------" );*/

  VnlVectorType preImage;
  preImage.set_size( m_inputMatrix.rows() );
  unsigned int compsUsed = 0;

  if( m_numOfComponentsUsed < 0 )
  {
    compsUsed = NumOfCompsNeededForVarThresh();

  }
  else
  {
    compsUsed = m_numOfComponentsUsed;
  }

  if( m_centered )
  {

    if( projection.size() < compsUsed )
    {
      compsUsed = projection.size();
    }
  }
  else
  {
    if( projection.size() - 1 < compsUsed )
    {
      compsUsed = projection.size() - 1;
    }
  }

  VnlMatrixType usedProjection;
  usedProjection.set_size( compsUsed, 1 );

  for( unsigned int i = 0; i < usedProjection.size(); i++ )
  {
    usedProjection[i][0] = projection[i]; //*(std::sqrt(m_eigenvalues[i]));
  }

  VnlMatrixType usedEigenvectors;
  usedEigenvectors.set_size( m_eigenvectors.rows(), compsUsed );

  VnlDiagMatrixType usedEigenvalues;
  usedEigenvalues.set_size( compsUsed );

  for( unsigned int i = 0; i < compsUsed; i++ )
  {
    //usedEigenvectors.set_column(i,m_eigenvectors.get_column(i)*(1/std::sqrt(m_eigenvalues[i])));
    usedEigenvalues[i] = m_eigenvalues[i];
    usedEigenvectors.set_column( i, m_eigenvectors.get_column( i ) );
  }

  imiDEBUGINFO( 5, "  Number of components used: "<<compsUsed<<"("<<projection.size()<<")" );

  /*
   //First step: compute the gamma vector (F.19, Canas)
   VnlMatrixType onesVec;
   onesVec.set_size(m_inputMatrix.cols(),1);
   onesVec.fill(1./m_inputMatrix.cols());

   VnlMatrixType identMat;
   identMat.set_size(m_inputMatrix.cols(),m_inputMatrix.cols());
   identMat.fill(0.0);
   identMat.fill_diagonal(1.0);

   VnlMatrixType onesMat;
   onesMat.set_size(m_inputMatrix.cols(),m_inputMatrix.cols());
   onesMat.fill(1./m_inputMatrix.cols());

   VnlMatrixType gammaVec;
   VnlMatrixType gammaVec2;
   VnlMatrixType gammaVecUncentered;

   gammaVecUncentered.set_size(m_inputMatrix.cols(),1);
   gammaVec2.set_size(m_inputMatrix.cols(),1);

   for (unsigned int i=0;i<gammaVecUncentered.rows();i++)
   {
   gammaVecUncentered[i][0]=.0;
   for (unsigned int j=0;j<usedEigenvectors.cols();j++)
   {
   gammaVecUncentered[i][0]+=usedProjection[j][0]*usedEigenvectors[j][i];
   }
   }

   MatrixValueType gammaSum=.0;

   for (unsigned int i=0;i<gammaVecUncentered.rows();i++)
   {
   gammaSum+=gammaVecUncentered[i][0];
   }

   for (unsigned int i=0;i<gammaVecUncentered.rows();i++)
   {
   //gammaSum+=gammaVecUncentered[i][0];
   gammaVec2[i][0]=gammaVecUncentered[i][0]+(1./m_inputMatrix.cols())*(1.-gammaSum);
   }




   gammaVec=onesVec+((identMat-onesMat)*(usedEigenvectors*usedProjection));

   imiDEBUGINFO( 5, "  Pre-image calculation..." );

   //Second step: pre-image calculation, closed-form solution based on Eq. 6 Rathi et al. (changed n to N)
   //only valid for gaussian kernel?

   VnlVectorType numerator;


   numerator.set_size(m_inputMatrix.rows());
   numerator.fill(0.);

   MatrixValueType denominator=0.;


   //compute feature space distance between x_i and the "projection" (F.22,Canas)
   VnlMatrixType tempMatrix=onesVec.transpose()*m_K*(2.0*gammaVec-onesVec);
   //const MatrixValueType featureDistanceConstPart=tempMatrix.get(0,0)+(usedProjection.get_column(0)).squared_magnitude();
   const MatrixValueType featureDistanceConstPart=tempMatrix.get(0,0)+(usedProjection.get_column(0)).squared_magnitude();


   for (unsigned int i=0; i<m_inputMatrix.cols();i++)
   {

   VnlMatrixType kernelEval;
   kernelEval.set_size(m_inputMatrix.cols(),1);

   for (unsigned int j=0; j<m_inputMatrix.cols();j++)
   {
   kernelEval[j][0]=m_kernel->Evaluate(m_inputMatrix.get_column(j),m_inputMatrix.get_column(i));
   }


   //std::cout<<"projSquared"<<(usedProjection.get_column(0)).squared_magnitude()<<std::endl;
   //std::cout<<"tempMatrix"<<tempMatrix.get(0,0)<<std::endl;


   tempMatrix=2.0*gammaVec.transpose()*kernelEval;
   MatrixValueType featureDistance=gammaVec[i][0]*(2.0-featureDistanceConstPart-tempMatrix[0][0]+m_kernel->Evaluate(m_inputMatrix.get_column(i),m_inputMatrix.get_column(i)));

   //std::cout<<"featureDist "<<featureDistance<<std::endl;

   numerator+=featureDistance*m_inputMatrix.get_column(i);
   denominator+=featureDistance;

   }



   preImage=numerator/denominator;*/

  /*std::cout<<std::endl;
   std::cout<<"numerator"<<numerator<<std::endl;
   std::cout<<"denominator"<<denominator<<std::endl;
   std::cout<<"preImage"<<preImage<<std::endl;*/

  preImage = m_kernel->ComputePreImageOfProjection( m_inputMatrix, usedProjection, usedEigenvectors, usedEigenvalues, m_centered );

  return preImage;
}

// ----------------------------
//   Protected / Private methods
// ----------------------------

unsigned int imiKernelPCA::NumOfCompsNeededForVarThresh() const
{
  MatrixValueType sum = 0.;
  unsigned int comps = 0;
  MatrixValueType parSum = 0.;

  unsigned int start;
  unsigned maxComps;

  if( m_centered )
  {
    start = 0;
    maxComps=m_eigenvalues.size()-1;
  }
  else
  {
    start = 1;
    maxComps=m_eigenvalues.size();
  }

  for( unsigned int i = start; i < maxComps; i++ )
  {
    sum += m_eigenvalues[i];
  }


  for( unsigned int i = start; comps == 0; i++ )
  {
    parSum += m_eigenvalues[i];
    if( parSum / sum > m_variabilityThreshold || i==maxComps )
    {
      comps = i;
    }
  }

  return comps;
}

} // namespace imi
