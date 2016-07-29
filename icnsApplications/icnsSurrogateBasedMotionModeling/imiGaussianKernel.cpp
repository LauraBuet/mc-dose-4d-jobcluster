/** \file imiGaussianKernel.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#include <cmath>


// Project includes
#include "imiGaussianKernel.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

imiGaussianKernel* imiGaussianKernel::New()
{
  return new imiGaussianKernel();
}

imiGaussianKernel::imiGaussianKernel()
{
  // Set default parameters:
  m_sigma=.5;
}

// imiGaussianKernel::~imiGaussianKernel()
imiGaussianKernel::~imiGaussianKernel()
{
}

// ----------------------------
//   Public methods
// ----------------------------

MatrixValueType imiGaussianKernel::Evaluate(const VnlVectorType& v1, const VnlVectorType& v2) const
{
  return std::exp(-vnl_vector_ssd(v1,v2)/(2.*m_sigma*m_sigma));
}

MatrixValueType imiGaussianKernel::Evaluate(MatrixValueType dist) const
{
  return std::exp(-(dist*dist)/(2.*m_sigma*m_sigma));///(m_sigma*vnl_math_sqr(2*vnl_math::pi));
}

MatrixValueType imiGaussianKernel::ComputeSigmaValueArias(VnlMatrixType& inputVectors) const
{
  MatrixValueType sigma;

  MatrixValueType meanNNDist=0.;

  //Canas F.10
  for(unsigned int i=0; i<inputVectors.cols();i++)
  {
    MatrixValueType NNDist=std::numeric_limits<MatrixValueType>::max();
    for (unsigned int j=0;j<inputVectors.cols();j++)
    {
      if(i!=j)
      {
        MatrixValueType tempDist=std::sqrt(vnl_vector_ssd(inputVectors.get_column(i),inputVectors.get_column(j)));
        if(tempDist<NNDist)
        {
          NNDist=tempDist;
        }
      }
    }
    meanNNDist+=1./inputVectors.cols()*NNDist;
  }

  sigma=meanNNDist;

  return sigma;
}

MatrixValueType imiGaussianKernel::ComputeSigmaValueCasero(VnlMatrixType& inputVectors) const
{
  MatrixValueType sigma;

  MatrixValueType meanNNDist=0.;

  //Canas F.10
  for(unsigned int i=0; i<inputVectors.cols();i++)
  {
    MatrixValueType NNDist=std::numeric_limits<MatrixValueType>::max();
    for (unsigned int j=0;j<inputVectors.cols();j++)
    {
      if(i!=j)
      {
        MatrixValueType tempDist=std::sqrt(vnl_vector_ssd(inputVectors.get_column(i),inputVectors.get_column(j)));
        if(tempDist<NNDist)
        {
          NNDist=tempDist;
        }
      }
    }
    meanNNDist+=1./inputVectors.cols()*NNDist;
  }

  sigma=1./std::sqrt(-2.*std::log(0.5))*meanNNDist;

  return sigma;
}

MatrixValueType imiGaussianKernel::ComputeSigmaValueCremers1(VnlMatrixType& inputVectors) const
{
  MatrixValueType sigma;

  MatrixValueType meanNNDist=0.;

  //Canas F.9b
  for(unsigned int i=0; i<inputVectors.cols();i++)
  {
    MatrixValueType NNDist=std::numeric_limits<MatrixValueType>::max();
    for (unsigned int j=0;j<inputVectors.cols();j++)
    {
      if(i!=j)
      {
        MatrixValueType tempDist=(vnl_vector_ssd(inputVectors.get_column(i),inputVectors.get_column(j)));
        if(tempDist<NNDist)
        {
          NNDist=tempDist;
        }
      }
    }
    meanNNDist+=1./inputVectors.cols()*NNDist;
  }

  sigma=std::sqrt(meanNNDist);

  return sigma;
}

MatrixValueType imiGaussianKernel::ComputeSigmaValueCremers2(VnlMatrixType& inputVectors) const
{
  return 1.5*ComputeSigmaValueCremers1(inputVectors);
}



VnlVectorType imiGaussianKernel::ComputePreImageOfProjection(VnlMatrixType& trainingData, VnlMatrixType& usedProjection,VnlMatrixType& usedEigenvectors,VnlDiagMatrixType& usedEigenvalues, bool centered) const
{
  //not centered version!!!!!!!!!!!!!!!!!!!!!!!!

    VnlVectorType numerator;


    numerator.set_size(trainingData.rows());
    numerator.fill(0.);

    MatrixValueType denominator=0.;

    /*VnlMatrixType onesVec;
    onesVec.set_size(trainingData.cols(),1);
    onesVec.fill(1./trainingData.cols());

    VnlMatrixType identMat;
    identMat.set_size(trainingData.cols(),trainingData.cols());
    identMat.fill(0.0);
    identMat.fill_diagonal(1.0);

    VnlMatrixType onesMat;
    onesMat.set_size(trainingData.cols(),trainingData.cols());
    onesMat.fill(1./trainingData.cols());

    VnlMatrixType gamma = onesVec+(identMat-onesMat)*usedEigenvectors*usedProjection;*/


    for (unsigned int i=0; i<trainingData.cols();i++)
    {
      MatrixValueType gamma=.0;

      MatrixValueType k_i=.0;

      for (unsigned int k=0; k<usedEigenvectors.cols();k++)
      {
        gamma+=usedEigenvectors.get_column(k)[i]/std::sqrt(usedEigenvalues[k])*usedProjection[k][0];//+(1./m_inputMatrix.cols()*(1.-gammaSum));
        k_i+=usedEigenvectors.get_column(k)[i]*std::sqrt(usedEigenvalues[k])*usedProjection[k][0];
        //gamma+=usedEigenvectors.get_column(k)[i]*usedProjection[k][0];
        //k_i+=usedEigenvectors.get_column(k)[i]*usedProjection[k][0];
      }

      numerator+=gamma*k_i*trainingData.get_column(i);
      denominator+=gamma*k_i;


    }

    VnlVectorType preImage=numerator/denominator;


    return preImage;
}

// ----------------------------
//   Protected / Private methods
// ----------------------------

}// namespace imi
