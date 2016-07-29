/** \file imiMotionPredictionLOOCV.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/
#include "imiMotionPredictionLOOCV.h"

#include "imiPCRMotionPrediction.h"
#include "imiPLSMotionPrediction.h"
#include "imiCCAMotionPrediction.h"
#include "imiKernelMLRMotionPrediction.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

imiMotionPredictionLOOCV* imiMotionPredictionLOOCV::New()
{
  return new imiMotionPredictionLOOCV();
}

imiMotionPredictionLOOCV::imiMotionPredictionLOOCV()
{
  // Set default parameters:
}

// imiPLSMotionPrediction::~imiPLSMotionPrediction()
imiMotionPredictionLOOCV::~imiMotionPredictionLOOCV()
{
}

// ----------------------------
//   Public methods
// ----------------------------

unsigned int imiMotionPredictionLOOCV::PerformLOOCVPCR(unsigned int minNumOfComps, unsigned int maxNumOfComps)
{
  double minValue=std::numeric_limits<double>::max();
  unsigned int minValueComps=0;

  imiPCRMotionPrediction* pcrPredictor=imiPCRMotionPrediction::New();
  m_predictor=pcrPredictor;

  std::vector<double> squaredDists;


  for (unsigned int i=minNumOfComps; i<=maxNumOfComps;i++)
  { 
    squaredDists.push_back(.0);
  }

  this->PerformLOOCV(minNumOfComps,maxNumOfComps,squaredDists);

  for (unsigned int index=0; index<squaredDists.size();index++)
  {
    double diff=squaredDists[index];
    if(diff<minValue)
    {
      minValue=diff;
      minValueComps=minNumOfComps+index;
    }
    imiINFO("LOOCV components: "<<minNumOfComps+index<<" diff: "<<diff);
  }

  return minValueComps;
}

unsigned int imiMotionPredictionLOOCV::PerformLOOCVPLS(unsigned int minNumOfComps, unsigned int maxNumOfComps)
{
  double minValue=std::numeric_limits<double>::max();
  unsigned int minValueComps=0;

  imiPLSMotionPrediction* plsPredictor=imiPLSMotionPrediction::New();
  plsPredictor->SetUseBothBases(false);
  m_predictor=plsPredictor;

  std::vector<double> squaredDists;


  for (unsigned int i=minNumOfComps; i<=maxNumOfComps;i++)
  { 
    squaredDists.push_back(.0);
  }

  this->PerformLOOCV(minNumOfComps,maxNumOfComps,squaredDists);

  for (unsigned int index=0; index<squaredDists.size();index++)
  {
    double diff=squaredDists[index];
    if(diff<minValue)
    {
      minValue=diff;
      minValueComps=minNumOfComps+index;
    }
    imiINFO("LOOCV components: "<<minNumOfComps+index<<" diff: "<<diff);
  }

  return minValueComps;
}

unsigned int imiMotionPredictionLOOCV::PerformLOOCVCCA(unsigned int minNumOfComps, unsigned int maxNumOfComps)
{
  double minValue=std::numeric_limits<double>::max();
  unsigned int minValueComps=0;

  imiCCAMotionPrediction* ccaPredictor=imiCCAMotionPrediction::New();
  ccaPredictor->SetUseBothBases(false);
  m_predictor=ccaPredictor;

  std::vector<double> squaredDists;


  for (unsigned int i=minNumOfComps; i<=maxNumOfComps;i++)
  { 
    squaredDists.push_back(.0);
  }

  this->PerformLOOCV(minNumOfComps,maxNumOfComps,squaredDists);

  for (unsigned int index=0; index<squaredDists.size();index++)
  {
    double diff=squaredDists[index];
    if(diff<minValue)
    {
      minValue=diff;
      minValueComps=minNumOfComps+index;
    }
    imiINFO("LOOCV components: "<<minNumOfComps+index<<" diff: "<<diff);
  }

  return minValueComps;
}

void imiMotionPredictionLOOCV::PerformLOOCV(unsigned int minNumOfComps, unsigned int maxNumOfComps, std::vector<double> &squaredDists)
{
  VnlMatrixType newRegressand;
  VnlMatrixType newRegressor;
  VnlMatrixType measurement;
  VnlMatrixType reference;



  newRegressand.set_size(m_regressandMatrix.rows(),m_regressandMatrix.cols()-1);
  newRegressor.set_size(m_regressorMatrix.rows(),m_regressorMatrix.cols()-1);
  measurement.set_size(m_regressorMatrix.rows(),1);
  reference.set_size(m_regressandMatrix.rows(),1);

  for (unsigned int i=0; i<m_regressandMatrix.cols();i++)
  {
    imiINFO("LOOCV index: "<<i);
    unsigned int index=0;
    for (unsigned int j=0; j<m_regressandMatrix.cols();j++)
    {
      if(i!=j)
      {
        newRegressand.set_column(index,m_regressandMatrix.get_column(j));
        newRegressor.set_column(index,m_regressorMatrix.get_column(j));
        index++;
      } else
      {
        measurement.set_column(0,m_regressorMatrix.get_column(j));
        reference.set_column(0,m_regressandMatrix.get_column(j));
      }
    }
    m_predictor->SetRegressandTrainingMatrix(newRegressand);
    m_predictor->SetRegressorTrainingMatrix(newRegressor);
    m_predictor->SetNumberOfComponents(maxNumOfComps);
    m_predictor->TrainEstimator();

    unsigned int idx=0;
    for (unsigned int comps=minNumOfComps; comps<=maxNumOfComps;comps++,idx++)
    {
      imiINFO("LOOCV prediction: "<<idx);
      m_predictor->SetNumberOfComponents(comps);
      VnlMatrixType predictedOutput;
      m_predictor->PredictOutput(measurement,predictedOutput);
      VnlMatrixType diffVec=(predictedOutput-reference);
      VnlMatrixType diff=diffVec.transpose()*diffVec;
      squaredDists[idx]+=diff[0][0];
    }
  }
  return;
}

float  imiMotionPredictionLOOCV::PerformLOOCVLinearKMLR(float ridgeParameterMin , float ridgeParameterMax,  unsigned int multiplicator)
{
  double minValue=std::numeric_limits<double>::max();
  float optRidgeParameter=0;

  std::vector<float> ridgeParameter;

  imiKernelMLRMotionPrediction* predictionFilter = imiKernelMLRMotionPrediction::New();
  imiLinearKernel* linearKernel;

  linearKernel=imiLinearKernel::New();
  predictionFilter->SetKernel(linearKernel);

  predictionFilter->SetMulticollinearityCheck( false );

  m_kernelPredictor=predictionFilter;

  std::vector<double> squaredDists;

  //ridgeParameter.push_back(.0);
  //squaredDists.push_back(.0);

  for (float i=ridgeParameterMin; i<=ridgeParameterMax;i*=multiplicator)
  {
    squaredDists.push_back(.0);
    ridgeParameter.push_back(i);
  }

  this->PerformLOOCV(ridgeParameter,squaredDists);

  for (unsigned int index=0; index<squaredDists.size();index++)
  {
    double diff=squaredDists[index];
    if(diff<minValue)
    {
      minValue=diff;
      optRidgeParameter=ridgeParameter[index];
    }
    imiINFO("LOOCV ridge parameter: "<<ridgeParameter[index]<<" diff: "<<diff);
  }

  return optRidgeParameter;

}

float  imiMotionPredictionLOOCV::PerformLOOCVGaussianKMLR(float ridgeParameterMin , float ridgeParameterMax,  unsigned int multiplicator)
{
  double minValue=std::numeric_limits<double>::max();
  float optRidgeParameter=0;

  std::vector<float> ridgeParameter;

  imiKernelMLRMotionPrediction* predictionFilter = imiKernelMLRMotionPrediction::New();
  //imiLinearKernel* linearKernel;
  imiGaussianKernel* gaussianKernel;


  //linearKernel=imiLinearKernel::New();
  gaussianKernel=imiGaussianKernel::New();
  gaussianKernel->SetSigma(gaussianKernel->ComputeSigmaValueArias(m_regressorMatrix));


  predictionFilter->SetKernel(gaussianKernel);

  predictionFilter->SetMulticollinearityCheck( false );

  m_kernelPredictor=predictionFilter;

  std::vector<double> squaredDists;

  //ridgeParameter.push_back(.0);
  //squaredDists.push_back(.0);

  for (float i=ridgeParameterMin; i<=ridgeParameterMax;i*=multiplicator)
  {
    squaredDists.push_back(.0);
    ridgeParameter.push_back(i);
  }

  this->PerformLOOCV(ridgeParameter,squaredDists);

  for (unsigned int index=0; index<squaredDists.size();index++)
  {
    double diff=squaredDists[index];
    if(diff<minValue)
    {
      minValue=diff;
      optRidgeParameter=ridgeParameter[index];
    }
    imiINFO("LOOCV ridge parameter: "<<ridgeParameter[index]<<" diff: "<<diff);
  }

  return optRidgeParameter;

}

void  imiMotionPredictionLOOCV::PerformLOOCV(std::vector<float> ridgeParameter, std::vector<double> &squaredDists)
{
  VnlMatrixType newRegressand;
  VnlMatrixType newRegressor;
  VnlMatrixType measurement;
  VnlMatrixType reference;



  newRegressand.set_size(m_regressandMatrix.rows(),m_regressandMatrix.cols()-1);
  newRegressor.set_size(m_regressorMatrix.rows(),m_regressorMatrix.cols()-1);
  measurement.set_size(m_regressorMatrix.rows(),1);
  reference.set_size(m_regressandMatrix.rows(),1);

  for (unsigned int i=0; i<m_regressandMatrix.cols();i++)
  {
    imiINFO("LOOCV index: "<<i);
    unsigned int index=0;
    for (unsigned int j=0; j<m_regressandMatrix.cols();j++)
    {
      if(i!=j)
      {
        newRegressand.set_column(index,m_regressandMatrix.get_column(j));
        newRegressor.set_column(index,m_regressorMatrix.get_column(j));
        index++;
      } else
      {
        measurement.set_column(0,m_regressorMatrix.get_column(j));
        reference.set_column(0,m_regressandMatrix.get_column(j));
      }
    }
    m_kernelPredictor->SetRegressandTrainingMatrix(newRegressand);
    m_kernelPredictor->SetRegressorTrainingMatrix(newRegressor);
    m_kernelPredictor->TrainLSEstimator();

    for (unsigned int idx=0; idx<ridgeParameter.size();idx++)
    {
      imiINFO("LOOCV prediction: "<<idx<<"->"<<i);
      m_kernelPredictor->SetTikhonovRegularizationParameter(ridgeParameter[idx]);
      VnlMatrixType predictedOutput;
      m_kernelPredictor->PredictOutput(measurement,predictedOutput);
      VnlMatrixType diffVec=(predictedOutput-reference);
      VnlMatrixType diff=diffVec.transpose()*diffVec;
      squaredDists[idx]+=diff[0][0];
    }
  }
  return;

}

}//namespace



