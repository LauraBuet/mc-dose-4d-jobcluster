/** \file imiMotionPredictionLOOCV.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/


#ifndef IMIMOTIONPREDICTIONLOOCV_H_
#define IMIMOTIONPREDICTIONLOOCV_H_

// Project includes:
#include "imiSubspaceMotionPrediction.h"
#include "imiKernelMLRMotionPrediction.h"




namespace imi{

class imiMotionPredictionLOOCV : public imiObject {

public:

  /** Class macro needed because class inherits from imiObject */
  imiClassMacro( imiMotionPredictionLOOCV, imiObject );

  /** \brief Constructor */
  imiMotionPredictionLOOCV();

  /** \brief Destructor */
  virtual ~imiMotionPredictionLOOCV();

  /** Create an instance of the object, use Delete() to destroy. */
  static imiMotionPredictionLOOCV* New();

// ----------------------------
//   Typedefs:
// ----------------------------

//----------------------
// Methods:
//----------------------

unsigned int  PerformLOOCVPCR(unsigned int minNumOfComps, unsigned int maxNumOfComps);
unsigned int  PerformLOOCVPLS(unsigned int minNumOfComps, unsigned int maxNumOfComps);
unsigned int  PerformLOOCVCCA(unsigned int minNumOfComps, unsigned int maxNumOfComps);
float  PerformLOOCVLinearKMLR(float ridgeParameterMin , float ridgeParameterMax,  unsigned int multiplicator);
float  PerformLOOCVGaussianKMLR(float ridgeParameterMin , float ridgeParameterMax,  unsigned int multiplicator);



//----------------------
// Setter:
//----------------------

void SetRegressandTrainingMatrix( VnlMatrixType& regressandObservations )
{
  m_regressandMatrix = regressandObservations;
}
void SetRegressorTrainingMatrix( VnlMatrixType& regressorObservations )
{
  m_regressorMatrix = regressorObservations;
}


//----------------------
// Variables:
//----------------------
protected:

void  PerformLOOCV(unsigned int minNumOfComps, unsigned int maxNumOfComps, std::vector<double> &squaredDists);
void  PerformLOOCV(std::vector<float> ridgeParameter, std::vector<double> &squaredDists);
void  PerformLOOCVGaussianKernel(std::vector<float> sigma, std::vector<double> &squaredDists);

VnlMatrixType m_regressorMatrix;
VnlMatrixType m_regressandMatrix;

imiSubspaceMotionPrediction* m_predictor;

//TODO Hack

imiKernelMLRMotionPrediction* m_kernelPredictor;
imiGaussianKernel* m_gaussianKernel;

};

} //endnamespace


#endif /* IMIMOTIONPREDICTIONLOOCV_H_ */
