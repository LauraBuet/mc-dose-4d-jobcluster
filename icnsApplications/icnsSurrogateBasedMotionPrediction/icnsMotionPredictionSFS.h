/** \file imiMotionPredictionSFSMethods.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef IMIMOTIONPREDICTIONSFS_H_
#define IMIMOTIONPREDICTIONSFS_H_

// Project includes:
#include "icnsSurrogateBasedMotionPredictionTypeDefinitions.h"
#include "icnsMLRMotionPrediction.h"

#include <set>
#include <vector>

namespace imi
{

/*
 * A class for sequential forward selection-based determination of optimal combinations
 * of surrogate signals for motion prediction
 */

class imiMotionPredictionSFS
{

public:

  /** Class macro needed because class inherits from imiObject */
  //imiClassMacro( imiMotionPredictionSFS, imiObject );

  /** \brief Constructor */
  imiMotionPredictionSFS();

  /** \brief Destructor */
  virtual ~imiMotionPredictionSFS();

  /** Create an instance of the object, use Delete() to destroy. */
  static imiMotionPredictionSFS* New();

  void SetRegressandMatrix(VnlMatrixType* regressandMatrix )
  {
    m_regressandMatrix = regressandMatrix;
  }

  void SetRidgeParameters( float ridgeParameterMin, float ridgeParameterMax, float ridgeMultiplicator )
  {
    m_ridgeParameterMin = ridgeParameterMin;
    m_ridgeParameterMax = ridgeParameterMax;
    m_ridgeMultiplicator = ridgeMultiplicator;
  }

  void SetMaxNumOfIterations( unsigned int maxNumOfIterations )
  {
    m_maxNumOfIterations=maxNumOfIterations;
  }

  //Each regressor matrix should represent one surrogate
  void AddRegressorMatrix(VnlMatrixType* regressorMatrix)
  {
    m_regressorMatrices.push_back(regressorMatrix);
  }

  std::set<unsigned int> GetOptimalIndices()
  {
    if(!m_SFSPerformed)
    {
    //  imiWARNING("Set of optimal indices is empty as SFS has not been performed so far. Please call this->PerformSFS() to compute optimal indices!")
    }
    /*std::cout << "Get Optimal indices: ";
    for( std::set<unsigned int>::iterator it2 = m_optimalIndices.begin(); it2 != m_optimalIndices.end(); ++it2 )
    {
      std::cout << *it2 << "," << std::endl;
    }*/

    return m_optimalIndices;
  }

  float GetRidgeParameter()
  {
    return m_ridgeParameter;
  }

  VnlMatrixType* GetOptimalSubset()
  {
    return m_optimalSubset;
  }

  void PerformSFS();

  static VnlMatrixType* CombineMatrices(VnlMatrixType* matrixA, VnlMatrixType* matrixB);

protected:
  double CalculateSSD(VnlMatrixType* predictions);

  std::vector<VnlMatrixType*> m_regressorMatrices;
  VnlMatrixType* m_regressandMatrix;
  std::set<unsigned int> m_optimalIndices;
  VnlMatrixType* m_optimalSubset;

  bool m_SFSPerformed;

  unsigned int m_maxNumOfIterations;

  double m_optimalSSD;

//TODO Hack; should be generalized to ensure that every 'prediction approach' can be used (abstract base class????)
  icnsMLRMotionPrediction* m_predictor;
  std::vector<icnsMLRMotionPrediction*> m_predictors;
  float m_ridgeParameterMin;
  float m_ridgeParameterMax;
  float m_ridgeMultiplicator;
  float m_ridgeParameter;

};

} //endnamespace

#endif /* IMIMOTIONPREDICTIONSFS_H_ */
