/** \file imiCCAMotionPrediction.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiCCAMotionPrediction_h
#define __imiCCAMotionPrediction_h





// Project includes:
#include "imiSubspaceMotionPrediction.h"
#include "imiKernelPCA.h"


namespace imi
{

class imiCCAMotionPrediction : public imiSubspaceMotionPrediction
{

public:



  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiCCAMotionPrediction,imiSubspaceMotionPrediction);

  /** Create an instance of the object, use Delete() to destroy. */
  static imiCCAMotionPrediction* New();



  // ----------------------------
  //   Public Methods:
  // ----------------------------

  // --- Setter / Getter --------



  void SetUseBothBases( bool useBothBases )
  {
    m_useBothBases=useBothBases;
  }

  virtual bool PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput );
  virtual bool TrainEstimator( void );


protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  imiCCAMotionPrediction();

  /** \brief Destructor */
  virtual ~imiCCAMotionPrediction();



  // ----------------------------
  //   Protected attributes.
  // ----------------------------




  //TODO: KernelPCA with linear kernel for regressor PCA
  VnlMatrixType m_regressorPCAComponents;
  VnlMatrixType m_regressorPCAmeanVector;
  imiRealType m_regressorPCAvariabilityThreshold;

  VnlMatrixType m_regressandPCAComponents;
  VnlMatrixType m_regressandPCAmeanVector;
  imiRealType m_regressandPCAvariabilityThreshold;


  unsigned int m_maxIters;
  imiRealType m_tolerance;
  bool m_useBothBases;


}; // end class

}

#endif /* __imiCCAMotionPrediction_h */

