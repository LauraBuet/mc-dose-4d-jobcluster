/** \file imiPCRMotionPrediction.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiPCRMotionPrediction_h
#define __imiPCRMotionPrediction_h


// Project includes:
#include "imiSubspaceMotionPrediction.h"


namespace imi
{

class imiPCRMotionPrediction : public imiSubspaceMotionPrediction
{

public:

  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiPCRMotionPrediction,imiSubspaceMotionPrediction);

  /** Create an instance of the object, use Delete() to destroy. */
  static imiPCRMotionPrediction* New();


  // ----------------------------
  //   Public Methods:
  // ----------------------------

  virtual bool PredictOutput( VnlMatrixType& regressorMeasurement, VnlMatrixType& predictedOutput );
  virtual bool TrainEstimator( void );

protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  imiPCRMotionPrediction();

  /** \brief Destructor */
  virtual ~imiPCRMotionPrediction();


  // ----------------------------
  //   Protected attributes.
  // ----------------------------


  VnlDiagMatrixType m_invEigenvalues;

  imiRealType m_variabilityThreshold;

}; // end class

}

#endif /* __imiPCRMotionPrediction_h */
