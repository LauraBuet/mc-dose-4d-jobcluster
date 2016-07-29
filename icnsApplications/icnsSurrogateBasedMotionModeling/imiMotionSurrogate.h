/** \file imiMotionSurrogate.h
 *
 *  \b Initial \b Author: blendowski \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef IMIMOTIONSURROGATE_H_
#define IMIMOTIONSURROGATE_H_

// Project includes:
#include "imiObject.h"
#include "imiITKImageTypeDefinitions.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"


namespace imi{

class imiMotionSurrogate : public imiObject {

  /**
   * Abstract base class to ensure that all inheriting classes implement functions to
   * compute regressor matrices and also implement getter functions for these matrices
   */


public:

  // Class macro needed because class inherits from imiObject
  imiClassMacro( imiMotionSurrogate, imiObject );

  imiMotionSurrogate(){ }

  virtual ~imiMotionSurrogate(){ }

  // ----------------------------
  //   Typedefs:
  // ----------------------------

  //----------------------
  // Methods:
  //----------------------

  virtual bool ComputeMeasurements()                                =     0;
  virtual bool GetMeasurementMatrix(VnlMatrixType& matrix)      =     0;


  //----------------------
  // Variables:
  //----------------------

protected:
  //cols=measurements
  //rows=signal dimensions
  VnlMatrixType m_measurementMatrix;

};

}


#endif /* IMIMOTIONSURROGATE_H_ */
