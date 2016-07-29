/** \file imiLinearKernel.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiLinearKernel_h
#define __imiLinearKernel_h


// vnl includes:
#include "vnl/vnl_vector.h"


// Project includes:
#include "imiObject.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"
#include "imiKernel.h"



namespace imi
{

class imiLinearKernel : public imiKernel
{

public:

  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiLinearKernel,imiKernel);

  /** Create an instance of the object, use Delete() to destroy. */
  static imiLinearKernel* New();

  // ----------------------------
  //   Public Methods:
  // ----------------------------

  virtual MatrixValueType Evaluate(const VnlVectorType& v1, const VnlVectorType& v2) const;

  VnlVectorType ComputePreImageOfProjection(VnlMatrixType& trainingData, VnlMatrixType& usedProjection,VnlMatrixType& usedEigenvectors,VnlDiagMatrixType& usedEigenvalues, bool centered) const;


protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  imiLinearKernel();

  /** \brief Destructor */
  virtual ~imiLinearKernel();


  // ----------------------------
  //   Protected attributes.
  // ----------------------------


}; // end class

}

#endif /* __imiLinearKernel_h */
