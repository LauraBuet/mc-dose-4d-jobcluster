/** \file imiKernel.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiKernel_h
#define __imiKernel_h


// vnl includes:
#include "vnl/vnl_vector.h"


// Project includes:
#include "imiObject.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"
//#include "imiImageThreads.h"


namespace imi
{

class imiKernel : public imiObject
{

public:

  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiKernel,imiObject);

  /** Create an instance of the object, use Delete() to destroy. */
  //static imiKernelPCA* New();

  // ----------------------------
  //   Public Methods:
  // ----------------------------

  virtual MatrixValueType Evaluate(const VnlVectorType& v1, const VnlVectorType& v2) const = 0;

  virtual VnlVectorType ComputePreImageOfProjection(VnlMatrixType& trainingData, VnlMatrixType& usedProjection,VnlMatrixType& usedEigenvectors, VnlDiagMatrixType& usedEigenvalues, bool centered) const = 0;


protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  imiKernel();

  /** \brief Destructor */
  virtual ~imiKernel();


  // ----------------------------
  //   Protected attributes.
  // ----------------------------

}; // end class

}

#endif /* __imiKernel_h */
