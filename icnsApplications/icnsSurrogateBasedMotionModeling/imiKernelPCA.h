/** \file imiKernelPCA.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiKernelPCA_h
#define __imiKernelPCA_h


// System includes:
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>

// vnl includes:
#include "vnl/vnl_vector.h"
#include "vnl/algo/vnl_svd.h"
#include "vnl/algo/vnl_symmetric_eigensystem.h"


// Project includes:
#include "imiObject.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"
#include "imiKernel.h"

/** \brief
 * Base class for kernel PCA.
 */


namespace imi
{

class imiKernelPCA : public imiObject
{

public:

  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiKernelPCA,imiObject);

  /** Create an instance of the object, use Delete() to destroy. */
  static imiKernelPCA* New();

  // ----------------------------
  //   Typedefs:
  // ----------------------------

  // ----------------------------
  //   Public Methods:
  // ----------------------------

  // --- Setter / Getter --------

  void SetInputMatrix( VnlMatrixType& inputMatrix )
  {
    m_inputMatrix = inputMatrix;
  }

  void SetKernel( imiKernel* kernel )
  {
    m_kernel = kernel;
  }

  void SetVariabilityThreshold(MatrixValueType thresh)
  {
    m_variabilityThreshold=thresh;
    m_numOfComponentsUsed=-1;
  }

  void SetNumOfComponentsUsed(MatrixValueType num)
  {
    m_variabilityThreshold=-1;
    m_numOfComponentsUsed=num;
  }

  void SetCentered(bool centered)
  {
    m_centered=centered;
  }

  VnlDiagMatrixType GetEigenvalues()
  {
    return m_eigenvalues;
  }

  bool ComputePrincipalComponents();

  VnlVectorType ProjectDataPointOntoEigenvectors(VnlVectorType& data);


  VnlVectorType ComputePreImageOfProjection(VnlVectorType& projection);

  unsigned int NumOfCompsNeededForVarThresh() const;



protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  imiKernelPCA();

  /** \brief Destructor */
  virtual ~imiKernelPCA();



  // ----------------------------
  //   Protected attributes.
  // ----------------------------



  VnlMatrixType m_inputMatrix;
  VnlMatrixType m_K; //Kernel matrix
  VnlMatrixType m_centeredK; //Centered kernel matrix
  VnlMatrixType m_H; //Centering matrix

  VnlMatrixType m_eigenvectors; //Eigenvectors of the centered kernel matrix
  VnlDiagMatrixType m_eigenvalues; //Eigenvalues of the centered kernel matrix

  bool m_centered; //Centered analysis?

  imiKernel* m_kernel;

  MatrixValueType m_variabilityThreshold;
  int m_numOfComponentsUsed;

}; // end class

}

#endif /* __imiKernelPCA_h */
