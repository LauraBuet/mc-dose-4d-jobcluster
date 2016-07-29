/** \file imiGaussianKernel.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiGaussianKernel_h
#define __imiGaussianKernel_h


// vnl includes:
#include "vnl/vnl_vector.h"


// Project includes:
#include "imiObject.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"
#include "imiKernel.h"

/***
 * Pre-image calculation using the method proposed by Rathi et al.
 *
 * Yogesh Rathi ; Samuel Dambreville ; Allen Tannenbaum;
 * Statistical shape analysis using kernel PCA.
 * Proc. SPIE 6064, Image Processing: Algorithms and Systems, Neural Networks, and Machine Learning, 60641B (February 17, 2006);
 * doi:10.1117/12.641417.
 *
 * The notation used here is from
 *
 * ROBUST TARGET LOCALIZATION AND SEGMENTATION USING STATISTICAL METHODS
 * Omar Arif, PhD Thesis, 2010, Georgia Institute of Technology
 *
 *
 *
 * Derivation of a feature space distance using the KPCA coefficients and an overview of sigma (for the gaussian kernel) determination methods:
 *
 * Appendix F of
 * Left ventricle functional analysis in 2D+t contrast echocardiography within an atlas-based deformable template model framework
 * Ramon Casero Canas, PhD Thesis, 2008, University of Oxford
 *
 *
 *
 *
 */



namespace imi
{

class imiGaussianKernel : public imiKernel
{

public:

  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiGaussianKernel,imiKernel);

  /** Create an instance of the object, use Delete() to destroy. */
  static imiGaussianKernel* New();

  // ----------------------------
  //   Public Methods:
  // ----------------------------

  virtual MatrixValueType Evaluate(const VnlVectorType& v1, const VnlVectorType& v2) const;

  MatrixValueType Evaluate(MatrixValueType dist) const;

  void SetSigma(MatrixValueType sigma) { m_sigma=sigma; };

  MatrixValueType GetSigma() {return m_sigma; };

  MatrixValueType ComputeSigmaValueArias(VnlMatrixType& inputVectors) const;

  MatrixValueType ComputeSigmaValueCasero(VnlMatrixType& inputVectors) const;

  MatrixValueType ComputeSigmaValueCremers1(VnlMatrixType& inputVectors) const;

  MatrixValueType ComputeSigmaValueCremers2(VnlMatrixType& inputVectors) const;

  VnlVectorType ComputePreImageOfProjection(VnlMatrixType& trainingData, VnlMatrixType& usedProjection,VnlMatrixType& usedEigenvectors,VnlDiagMatrixType& usedEigenvalues, bool centered) const;


protected:

  // ----------------------------
  //   Protected methods.
  // ----------------------------

  /** \brief Constructor */
  imiGaussianKernel();

  /** \brief Destructor */
  virtual ~imiGaussianKernel();


  // ----------------------------
  //   Protected attributes.
  // ----------------------------

  MatrixValueType m_sigma;

}; // end class

}

#endif /* __imiGaussianKernel_h */
