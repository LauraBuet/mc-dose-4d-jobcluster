/** \file imiPLSMotionPrediction.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 *  This code is partially based on (Python) code from
 *
 *  http://scikit-learn.org,
 *
 *  especially https://github.com/scikit-learn/scikit-learn/blob/master/sklearn/pls.py,
 *
 *  and implements Partial Least Squares (PLS).
 *
 *  @article{scikit-learn,
 *  title={{Scikit-learn: Machine Learning in Python }},
 *  author={Pedregosa, F. and Varoquaux, G. and Gramfort, A. and Michel, V.
 *         and Thirion, B. and Grisel, O. and Blondel, M. and Prettenhofer, P.
 *        and Weiss, R. and Dubourg, V. and Vanderplas, J. and Passos, A. and
 *        Cournapeau, D. and Brucher, M. and Perrot, M. and Duchesnay, E.},
 *  journal={Journal of Machine Learning Research},
 *  volume={12},
 *  pages={2825--2830},
 *  year={2011}
 *  }
 *
 ****************************************************************************/

#ifndef __imiPLSMotionPrediction_h
#define __imiPLSMotionPrediction_h

// Project includes:
#include "imiSubspaceMotionPrediction.h"


namespace imi
{

class imiPLSMotionPrediction : public imiSubspaceMotionPrediction
{

public:



  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiPLSMotionPrediction,imiSubspaceMotionPrediction);

  /** Create an instance of the object, use Delete() to destroy. */
  static imiPLSMotionPrediction* New();


  enum PLSMode
  {
    PLSModeA, // "classical" PLS
    PLSModeB  // CCA
  };

  enum PLSDeflationMode
  {
    PLSDeflationModeRegression, //asymmetric
    PLSDeflationModeCanonical //symmetric
  };



  // ----------------------------
  //   Public Methods:
  // ----------------------------


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
  imiPLSMotionPrediction();

  /** \brief Destructor */
  virtual ~imiPLSMotionPrediction();

  bool NIPALSInnerLoop(VnlMatrixType &regressorMatrix, VnlMatrixType &regressandMatrix,VnlVectorType &regressorWeights,VnlVectorType &regressandWeights,VnlVectorType &regressorScores, VnlVectorType &regressandScores);


  // ----------------------------
  //   Protected attributes.
  // ----------------------------

  unsigned int m_maxIters;
  imiRealType m_tolerance;
  bool m_useBothBases;


  PLSMode m_mode;
  PLSDeflationMode m_deflationMode;

}; // end class

}

#endif /* __imiPLSMotionPrediction_h */

