/** \file imiSubspaceMotionPrediction.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/


#ifndef IMISUBSPACEMOTIONPREDICTION_H_
#define IMISUBSPACEMOTIONPREDICTION_H_

// Project includes:
#include "imiMotionPrediction.h"


namespace imi
{

class imiSubspaceMotionPrediction : public imiMotionPrediction
{

public:



  /* Class macro needed because class inherits from imiObject. */
  imiClassMacro(imiSubspaceMotionPrediction,imiMotionPrediction);

  imiSubspaceMotionPrediction(){ }

  virtual ~imiSubspaceMotionPrediction(){ }

  // ----------------------------
  //   Public Methods:
  // ----------------------------

  // --- Setter / Getter --------

  void SetNumberOfComponents( unsigned int numOfComponents )
  {
    m_numOfComponents=numOfComponents;
  }

protected:


  // ----------------------------
  //   Protected attributes.
  // ----------------------------

  VnlMatrixType m_regressorComponents;
  VnlMatrixType m_regressandComponents;
  VnlMatrixType m_regressorScores;
  VnlMatrixType m_regressandScores;
  VnlDiagMatrixType m_regressorEigenvalues;
  VnlDiagMatrixType m_regressandEigenvalues;

  unsigned int m_numOfComponents;

}; // end class

}



#endif /* IMISUBSPACEMOTIONPREDICTION_H_ */
