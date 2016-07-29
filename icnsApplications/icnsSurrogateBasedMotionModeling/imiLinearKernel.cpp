/** \file imiLinearKernel.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#include <cmath>


// Project includes
#include "imiLinearKernel.h"

namespace imi
{

// ----------------------------
//   Constructor / Destructor
// ----------------------------

imiLinearKernel* imiLinearKernel::New()
{
  return new imiLinearKernel();
}

imiLinearKernel::imiLinearKernel()
{
  // Set default parameters:
}

// imiLinearKernel::~imiLinearKernel()
imiLinearKernel::~imiLinearKernel()
{
}

// ----------------------------
//   Public methods
// ----------------------------

MatrixValueType imiLinearKernel::Evaluate(const VnlVectorType& v1, const VnlVectorType& v2) const
{
  return dot_product(v1,v2);
}


VnlVectorType imiLinearKernel::ComputePreImageOfProjection(VnlMatrixType& trainingData, VnlMatrixType& usedProjection,VnlMatrixType& usedEigenvectors,VnlDiagMatrixType& usedEigenvalues, bool centered) const
{
  //only centered version!!!!!!!!!!!!!!!!!!!!!!!!

    VnlVectorType numerator;


    numerator.set_size(trainingData.rows());
    numerator.fill(0.);


    for (unsigned int i=0;i<usedEigenvectors.cols();i++)
    {
        usedEigenvectors.set_column(i,(MatrixValueType)1./std::sqrt(usedEigenvalues[i])*usedEigenvectors.get_column(i));
    }

    VnlMatrixType onesVec;
    onesVec.set_size(trainingData.cols(),1);
    onesVec.fill(1./trainingData.cols());

    VnlMatrixType identMat;
    identMat.set_size(trainingData.cols(),trainingData.cols());
    identMat.fill(0.0);
    identMat.fill_diagonal(1.0);

    VnlMatrixType onesMat;
    onesMat.set_size(trainingData.cols(),trainingData.cols());
    onesMat.fill(1./trainingData.cols());

    VnlMatrixType gamma = onesVec+(identMat-onesMat)*usedEigenvectors*usedProjection;

    //VnlVectorType preImage=numerator;

    for (unsigned int i=0; i<trainingData.cols();i++)
        {

          numerator+=gamma[i][0]*trainingData.get_column(i);


        }

    VnlVectorType preImage=numerator;


    return preImage;
}

// ----------------------------
//   Protected / Private methods
// ----------------------------

}// namespace imi
