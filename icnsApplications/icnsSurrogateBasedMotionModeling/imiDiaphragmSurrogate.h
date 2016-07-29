/** \file imiDiaphragmSurrogate.h
 *
 *  \b Initial \b Author: Matthias Wilms\n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef IMIDIAPHRAGMSURROGATE_H_
#define IMIDIAPHRAGMSURROGATE_H_

// Project includes:
#include "imiObject.h"
#include "imiITKImageTypeDefinitions.h"
#include "imiImageReader.h"
#include "imiImageWriter.h"
#include "imiMotionSurrogate.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"
#include <vector>

namespace imi
{

class imiDiaphragmSurrogate: public imiMotionSurrogate
{

public:

  /** Class macro needed because class inherits from imiObject */
  imiClassMacro( imiDiaphragmSurrogate, imiObject )
  ;

  /** \brief Constructor */
  imiDiaphragmSurrogate();

  /** \brief Destructor */
  virtual ~imiDiaphragmSurrogate();

  /** Create an instance of the object, use Delete() to destroy. */
  static imiDiaphragmSurrogate* New();

// ----------------------------
//   Typedefs:
// ----------------------------

//----------------------
// Methods:
//----------------------

//inherited methods needed to interact with MLR
virtual bool ComputeMeasurements();
virtual bool GetMeasurementMatrix(VnlMatrixType& matrix)
{
  matrix = m_measurementMatrix;
  return true;
};


//----------------------
// Setter:
//----------------------

void   SetInputImages(std::vector<ImageType::Pointer> inputImages) {m_inputImages=inputImages;};
void   SetMaskImages(std::vector<ImageType::Pointer> maskImages) {m_maskImages=maskImages;};
void   SetCoronalSlicePosition(unsigned int coroSlice){m_coroSlice=coroSlice;};
void   SetStartingPoints(std::vector<unsigned int> x, std::vector<unsigned int> z){m_xStart=x;m_zStart=z;};
void   SetXSize(unsigned int xSize){m_xSize=xSize;};
void   SetLungThreshold(int threshold){m_threshold=threshold;};
void   SetHeadFirst(bool headFirst){m_headFirst=headFirst;};

//----------------------
// Variables:
//----------------------
protected:
std::vector<ImageType::Pointer> m_inputImages;
std::vector<ImageType::Pointer> m_maskImages;
unsigned int m_coroSlice;
std::vector<unsigned int> m_xStart;
std::vector<unsigned int> m_zStart;
unsigned int m_xSize;
unsigned int m_beltSize;
int m_threshold;
bool m_headFirst;


};

}

#endif /* IMIDIAPHRAGMSURROGATE_H_ */
