/** \file imiBeltSurrogate.h
 *
 *  \b Initial \b Author: blendowski \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef IMICHESTBELTSURROGATE_H_
#define IMICHESTBELTSURROGATE_H_

// Project includes:
#include "imiObject.h"
#include "imiITKImageTypeDefinitions.h"
#include "imiMotionSurrogate.h"
#include <vector>



namespace imi{

class imiBeltSurrogate : public imiMotionSurrogate {

public:

  /** Class macro needed because class inherits from imiObject */
  imiClassMacro( imiBeltSurrogate, imiObject );

  /** \brief Constructor */
  imiBeltSurrogate();

  /** \brief Destructor */
  virtual ~imiBeltSurrogate();

  /** Create an instance of the object, use Delete() to destroy. */
  static imiBeltSurrogate* New();

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

void   SetInputImages(std::vector<BinarySegmentationPointerType> inputImages) {m_inputImages=inputImages;};
void   SetBeltPosition(unsigned int beltPosition){m_beltPosition=beltPosition;};
void   SetBeltSize(unsigned int beltSize){m_beltSize=beltSize;};

void SetROI( unsigned int lowX, unsigned int highX, unsigned int lowY, unsigned int highY )
{
  m_roiXLow = lowX;
  m_roiXHigh = highX;
  m_roiYLow = lowY;
  m_roiYHigh = highY;
};

//----------------------
// Variables:
//----------------------
protected:
std::vector<BinarySegmentationPointerType> m_inputImages;
unsigned int m_beltPosition;
unsigned int m_beltSize;
unsigned int m_roiXLow;
unsigned int m_roiXHigh;
unsigned int m_roiYLow;
unsigned int m_roiYHigh;

};

} //endnamespace

#endif /* IMICHESTBELTSURROGATE_H_ */
