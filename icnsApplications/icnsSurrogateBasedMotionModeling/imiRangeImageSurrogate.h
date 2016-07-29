/** \file imiRangeImageSurrogate.h
 *
 *  \b Initial \b Author: Matthias Wilms, Maximilian Blendowski \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef IMIRANGEIMAGESURROGATE_H_
#define IMIRANGEIMAGESURROGATE_H_

// Standard includes:
#include <vector>

//ITK includes
#include "itkImageSliceIteratorWithIndex.h"      //iteration (outer loop: x, middle: z, inner: y)
#include "itkImageLinearIteratorWithIndex.h"     //parallel movement through the plane to store the y values

// Project includes:
#include "imiObject.h"
#include "imiITKImageTypeDefinitions.h"
#include "imiImageReader.h"
#include "imiImageWriter.h"
#include "imiMotionSurrogate.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"


namespace imi
{

class imiRangeImageSurrogate: public imiMotionSurrogate
{

public:

  /** Class macro needed because class inherits from imiObject */
  imiClassMacro( imiRangeImageSurrogate, imiObject )
  ;

  /** \brief Constructor */
  imiRangeImageSurrogate();

  /** \brief Destructor */
  virtual ~imiRangeImageSurrogate();

  /** Create an instance of the object, use Delete() to destroy. */
  static imiRangeImageSurrogate* New();

// ----------------------------
//   Typedefs:
// ----------------------------

//----------------------
// Methods:
//----------------------

//inherited methods needed to interact with MLR
  virtual bool ComputeMeasurements();
  virtual bool GetMeasurementMatrix( VnlMatrixType& matrix );

//methods needed to compute the matrices for MLR

//----------------------
// Setter:
//----------------------

  void SetImageListForChestComputing( std::vector<ImageType::Pointer> inputVector );
  void SetSkinThreshold( ImagePixelType skinThreshold )
  {
    m_skinThreshold = skinThreshold;
  }
  ;
  void SetROI( float lowX, float highX, float lowZ, float highZ )
  {
    m_roiXLow = lowX;
    m_roiXHigh = highX;
    m_roiZLow = lowZ;
    m_roiZHigh = highZ;
  }
  ;
  void SetSampling( unsigned int x, unsigned int z )
  {
    m_samplingX = x;
    m_samplingZ = z;
  }
  ;
  void SetOrientationAP()
  {
    m_orientationAP = true;
  }
  ;
  void SetOrientationPA()
  {
    m_orientationAP = false;
  }
  ;

  void SetUseBilateralFilter(bool useBilateralFilter)
  {
    m_useBilateralFilter=useBilateralFilter;
  }

  void SetAddKinectNoise(bool addKinectNoise)
  {
    m_addKinectNoise = addKinectNoise;
  }

  void SetSigmaSpatial(float sigmaSpatial)
  {
    m_sigmaSpatial=sigmaSpatial;
  }

  void SetSigmaIntens(float sigmaIntens)
  {
    m_sigmaIntens=sigmaIntens;
  }

  void SetCenterSamplingPoints(bool centerSamplingPoints)
  {
    m_centerSamplingPoints = centerSamplingPoints;
  }
  ;

  //camera distance in mm
  void SetCameraDistance(float cameraDistance)
  {
    m_cameraDistance=cameraDistance;
  }

//----------------------
// Getter:
//----------------------

  std::vector<VnlMatrixType> GetRangeImages()
  {
    return m_rangeImages;
  }
  ;

//----------------------
// Internal handling:
//----------------------


//methods needed to extract the chest border movement
  void PerformSurfaceTracking();

private:
  VnlMatrixType GenerateRangeImage( ImageType::Pointer inputImage, unsigned int phase );
  static MatrixValueType AddKinectNoise(MatrixValueType dist, MatrixValueType radius=0);

//methods needed to compute the

//----------------------
// Variables:
//----------------------
protected:

  ImageFloat3DTypePointer m_surfaceImage;
  std::vector<ImageType::Pointer> m_inputImages;
  int m_inputNumber;
  ImagePixelType m_skinThreshold;

  std::vector<VnlMatrixType> m_rangeImages;

  bool m_orientationAP;
  bool m_centerSamplingPoints;
  unsigned int m_samplingX;
  unsigned int m_samplingZ;
  float m_roiXLow;
  float m_roiXHigh;
  float m_roiZLow;
  float m_roiZHigh;

  bool m_useBilateralFilter;
  float m_sigmaSpatial;
  float m_sigmaIntens;

  bool m_addKinectNoise;
  float m_cameraDistance;

};

}

#endif /* IMIRANGEIMAGESURROGATE_H_ */
