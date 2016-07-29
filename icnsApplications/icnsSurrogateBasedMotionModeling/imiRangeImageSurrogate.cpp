/** \file imiRangeImageSurrogate.cpp
 *
 *  \b Initial \b Author: Matthias Wilms, Maximilian Blendowski \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// Project includes
#include "imiRangeImageSurrogate.h"

#include <limits>

#include "vnl/vnl_sample.h"

//#include "imiFastBilateralImageFilter/itkFastBilateralImageFilter.h"
#include "itkBilateralImageFilter.h"
#include "itkCastImageFilter.h"

using namespace imi;

// --------------------------------------------------------
//   Constructor / Destructor / No Copy-Constructor needed
// --------------------------------------------------------

imiRangeImageSurrogate::imiRangeImageSurrogate()
{
  m_centerSamplingPoints = false;
  m_skinThreshold = 500;
  m_orientationAP = true;
  m_samplingX = 1;
  m_samplingZ = 1;
  m_roiXLow = 0;
  m_roiXHigh = 1;
  m_roiZLow = 0;
  m_roiZHigh = 1;
  m_sigmaIntens = 100;
  m_sigmaSpatial = 1.0;
  m_useBilateralFilter = false;
  m_addKinectNoise = false;
  m_cameraDistance = 0;
}

imiRangeImageSurrogate::~imiRangeImageSurrogate()
{
}

imiRangeImageSurrogate* imiRangeImageSurrogate::New()
{
  return new imiRangeImageSurrogate();
}

//--------------------------------------------------------
// Methods:
//--------------------------------------------------------

//inherited methods:
bool imiRangeImageSurrogate::ComputeMeasurements()
{
  PerformSurfaceTracking();
  //Filling m_regressorMatrix (first x, then z ->0-511,0; 0-511,1;... row-wise
  // and then breathing phases column-wise)
  int x_dim = m_surfaceImage->GetLargestPossibleRegion().GetSize()[0];
  int z_dim = m_surfaceImage->GetLargestPossibleRegion().GetSize()[2];
  m_measurementMatrix = VnlMatrixType( x_dim * z_dim, m_inputNumber );
  ImageType::IndexType valueIndex;

  for( int z_iter = 0; z_iter < z_dim; z_iter++ )
  {
    valueIndex[2] = z_iter;
    for( int x_iter = 0; x_iter < x_dim; x_iter++ )
    {
      valueIndex[0] = x_iter;
      for( int y_iter = 0; y_iter < m_inputNumber; y_iter++ )
      {
        valueIndex[1] = y_iter;
        m_measurementMatrix[z_iter * x_dim + x_iter][y_iter] = m_surfaceImage->GetPixel( valueIndex );
      }
    }
  }
  return true;
}

bool imiRangeImageSurrogate::GetMeasurementMatrix( VnlMatrixType& regressorMatrix )
{
  regressorMatrix = m_measurementMatrix;
  return true;
}

//----------------------
// Setter:
//----------------------

void imiRangeImageSurrogate::SetImageListForChestComputing( std::vector<ImageType::Pointer> inputVector )
{
  /**
   * Sets the list of images for which the surface has to be computed and
   * allocates an output image with a fitting size
   */
  m_inputImages = inputVector;
  m_inputNumber = m_inputImages.size();
}

//----------------------
// Internal handling:
//----------------------
void imiRangeImageSurrogate::PerformSurfaceTracking()
{
  for( int z = 0; z < m_inputNumber; z++ )
  {
    m_rangeImages.push_back( GenerateRangeImage( m_inputImages[z], z ) );
  }

  /** again an iteration over the output image to check
   *  if there are valid y values for all breathing phases
   */

  for( unsigned int x = 0; x < m_rangeImages[0].cols(); x++ )
  {
    for( unsigned int z = 0; z < m_rangeImages[0].rows(); z++ )
    {
      bool negFound = false;
      for( unsigned int t = 0; t < m_rangeImages.size() && !negFound; t++ )
      {
        if( m_rangeImages[t].get( z, x ) < 0 )
        {
          negFound = true;
        }
      }

      if( negFound )
      {
        for( unsigned int t = 0; t < m_rangeImages.size() && !negFound; t++ )
        {
          m_rangeImages[t].put( z, x, 0 );
        }
      }
    }
  }

  //account for distance between camera and body
  if(m_cameraDistance >=0)
  {
    MatrixValueType minDist = m_inputImages[0]->GetLargestPossibleRegion().GetSize( 1 ) * m_inputImages[0]->GetSpacing()[1];
    for( unsigned int t = 0; t < m_rangeImages.size(); t++ )
    {
      for( unsigned int x = 0; x < m_rangeImages[t].cols(); x++ )
      {
        for( unsigned int z = 0; z < m_rangeImages[t].rows(); z++ )
        {
          if( m_rangeImages[t].get( z, x ) < minDist )
          {
            minDist = m_rangeImages[t].get( z, x );
          }
        }
      }
    }

    for( unsigned int t = 0; t < m_rangeImages.size(); t++ )
    {
      for( unsigned int x = 0; x < m_rangeImages[t].cols(); x++ )
      {
        for( unsigned int z = 0; z < m_rangeImages[t].rows(); z++ )
        {
            m_rangeImages[t].put( z, x, m_rangeImages[t].get( z, x ) - minDist + m_cameraDistance );
        }
      }
    }
  }

  if( m_addKinectNoise )
  {
    vnl_sample_reseed();
    for( unsigned int t = 0; t < m_rangeImages.size(); t++ )
    {
      for( unsigned int x = 0; x < m_rangeImages[t].cols(); x++ )
      {
        for( unsigned int z = 0; z < m_rangeImages[t].rows(); z++ )
        {
            MatrixValueType noise=this->AddKinectNoise(m_rangeImages[t].get( z, x ));
            //std::cout<<"Distance: "<<m_rangeImages[t].get( z, x )<<" Noise: "<<noise<<std::endl;
            m_rangeImages[t].put( z, x,  m_rangeImages[t].get( z, x )+noise);
        }
      }
    }
  }
}

VnlMatrixType imiRangeImageSurrogate::GenerateRangeImage( ImageType::Pointer inputImage, unsigned int plane )
{
  /**
   * Samples the skin surface according to the specified sampling rates within a region of interest
   */

  VnlMatrixType rangeImage;

  float roiXLow, roiXHigh, roiZLow, roiZHigh;
  roiXLow = roiXHigh = roiZLow = roiZHigh = 0;

  if( m_centerSamplingPoints )
  {
    unsigned int origRangeX = m_roiXHigh - m_roiXLow + 1;
    unsigned int origRangeZ = m_roiZHigh - m_roiZLow + 1;
    float corrX = origRangeX / ((float) 2 + m_samplingX);
    float corrZ = origRangeZ / ((float) 2 + m_samplingZ);
    roiXLow = m_roiXLow + corrX;
    roiXHigh = m_roiXHigh - corrX;
    roiZLow = m_roiZLow + corrZ;
    roiZHigh = m_roiZHigh - corrZ;
  }
  else
  {
    roiXLow = m_roiXLow;
    roiXHigh = m_roiXHigh;
    roiZLow = m_roiZLow;
    roiZHigh = m_roiZHigh;
  }

  if( m_useBilateralFilter )
  {

    typedef itk::CastImageFilter<ImageType, ImageFloat3DType> CastFilterType;
    CastFilterType::Pointer castFilter = CastFilterType::New();

    //castFilter->SetInput(inputImage);

    //typedef itk::FastBilateralImageFilter<ImageType,ImageType>  BilateralFilterType;
    typedef itk::BilateralImageFilter<ImageType, ImageType> BilateralFilterType;

    BilateralFilterType::Pointer bilateralFilter = BilateralFilterType::New();

    bilateralFilter->SetInput( inputImage );
    bilateralFilter->SetDomainSigma( m_sigmaSpatial );
    bilateralFilter->SetRangeSigma( m_sigmaIntens );
    bilateralFilter->SetAutomaticKernelSize( true );
    bilateralFilter->SetNumberOfRangeGaussianSamples( 100 );
    bilateralFilter->Update();

    ImageFloat3DTypePointer floatImage;
    inputImage = bilateralFilter->GetOutput();
    inputImage->DisconnectPipeline();
  }

  unsigned int rangeX = roiXHigh - roiXLow + 1;
  unsigned int rangeZ = roiZHigh - roiZLow + 1;

  /*std::cout<<"m_roiZHigh: "<<m_roiZHigh<<std::endl;
   std::cout<<"m_roiZLow: "<<m_roiZLow<<std::endl;
   std::cout<<"rangeZ: "<<rangeZ<<std::endl;*/

  rangeImage.set_size( m_samplingZ, m_samplingX );
  rangeImage.fill( -1 );

  float stepX = rangeX / (static_cast<float>( m_samplingX ) - 1);
  float stepZ = rangeZ / (static_cast<float>( m_samplingZ ) - 1);

  if (m_samplingX==1)
  {
    stepX=1;
  }

  if (m_samplingZ==1)
  {
    stepZ=1;
  }

  unsigned int imageMaxX = inputImage->GetLargestPossibleRegion().GetSize( 0 ) - 1;
  unsigned int imageMaxY = inputImage->GetLargestPossibleRegion().GetSize( 1 ) - 1;
  unsigned int imageMaxZ = inputImage->GetLargestPossibleRegion().GetSize( 2 ) - 1;

  ImageType::ValueType predValue = std::numeric_limits<ImageType::ValueType>::min();

  for( unsigned int x = 0; x < m_samplingX; x++ )
  {
    //std::cout<<"x: "<<x<<std::endl;
    for( unsigned int z = 0; z < m_samplingZ; z++ )
    {
      //std::cout<<x<<","<<z<<std::endl;
      //bilinear interpolation of the sampling position
      float realX = roiXLow + x * stepX;
      unsigned int xLeft = realX;
      unsigned int xRight = realX + 1;

      float realZ = roiZLow + z * stepZ;
      unsigned int zLower = realZ;
      unsigned int zUpper = realZ + 1;

      //out of bounds checks
      xRight = (xRight > imageMaxX) ? imageMaxX : xRight;
      zUpper = (zUpper > imageMaxZ) ? imageMaxZ : zUpper;

      bool notFound = true;

      for( unsigned int y = 0; y <= imageMaxY && notFound; y++ )
      {
        ImageType::IndexType leftUpperIndex, rightUpperIndex, leftLowerIndex, rightLowerIndex;

        leftUpperIndex[0] = xLeft;
        leftUpperIndex[1] = y;
        leftUpperIndex[2] = zUpper;

        rightUpperIndex[0] = xRight;
        rightUpperIndex[1] = y;
        rightUpperIndex[2] = zUpper;

        leftLowerIndex[0] = xLeft;
        leftLowerIndex[1] = y;
        leftLowerIndex[2] = zLower;

        rightLowerIndex[0] = xRight;
        rightLowerIndex[1] = y;
        rightLowerIndex[2] = zLower;

        ImageType::ValueType leftUpperValue, rightUpperValue, leftLowerValue, rightLowerValue;

        leftUpperValue = inputImage->GetPixel( leftUpperIndex );
        rightUpperValue = inputImage->GetPixel( rightUpperIndex );
        leftLowerValue = inputImage->GetPixel( leftLowerIndex );
        rightLowerValue = inputImage->GetPixel( rightLowerIndex );

        /*std::cout<<"leftUpperValue"<<leftUpperValue<<std::endl;
        std::cout<<"rightUpperValue"<<rightUpperValue<<std::endl;
        std::cout<<"leftLowerValue"<<leftLowerValue<<std::endl;
        std::cout<<"m_samplingZ"<<m_samplingZ<<std::endl;*/

        ImageType::ValueType upperInterpolatedValue, lowerInterpolatedValue;
        float subPositionX = -(realX - xLeft);

        //std::cout<<"subPositionX: "<<subPositionX<<std::endl;

        upperInterpolatedValue = leftUpperValue + (leftUpperValue - rightUpperValue) * subPositionX;
        lowerInterpolatedValue = leftLowerValue + (leftLowerValue - rightLowerValue) * subPositionX;

        float subPositionZ = -(realZ - zLower);

        //std::cout<<"subPositionZ: "<<subPositionZ<<std::endl;

        ImageType::ValueType currValue;

        currValue = lowerInterpolatedValue + (lowerInterpolatedValue - upperInterpolatedValue) * subPositionZ;
        //std::cout<<"currValue"<<currValue<<std::endl;

        if( currValue >= m_skinThreshold )
        {
          notFound = false;

          if( y > 0 )
          {

            //linear interpolation between pred-position and curr-position value
            /*std::cout<<"y: "<<y<<std::endl;
             std::cout<<"m_skinThreshold: "<<m_skinThreshold<<std::endl;
             std::cout<<"currValue: "<<currValue<<std::endl;
             std::cout<<"predValue: "<<predValue<<std::endl;*/
            /*std::cout<<"leftUpperValue: "<<leftUpperValue<<std::endl;
             std::cout<<"leftLowerValue: "<<leftLowerValue<<std::endl;
             std::cout<<"zLower: "<<zLower<<std::endl;
             std::cout<<"zUpper: "<<zUpper<<std::endl;
             std::cout<<"m_skinThreshold: "<<m_skinThreshold<<std::endl;*/
            MatrixValueType skinPos = (y - 1) + (static_cast<float>( m_skinThreshold ) - predValue) / (currValue - predValue);
            //std::cout<<"skinPos: "<<skinPos<<std::endl;
             /*std::cout<<"currValue: "<<currValue<<std::endl;
             std::cout<<"predValue: "<<predValue<<std::endl;*/

            rangeImage.put( z, x, skinPos * inputImage->GetSpacing()[1] );
            //std::cout<<"pos: "<<skinPos * inputImage->GetSpacing()[1]<<std::endl;
            //rangeImage.put( z, x, skinPos );
          }
          else
          {
            rangeImage.put( z, x, 0.0 );
          }

        }
        else
        {
          predValue = currValue;
        }

      }
    }
  }

  return rangeImage;
}

MatrixValueType imiRangeImageSurrogate::AddKinectNoise( MatrixValueType dist, MatrixValueType radius)
{
  //Noise model proposed by Olesen et al. in Real-time extraction of surface patches with associated uncertainties
  //by means of Kinect cameras in J Real-Time Image Pro
  //DOI 10.007/s11554-012-0261-x
  const MatrixValueType p00 = 2.344;
  const MatrixValueType p10 = -1.202e-2;
  const MatrixValueType p01 = -2.734e-3;
  const MatrixValueType p20 = 1.818e-5;
  const MatrixValueType p11 = 6.516e-6;
  const MatrixValueType p02 = 1.233e-6;
  MatrixValueType sigma=p00+(p10*radius)+(p01*dist)+(p20*radius*radius)+(p11*radius*dist)+(p02*dist*dist);
  //std::cout<<"Sigma: "<<sigma<<std::endl;
  MatrixValueType noise=vnl_sample_normal(0.0,sigma);
  return noise;
}
