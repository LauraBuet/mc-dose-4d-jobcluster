/** \file icnsImageToImageDistanceComputation.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2016 Department of Computational Neuroscience,
 *     University Medical Center Hamburg-Eppendorf
 *
 ****************************************************************************/

// System includes:
#include <iostream>
extern "C"
{
#include "getopt.h"
}

// ITK includes:
//#include <itkAffineTransform.h>
//#include <itkBSplineInterpolateImageFunction.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageRegionIterator.h>
#include <itkImageToImageMetric.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkTranslationTransform.h>
#include <itkImageMaskSpatialObject.h>
//#include <itkNearestNeighborInterpolateImageFunction.h>
//#include <itkResampleImageFilter.h>
//#include <itkTransformFileReader.h>

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 3>                                               ImageType;
typedef itk::ImageFileReader<ImageType>                                    ImageReaderType;
typedef itk::ImageFileWriter<ImageType>                                    ImageWriterType;
typedef itk::ImageMaskSpatialObject<3>                                     MaskType;
typedef itk::NormalizedCorrelationImageToImageMetric<ImageType, ImageType> MutualInformationMetricType;
typedef MaskType::ImageType                                                SOImageType;
typedef itk::LinearInterpolateImageFunction<ImageType, double>             InterpolatorType;
typedef itk::TranslationTransform<double, 3>                               TransformType;
//typedef itk::BSplineInterpolateImageFunction<ImageType>         BSplineInterpolatorType;


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
    std::cout << "\n";
    std::cout << "Usage :\n";
    std::cout << "icnsImageToImageDistanceComputation -I <input image 1> <input image 2> -M <mask image> -o <metric type>\n\n";
    
    std::cout << "-I <image 1> <image 2>         Filenames of images to be compared.\n";
    std::cout << "-M <mask image>                Filename of mask image to restrict copmutation to.\n";
    std::cout << "-O <output image>              Filename of the output (=transformed) image.\n";
    std::cout << "-o [0]                         Comparison distance measure option.\n";
    std::cout << "                                 0: mutual information (default).\n";
    std::cout << "-l <logfile>                   Filename of logfile.\n";
  
    std::cout << "-h                             Print this help.\n";
    std::cout << "\n";
}


// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 2 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  std::cout << "========================================" << std::endl;
  std::cout << "icnsImageToImageDistanceComputation" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputImage1Filename      = NULL;
  char* inputImage2Filename      = NULL;
  char* maskImageFilename        = NULL;
  char* logfileFilename          = NULL;
  
  bool useMask                   = false;
  
  short distanceMeasureChoice    = 0;
  short backgroundIntensityValue = 0;
  short interpolationApproach    = 1; // 0: NN; 1: linear; 2: cubic B-Splines
  
  // Reading parameters: 
  while( (c = getopt( argc, argv, "I::M:o:?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
#ifdef __APPLE__
        std::cout << "  -- getopt with APPLE system --" << std::endl;
        inputImage1Filename = argv[optind-1];
        inputImage2Filename = argv[optind];
        optind = optind+1;
#elif __unix
        std::cout << "  -- getopt with unix system --" << std::endl;
        inputImage1Filename = argv[optind];
        inputImage2Filename = argv[optind+1];
        optind = optind+2;
#endif
        std::cout << "  Input image 1 filename:            " << inputImage1Filename << std::endl;
        std::cout << "  Input image 2 filename:            " << inputImage2Filename << std::endl;
        break;
      case 'M':
        maskImageFilename = optarg;
        useMask = true;
        std::cout << "  Mask image filename:               " << maskImageFilename << std::endl;
        break;
      case 'l':
        logfileFilename = optarg;
        std::cout << "  Logfile filename:                  " << logfileFilename << std::endl;
        break;
      case 'o':
        distanceMeasureChoice = atoi( optarg );
        if( distanceMeasureChoice == 0 )
        {
          std::cout << "  Distance measure:                  MUTUAL INFORMATION" << std::endl;
        }
        else
        {
          std::cout << "  Distance measure:                  MUTUAL INFORMATION" << std::endl;
        }
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_FAILURE;
        break;
      default:
        std::cout << "  Argument " << (char)c << " not processed!\n" << std::endl;
        return EXIT_FAILURE;
    }
  }
  
  // Check validity of arguments:
  
  if( inputImage1Filename == NULL )
  {
    std::cerr << "ERROR: No first input image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( inputImage2Filename == NULL )
  {
    std::cerr << "ERROR: No second input image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input data:"                      << std::endl;
  
  // Loading input intensity images:
  
  std::cout << "  First input image ... " << std::flush;
  
  ImageReaderType::Pointer inputImage1Reader = ImageReaderType::New();
  inputImage1Reader->SetFileName( inputImage1Filename );
  try
  {
    inputImage1Reader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading first input image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  std::cout << "  Second input image ... " << std::flush;
  
  ImageReaderType::Pointer inputImage2Reader = ImageReaderType::New();
  inputImage2Reader->SetFileName( inputImage2Filename );
  try
  {
    inputImage2Reader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading second input image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  
  // Loading mask image (if specified) and converting it into spatial object
  // mask as required for itkImageToImageMetrics:

  ImageReaderType::Pointer maskImageReader;
  SOImageType::Pointer maskImageSO;
  MaskType::Pointer maskSO;
  
  if( useMask )
  {
    // Reading mask image:
    
    std::cout << "  Mask image ... " << std::flush;
    
    maskImageReader = ImageReaderType::New();
    maskImageReader->SetFileName( maskImageFilename );
    try
    {
      maskImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading mask image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    
    std::cout << "OK." << std::endl;
    
    // Converting ImageType-mask image to ImageMaskSpatialObject::PixelType-mask image:
    
    maskImageSO = SOImageType::New();
    SOImageType::SizeType size = maskImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
    SOImageType::IndexType index = {{ 0, 0, 0 }};
    
    SOImageType::RegionType region;
    region.SetSize( size );
    region.SetIndex( index );
    
    maskImageSO->SetRegions( region );
    maskImageSO->Allocate( true ); // initialize buffer to zero
    
    typedef itk::ImageRegionIterator< SOImageType > imageSOIterator;
    
    imageSOIterator itSO( maskImageSO, region );
    itSO.GoToBegin();
    
    while( !itSO.IsAtEnd() )
    {
      itSO.Set( (unsigned char)(maskImageReader->GetOutput()->GetPixel( itSO.GetIndex() )) );
      ++itSO;
    }
    maskSO = MaskType::New();
    maskSO->SetImage( maskImageSO );
  }
  
  // -------------------------------------------------------------
  // Computing distance between images:
  
  // For mutual information: Acc. to the ITK documentation, the
  // images to be compared should be normalized to zero mean and
  // unit standard deviation.
  // TODO: do so!

  MutualInformationMetricType::Pointer mutualInformationMetric = MutualInformationMetricType::New();
  
  //mutualInformationMetric->SetFixedImageStandardDeviation(  0.4 );
  //mutualInformationMetric->SetMovingImageStandardDeviation( 0.4 );
  
  TransformType::Pointer transform = TransformType::New();
  
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( inputImage1Reader->GetOutput() );
  
  mutualInformationMetric->SetFixedImage( inputImage1Reader->GetOutput() );
  mutualInformationMetric->SetMovingImage( inputImage2Reader->GetOutput() );

  mutualInformationMetric->SetFixedImageRegion( inputImage1Reader->GetOutput()->GetLargestPossibleRegion() );
  mutualInformationMetric->SetTransform( transform );
  mutualInformationMetric->SetInterpolator( interpolator );
  if( useMask )
  {
    mutualInformationMetric->SetFixedImageMask( maskSO );
  }
  
  TransformType::ParametersType params(transform->GetNumberOfParameters());
  params.Fill(0.0);
  
  mutualInformationMetric->Initialize();
  
  
  std::cout << "  MI VALUE: " << mutualInformationMetric->GetValue( params ) << std::endl;
  
  /*
  // -------------------------------------------------------------
  // Applying transformation to input image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Applying transformation ... " << std::flush;

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  
  typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetTransform( affineTransform );
  resample->SetInput( inputImageReader->GetOutput() );
  resample->SetReferenceImage( inputImageReader->GetOutput() );
  resample->UseReferenceImageOn();
  resample->SetDefaultPixelValue( backgroundIntensityValue );

  if( interpolationApproach == 0 )
  {
    // generate a nearest neighbor interpolator
    NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
    resample->SetInterpolator( interpolator );
  }
  else if( interpolationApproach == 2 )
  {
    // generate a B-spline interpolator
    BSplineInterpolatorType::Pointer interpolator = BSplineInterpolatorType::New();
    resample->SetInterpolator( interpolator );
  }
  else
  {
    // generate a linear interpolator
    LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
    resample->SetInterpolator( interpolator );
  }

  resample->Update();
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetInput( resample->GetOutput() );
  imageWriter->SetFileName( outputImageFilename );
  try
  {
    imageWriter->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "   ERROR while writing warped target image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;*/
  
  return EXIT_SUCCESS;
}
