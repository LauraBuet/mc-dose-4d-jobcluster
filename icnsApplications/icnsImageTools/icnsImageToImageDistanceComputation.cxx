/** \file icnsImageToImageDistanceComputation.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2017 Department of Computational Neuroscience,
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
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>
#include <itkImageMaskSpatialObject.h>

#include <itkMeanSquaresImageToImageMetric.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
#include <itkMutualInformationImageToImageMetric.h>
#include <itkNormalizedMutualInformationHistogramImageToImageMetric.h>

#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkTranslationTransform.h>
#include <itkNormalizeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>

#include <itkMaskImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkNormalizeImageFilter.h>

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 3>                                   ImageType;
typedef itk::ImageFileReader<ImageType>                        ImageReaderType;
typedef itk::ImageFileWriter<ImageType>                        ImageWriterType;
typedef itk::Image<unsigned char, 3>                           MaskImageType;
typedef itk::Image<double, 3>                                  InternalImageType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> InterpolatorType;
typedef itk::NearestNeighborInterpolateImageFunction<InternalImageType, double> InternalImageTypeInterpolatorType;
typedef itk::TranslationTransform<double, 3>                   TransformType;
typedef itk::CastImageFilter<ImageType, MaskImageType>         CastFilterType;


// ---------------------------------------------------------------
// Declaration of additional routines:
// ---------------------------------------------------------------

double ComputeMSD( const ImageType::Pointer image1, const ImageType::Pointer image2 );

double ComputeMSD( const ImageType::Pointer image1,
                   const ImageType::Pointer image2,
                   const MaskImageType::Pointer mask,
                   const bool useMask );

double ComputeMI( const ImageType::Pointer image1,
                  const ImageType::Pointer image2,
                  const MaskImageType::Pointer mask,
                  const bool useMask );

double ComputeNCC( const ImageType::Pointer image1,
                   const ImageType::Pointer image2,
                   const MaskImageType::Pointer mask,
                   const bool useMask );

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsImageToImageDistanceComputation -I <input image 1> <input image 2> -M <mask image> -m <metric type>" << std::endl;
  std::cout << std::endl;
  std::cout << "-I <image 1> <image 2>         Filenames of images to be compared." << std::endl;
  std::cout << "-M <mask image>                Filename of mask image to restrict computation to." << std::endl;
  std::cout << "-m [0]                         Comparison distance measure option." << std::endl;
  std::cout << "                                 0: all (default)." << std::endl;
  std::cout << "                                 1: MSD." << std::endl;
  std::cout << "                                 2: NCC." << std::endl;
  std::cout << "                                 3: NMI." << std::endl;
  std::cout << "-l <logfile>                   Filename of logfile." << std::endl;
  std::cout << "-h                             Print this help." << std::endl;
  std::cout << std::endl;
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
  std::cout << "icnsImageToImageDistanceComputation"      << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputImage1Filename      = NULL;
  char* inputImage2Filename      = NULL;
  char* maskImageFilename        = NULL;
  char* logfileFilename          = NULL;
  
  bool useMask                   = false;
  
  short distanceMeasureChoice    = 0;
  bool computeMSD                = true;
  bool computeNCC                = true;
  bool computeMI                 = true;
  
  // Reading parameters: 
  while( (c = getopt( argc, argv, "I::M:m:l:?" )) != -1 )
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
      case 'm':
        distanceMeasureChoice = atoi( optarg );
        if( distanceMeasureChoice == 1 )
        {
          std::cout << "  Distance measure:                  only MSD" << std::endl;
          computeNCC = false;
          computeMI  = false;
        }
        else if( distanceMeasureChoice == 2 )
        {
          std::cout << "  Distance measure:                  only NCC" << std::endl;
          computeMSD = false;
          computeMI  = false;
        }
        else if( distanceMeasureChoice == 3 )
        {
          std::cout << "  Distance measure:                  only NMI" << std::endl;
          computeMSD = false;
          computeNCC  = false;
        }
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_SUCCESS;
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
  if( logfileFilename == NULL )
  {
    std::cerr << "No logfile given!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input data:"                      << std::endl;
  
  // Loading input intensity images:
  
  std::cout << "  First input image ... " << std::flush;
  
  ImageType::Pointer image1;
  ImageReaderType::Pointer inputImage1Reader = ImageReaderType::New();
  inputImage1Reader->SetFileName( inputImage1Filename );
  try
  {
    inputImage1Reader->Update();
    image1 = inputImage1Reader->GetOutput();
    
    ImageType::SpacingType imageSpacing;
    imageSpacing.Fill(1.0);
    
    image1->SetSpacing(imageSpacing);
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading first input image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;
  
  std::cout << "  Second input image ... " << std::flush;
  ImageType::Pointer image2;
  ImageReaderType::Pointer inputImage2Reader = ImageReaderType::New();
  inputImage2Reader->SetFileName( inputImage2Filename );
  try
  {
    inputImage2Reader->Update();
    image2 = inputImage2Reader->GetOutput();
    
    ImageType::SpacingType imageSpacing;
    imageSpacing.Fill(1.0);
    
    image2->SetSpacing(imageSpacing);
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading second input image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;
  
  // Loading mask image (if specified) and converting it to maskImageType
  // as required for itkImageToImageMetrics:

  ImageReaderType::Pointer maskImageReader;
  MaskImageType::Pointer maskImage;
  
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
    
    // Converting mask image to unsigned char-type:
    CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput( maskImageReader->GetOutput() );
    castFilter->Update();
    
    maskImage = castFilter->GetOutput();
    
    MaskImageType::SpacingType maskImageSpacing;
    maskImageSpacing.Fill(1.0);
    maskImage->SetSpacing(maskImageSpacing);
  }
  
  // If mask is used, use mask to mask input images:
  typedef itk::MaskImageFilter< ImageType, MaskImageType, ImageType > MaskFilterType;
  
  MaskFilterType::Pointer maskFilter1;
  MaskFilterType::Pointer maskFilter2;
  double voxelRatio = 1.0;
  
  if( useMask )
  {
    // Compute images:
    
    maskFilter1= MaskFilterType::New();
    maskFilter1->SetInput( image1 );
    maskFilter1->SetMaskImage( maskImage );
    maskFilter1->Update();
    image1 = maskFilter1->GetOutput();
  
    maskFilter2 = MaskFilterType::New();
    maskFilter2->SetInput( image2 );
    maskFilter2->SetMaskImage( maskImage );
    maskFilter2->Update();
    image2 = maskFilter2->GetOutput();
  
    // Compute ratio:
    double nTotalVoxel = 0;
    double nMaskVoxel = 0;
    
    itk::ImageRegionIterator<MaskImageType> maskIt( maskImage, maskImage->GetLargestPossibleRegion() );
    for( maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt )
    {
      nTotalVoxel++;
      if( maskIt.Get() != 0 ) nMaskVoxel++;
    }
    voxelRatio =  nTotalVoxel/nMaskVoxel;
  }
  
  // Generate logfile object:
  
  std::ofstream logfile;
  logfile.open( logfileFilename, std::ios::app );
  
  logfile << "ICNS IMAGE TO IMAGE METRIC STATISTICS" << std::endl;
  logfile << "Image 1:      " << inputImage1Filename << std::endl;
  logfile << "Image 2:      " << inputImage2Filename << std::endl;
  if( useMask ) logfile << "Mask image:   " << maskImageFilename << std::endl;
  
  // -------------------------------------------------------------
  // Computing distance between images:
  
  // MSD:
  if( computeMSD )
  {
    //double MSDValue = ComputeMSD_WO( image1, image2, maskImage, useMask);
    double MSDValue = ComputeMSD( image1, image2);
    if( useMask ) MSDValue *= voxelRatio;
    logfile << "MSD: " << MSDValue << std::endl;
    std::cout << "MSD: " << MSDValue << std::endl;
  }
  if( computeNCC )
  {
    double NCCValue = ComputeNCC( image1, image2, maskImage, useMask);
    logfile << "NCC: " << NCCValue << std::endl;
    std::cout << "NCC: " << NCCValue << std::endl;
  }
  if( computeMI )
  {
    double MIValue  = ComputeMI( image1, image2, maskImage, useMask);
    logfile << "MI: " << MIValue << std::endl;
    std::cout << "MI: " << MIValue << std::endl;
  }
  
  // -------------------------------------------------------------
  // Show logfile content and close logfile:
  
  logfile.close();
  
  std::ifstream writtenLogfile( logfileFilename );
  std::string str;
  std::string fileContents;
  while (std::getline( writtenLogfile, str ))
  {
    fileContents += str;
    fileContents.push_back('\n');
  }
  std::cout << fileContents;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Finished." << std::endl;
  std::cout << "==========================================" << std::endl;
  
  return EXIT_SUCCESS;
}

// Implementation of additional routines:

double ComputeMSD( const ImageType::Pointer image1, const ImageType::Pointer image2 )
{
  typedef itk::MeanSquaresImageToImageMetric < ImageType , ImageType > MetricType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double > InterpolatorType;
  typedef itk::TranslationTransform < double , 3 > TransformType;
  
  MetricType::Pointer metric = MetricType::New();
  TransformType::Pointer transform = TransformType::New();
  
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( image1 );
  
  metric->SetFixedImage( image1 );
  metric->SetMovingImage( image2 );
  metric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
  metric->SetTransform( transform );
  metric->SetInterpolator( interpolator );
  
  TransformType::ParametersType params(transform->GetNumberOfParameters());
  params.Fill(0.0);
  metric->Initialize();
  
  return metric->GetValue( params );
}


double ComputeMSD( const ImageType::Pointer image1,
                   const ImageType::Pointer image2,
                   const MaskImageType::Pointer mask,
                   const bool useMask )
{
  typedef itk::MeanSquaresImageToImageMetric < ImageType , ImageType > MetricType;
  typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double > InterpolatorType;
  typedef itk::TranslationTransform < double , 3 > TransformType;
  
  MetricType::Pointer metric = MetricType::New();
  TransformType::Pointer transform = TransformType::New();
  
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( image1 );
  
  metric->SetFixedImage( image1 );
  metric->SetMovingImage( image2 );
  metric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
  metric->SetTransform( transform );
  metric->SetInterpolator( interpolator );
  
  if( useMask )
  {
    typedef itk::ImageMaskSpatialObject<3>   MaskType;
    MaskType::Pointer spatialObjectMask = MaskType::New();
    spatialObjectMask->SetImage( mask ); // mask has to be unsigned char
    metric->SetFixedImageMask( spatialObjectMask );
  }
  
  TransformType::ParametersType params(transform->GetNumberOfParameters());
  params.Fill(0.0);
  metric->Initialize();
  
  return metric->GetValue( params );
} // end of ComputeMSD()


double ComputeMI( const ImageType::Pointer image1,
                  const ImageType::Pointer image2,
                  const MaskImageType::Pointer mask,
                  const bool useMask )
{
  typedef itk::NormalizedMutualInformationHistogramImageToImageMetric<InternalImageType, InternalImageType >    NormalizedMetricType;
  
  NormalizedMetricType::Pointer metric = NormalizedMetricType::New();
  
  typedef itk::NormalizeImageFilter< ImageType, InternalImageType > NormalizeFilterType;
  
  NormalizeFilterType::Pointer normalizeFilter1 = NormalizeFilterType::New();
  normalizeFilter1->SetInput(image1);
  normalizeFilter1->Update();
  
  NormalizeFilterType::Pointer normalizeFilter2 = NormalizeFilterType::New();
  normalizeFilter2->SetInput(image2);
  normalizeFilter2->Update();
  
  InternalImageTypeInterpolatorType::Pointer interpolator = InternalImageTypeInterpolatorType::New();
  //interpolator->SetInputImage( image1 );
  interpolator->SetInputImage( normalizeFilter1->GetOutput() );
  metric->SetInterpolator(interpolator);
  
  unsigned int numberOfHistogramBins = 6;
  NormalizedMetricType::HistogramType::SizeType histogramSize;
  histogramSize.SetSize(2);
  histogramSize[0] = numberOfHistogramBins;
  histogramSize[1] = numberOfHistogramBins;
  metric->SetHistogramSize( histogramSize );
  
  typedef itk::TranslationTransform<double , 3 > TransformType;
  TransformType::Pointer transform = TransformType::New();
  metric->SetTransform(transform);
  
  metric->SetFixedImage( normalizeFilter1->GetOutput() );
  metric->SetMovingImage( normalizeFilter2->GetOutput() );
  metric->SetFixedImageRegion( normalizeFilter1->GetOutput()->GetLargestPossibleRegion());
  
  TransformType::ParametersType parameters(transform->GetNumberOfParameters());
  parameters.Fill(0);
  
  if( useMask )
  {
    typedef itk::ImageMaskSpatialObject<3>   MaskType;
    MaskType::Pointer  spatialObjectMask = MaskType::New();
    spatialObjectMask->SetImage( mask ); // mask has to be unsigned char
    metric->SetFixedImageMask( spatialObjectMask );
  }

  metric->Initialize();
  return 1-metric->GetValue( parameters );
} // end of ComputeMI()


double ComputeNCC( const ImageType::Pointer image1,
                   const ImageType::Pointer image2,
                   const MaskImageType::Pointer mask,
                   const bool useMask )
{
  typedef itk::NormalizedCorrelationImageToImageMetric < ImageType , ImageType > MetricType;
  typedef itk::LinearInterpolateImageFunction<ImageType, double > InterpolatorType;
  typedef itk::TranslationTransform < double , 3 > TransformType;
  
  MetricType::Pointer metric = MetricType::New();
  TransformType::Pointer transform = TransformType::New();
  
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( image1 );
  
  metric->SetFixedImage( image1 );
  metric->SetMovingImage( image2 );
  metric->SetFixedImageRegion( image1->GetLargestPossibleRegion() );
  metric->SetTransform( transform );
  metric->SetInterpolator( interpolator );
  
  if( useMask )
  {
    typedef itk::ImageMaskSpatialObject<3>   MaskType;
    MaskType::Pointer  spatialObjectMask = MaskType::New();
    spatialObjectMask->SetImage( mask ); // mask has to be unsigned char
    metric->SetFixedImageMask( spatialObjectMask );
  }
  
  TransformType::ParametersType params(transform->GetNumberOfParameters());
  params.Fill(0.0);
  metric->Initialize();
  
  return metric->GetValue( params );
} // end of ComputeNCC()
