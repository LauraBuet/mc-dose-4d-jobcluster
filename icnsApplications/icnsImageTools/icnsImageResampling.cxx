/** \file icnsImageResampling.cxx
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
#include <itkBinaryThresholdImageFilter.h>
#include <itkIdentityTransform.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkLinearInterpolateImageFunction.h>

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 3>                 ImageType;
typedef itk::ImageFileReader<ImageType>      ImageReaderType;
typedef itk::ImageFileWriter<ImageType>      ImageWriterType;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> BinaryThresholdFilterType;
typedef itk::IdentityTransform<double, ImageType::ImageDimension> TransformType;
typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NNInterpolatorType;
typedef itk::BSplineInterpolateImageFunction<ImageType, double> BSplineInterpolatorType;
typedef itk::LinearInterpolateImageFunction<ImageType, double> LinearInterpolatorType;

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsImageResampling -I <input image> -O <output image> [-x  <x spacing>] [-y <y spacing>] [-z <z spacing>]" << std::endl;
  std::cout << std::endl;
  std::cout << "-I <input image>               Filename of the input image." << std::endl;
  std::cout << "-R <reference image>           Filename of optional reference image." << std::endl;
  std::cout << "-O <output image>              Filename of the output image." << std::endl;
  std::cout << "-x <x spacing>                 X spacing of output image (default = input x spacing)." << std::endl;
  std::cout << "-y <y spacing>                 Y spacing of output image (default = input y spacing)." << std::endl;
  std::cout << "-z <z spacing>                 Z spacing of output image (default = input z spacing)." << std::endl;
  std::cout << "-b                             Process binary input image." << std::endl;
  std::cout << "-l                             Use linear interpolation (default)." << std::endl;
  std::cout << "-n                             Use nearest neighbor interpolation." << std::endl;
  std::cout << "-s                             Use cubic B-spline interpolation." << std::endl;
  std::cout << "-h                             Print this help." << std::endl;
  std::cout << std::endl;
}


// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  std::cout << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsImageResampling"                        << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputImageFilename     = NULL;
  char* referenceImageFilename = NULL;
  char* outputImageFilename    = NULL;
  
  float outputXSpacing         = -1.0;
  float outputYSpacing         = -1.0;
  float outputZSpacing         = -1.0;
  
  int interpolationType        = 1; // 0: NN; 1: linear; 2: cubic B-Splines
  bool processBinaryImage      = false;
  
  
  // Reading parameters: 
  
  while( (c = getopt( argc, argv, "I:O:R:blnsx:y:z:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputImageFilename = optarg;
        std::cout << "  Input image filename:              " << inputImageFilename << std::endl;
        break;
      case 'R':
        referenceImageFilename = optarg;
        std::cout << "  Reference image filename:          " << referenceImageFilename << std::endl;
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << "  Output image filename:             " << outputImageFilename << std::endl;
        break;
      case 'x':
        outputXSpacing = atof( optarg );
        std::cout << "  Output x spacing:                  " << outputXSpacing << std::endl;
        break;
      case 'y':
        outputYSpacing = atof( optarg );
        std::cout << "  Output y spacing:                  " << outputYSpacing << std::endl;
        break;
      case 'z':
        outputZSpacing = atoi( optarg );
        std::cout << "  Output z spacing:                  " << outputZSpacing << std::endl;
        break;
      
      // additional parameters:
        
      case 'b':
        processBinaryImage = true;
        std::cout << "  Processing binary image:           TRUE" << std::endl;
        break;
      case 'l':
        interpolationType = 1;
        std::cout << "  Interpolation type:                LINEAR" << std::endl;
        break;
      case 'n':
        interpolationType = 0;
        std::cout << "  Interpolation type:                NEAREST NEIGHBOR" << std::endl;
        break;
      case 's':
        interpolationType = 2;
        std::cout << "  Interpolation type:                CUBIC B-SPLINES" << std::endl;
        break;
        
      // help data:
        
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
  
  if( inputImageFilename == NULL )
  {
    std::cerr << "ERROR: No input image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputImageFilename == NULL )
  {
    std::cerr << "ERROR: No output image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( (outputXSpacing == -1 ) && (outputYSpacing == -1 ) && (outputZSpacing == -1 ) && ( referenceImageFilename == NULL) )
  {
    std::cerr << "ERROR: no output spacing defined!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input images:
  // (1/2) Image to be resampled:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input image ... " << std::flush;

  ImageReaderType::Pointer inputImageReader = ImageReaderType::New();
  inputImageReader->SetFileName( inputImageFilename );
  try
  {
    inputImageReader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading input image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  
  // (2/2) Optional reference image:
  
  ImageReaderType::Pointer referenceImageReader;
  if( referenceImageFilename != NULL )
  {
    std::cout << "Loading reference image ... " << std::flush;
    
    referenceImageReader = ImageReaderType::New();
    referenceImageReader->SetFileName( referenceImageFilename );
    try
    {
      referenceImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading reference image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    
    std::cout << "OK." << std::endl;
  }
  
  // -------------------------------------------------------------
  // Resampling input image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Resampling input image ... " << std::endl;
  
  // Get input image spacing:
  
  ImageType::SpacingType inputSpacing = inputImageReader->GetOutput()->GetSpacing();
  ImageType::SizeType inputSize = inputImageReader->GetOutput()->GetLargestPossibleRegion().GetSize();
  
  // Set output image spacing, resampling factor and output size:
  
  ImageType::SpacingType outputSpacing( 1.0 );
  if( referenceImageFilename != NULL )
  {
    outputSpacing = referenceImageReader->GetOutput()->GetSpacing();
  }
  else
  {
    ( outputXSpacing > 0 ) ? ( outputSpacing[0] = outputXSpacing ) : ( inputImageReader->GetOutput()->GetSpacing()[0] );
    ( outputXSpacing > 0 ) ? ( outputSpacing[1] = outputYSpacing ) : ( inputImageReader->GetOutput()->GetSpacing()[1] );
    ( outputXSpacing > 0 ) ? ( outputSpacing[2] = outputZSpacing ) : ( inputImageReader->GetOutput()->GetSpacing()[2] );
  }
  
  float resamplingFactor[3];
  resamplingFactor[0] = inputSpacing[0] / outputSpacing[0];
  resamplingFactor[1] = inputSpacing[1] / outputSpacing[1];
  resamplingFactor[2] = inputSpacing[2] / outputSpacing[2];
  
  ImageType::SizeType outputSize;
  for( unsigned int i = 0; i < 3; i++ )
  {
    outputSize[i] = static_cast<ImageType::SizeValueType>( inputSize[i] * resamplingFactor[i] );
  }
  
  ImageType::PointType outputOrigin;
  outputOrigin = inputImageReader->GetOutput()->GetOrigin();

  // Pre-processing binary input data if requested:

  BinaryThresholdFilterType::Pointer thresholdFilter;
  if( processBinaryImage )
  {
    std::cout << "  Pre-processing binary image (rescaling dynamics to [0,128]) ... "  << std::flush;
    
    thresholdFilter = BinaryThresholdFilterType::New();
    thresholdFilter->SetOutsideValue( 0 );
    thresholdFilter->SetInsideValue( 128 );
    thresholdFilter->SetLowerThreshold( 1 );
    thresholdFilter->SetUpperThreshold( 255 );
    thresholdFilter->SetInput( inputImageReader->GetOutput() );
    
    try
    {
      thresholdFilter->Update();
    }
    catch( itk::ExceptionObject &err )
    {
      std::cerr << "Error during thresholding binary input image! Exception error: " << err << std::endl;
      return EXIT_FAILURE;
    }
    
    std::cout << "OK." << std::endl;
  }
  
  // Preparation done:
  
  std::cout << "  Resample factor: [" << resamplingFactor[0] << ", " <<resamplingFactor[1] << ", " << resamplingFactor[2] << "]" << std::endl;
  std::cout << "  Input spacing:   " << inputSpacing << std::endl;
  std::cout << "  Output spacing:  " << outputSpacing << std::endl;

  // Generate resampler itself:
  
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  resampler->SetDefaultPixelValue( 0 );
  
  TransformType::Pointer transform = TransformType::New();
  transform->SetIdentity();
  
  resampler->SetTransform( transform );
  
  // Set up interpolator:
  
  if( interpolationType == 0 )
  {
    // generate a nearest neighbor interpolator
    NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
    resampler->SetInterpolator( interpolator );
  }
  else if( interpolationType == 2 )
  {
    // generate a B-spline interpolator
    BSplineInterpolatorType::Pointer interpolator = BSplineInterpolatorType::New();
    resampler->SetInterpolator( interpolator );
  }
  else
  {
    // generate a linear interpolator
    LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
    resampler->SetInterpolator( interpolator );
  }
  
  // set (calculated) output spacing:
  resampler->SetOutputSpacing( outputSpacing );
  resampler->SetSize( outputSize );
  resampler->SetOutputOrigin( outputOrigin );
  resampler->SetOutputDirection( inputImageReader->GetOutput()->GetDirection() );

  if( processBinaryImage )
  {
    resampler->SetInput( thresholdFilter->GetOutput() );
  }
  else
  {
    resampler->SetInput( inputImageReader->GetOutput() );
  }
  
  try
  {
    resampler->Update();
  }
  catch( itk::ExceptionObject &err )
  {
    std::cerr << "ERROR: Cannot resample image! Exception error: " << err << std::endl;
    return EXIT_FAILURE;
  }
   
  // Post-process resampled binary image:
  
  BinaryThresholdFilterType::Pointer rescalingFilter;
  if( processBinaryImage )
  {
    std::cout << "Post-processing resampled binary segmentation image: rescaling to [0,1] ... "  << std::flush;
    
    rescalingFilter = BinaryThresholdFilterType::New();
    rescalingFilter->SetOutsideValue( 0 );
    rescalingFilter->SetInsideValue( 1 );
    rescalingFilter->SetLowerThreshold( 64 );
    rescalingFilter->SetUpperThreshold( 255 );
    rescalingFilter->SetInput( resampler->GetOutput() );
    
    try
    {
      rescalingFilter->Update();
    }
    catch( itk::ExceptionObject &err )
    {
      std::cerr << "Cannot execute rescaling! Exception error: " << err << std::endl;
      return EXIT_FAILURE;
    }
      
    std::cout << "OK." << std::endl;
  }
  
  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetFileName( outputImageFilename );
  
  if( processBinaryImage )
  {
    imageWriter->SetInput( rescalingFilter->GetOutput() );
  }
  else
  {
    imageWriter->SetInput( resampler->GetOutput() );
  }
  
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
  std::cout << "OK." << std::endl;
  
  return EXIT_SUCCESS;
}
