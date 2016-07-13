/** \file icnsImageHistogramMatching.cxx
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
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkHistogramMatchingImageFilter.h>

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 3>                 ImageType;
typedef itk::ImageFileReader<ImageType>      ImageReaderType;
typedef itk::ImageFileWriter<ImageType>      ImageWriterType;
typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> HistogramMatchingFilterType;

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
    std::cout << "\n";
    std::cout << "Usage :\n";
    std::cout << "icnsImageHistogramMatching -R <reference image> -T <target image> -O <output image> [...]\n\n";
    
    std::cout << "-R <reference image>           Filename of input reference image.\n";
    std::cout << "-T <target image>              Filename of input target image.\n";
    std::cout << "-O <output image>              Filename of the transformed target image (= output image).\n";

    std::cout << "-h                             Print this help.\n";
    std::cout << "\n";
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
  
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsImageHistogramMatching" << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "Reading parameters...\n" << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputReferenceImageFilename  = NULL;
  char* inputTargetImageFilename  = NULL;
  char* outputImageFilename = NULL;
  
  // Reading parameters: 
  while( (c = getopt( argc, argv, "R:T:O:h?" )) != -1 )
  {
    switch( c )
    {
      case 'R':
        inputReferenceImageFilename = optarg;
        std::cout << "  Input reference image filename:    " << inputReferenceImageFilename << std::endl;
        break;
      case 'T':
        inputTargetImageFilename = optarg;
        std::cout << "  Input target image filename:       " << inputTargetImageFilename << std::endl;
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << "  Output image filename:             " << outputImageFilename << std::endl;
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
  
  if( inputReferenceImageFilename == NULL )
  {
    std::cerr << "ERROR: No input reference image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( inputTargetImageFilename == NULL )
  {
    std::cerr << "ERROR: No input target image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputImageFilename == NULL )
  {
    std::cerr << "ERROR: No output image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input image:
  // (1/2): reference image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input reference image ... " << std::flush;

  ImageReaderType::Pointer inputReferenceImageReader = ImageReaderType::New();
  inputReferenceImageReader->SetFileName( inputReferenceImageFilename );
  try
  {
    inputReferenceImageReader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading input reference image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  
  // (2/2): target image
  
  std::cout << "Loading input target image ... " << std::flush;
  
  ImageReaderType::Pointer inputTargetImageReader = ImageReaderType::New();
  inputTargetImageReader->SetFileName( inputTargetImageFilename );
  try
  {
    inputTargetImageReader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading input target image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Histogram matching:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Performing histogram matching ... " << std::flush;
  
  HistogramMatchingFilterType::Pointer histogramMatchingFilter = HistogramMatchingFilterType::New();
  histogramMatchingFilter->SetSourceImage( inputTargetImageReader->GetOutput() );
  histogramMatchingFilter->SetReferenceImage( inputReferenceImageReader->GetOutput() );
  histogramMatchingFilter->SetNumberOfHistogramLevels( 1024 );
  histogramMatchingFilter->SetNumberOfMatchPoints( 7 );
  histogramMatchingFilter->ThresholdAtMeanIntensityOn();
  
  try
  {
    histogramMatchingFilter->Update();
  }
  catch( itk::ExceptionObject& err )
  {
    std::cerr << "ERROR: Could not match input images! ITK error: " << err;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  std::cout << "Histogram-Matching-Filter:" << std::endl;
  histogramMatchingFilter->Print( std::cout );

  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetInput( histogramMatchingFilter->GetOutput() );
  imageWriter->SetFileName( outputImageFilename );
  try
  {
    imageWriter->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "   ERROR while writing transformed target image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;
  
  return EXIT_SUCCESS;
}
