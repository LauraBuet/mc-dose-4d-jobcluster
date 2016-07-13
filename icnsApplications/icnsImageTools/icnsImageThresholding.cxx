/** \file icnsImageThresholding.cxx
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

// ITK incldues:
#include <itkBinaryThresholdImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 3>                 ImageType;
typedef itk::ImageFileReader<ImageType>      ImageReaderType;
typedef itk::ImageFileWriter<ImageType>      ImageWriterType;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> BinaryThresholdFilterType;


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
    std::cout << "\n";
    std::cout << "Usage :\n";
    std::cout << "icnsImageThresholding -I <input image> -O <output image> -l  <lower threshold> -u <upper threshold> -f <inside value> -b <outside value> \n\n";
    
    std::cout << "-I <input image>               Filename of the input image.\n";
    std::cout << "-O <output image>              Filename of the output image.\n";
    std::cout << "-l <lower threshold>           Lower threshold (default = min value).\n";
    std::cout << "-u <upper threshold>           Upper threshold (default = max value).\n";
    std::cout << "-f <inside>                    Inside output value (default = 1).\n";
    std::cout << "-b <outside>                   Outside output value (default = 0).\n";
  
    std::cout << "-h                             Print this help.\n";
    std::cout << "\n";
}


// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 4 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  std::cout << "==========================================" << std::endl;
  std::cout << "=======    icnsImageThresholding    ======" << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "READING parameters...\n" << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputImageFilename  = NULL;
  char* outputImageFilename = NULL;
  
  short outsideValue        = 0;
  short insideValue         = 1;
  short lowerThreshold      = -1;
  short upperThreshold      = -1;
  
  bool lowerThresholdSet    = false;
  bool upperThresholdSet    = false;
  
  // Reading parameters: 
  while( (c = getopt( argc, argv, "I:O:l:u:f:b:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputImageFilename = optarg;
        std::cout << "  Input image filename:              " << inputImageFilename << std::endl;
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << "  Output image filename:             " << outputImageFilename << std::endl;
        break;
      case 'l':
        lowerThreshold = atoi( optarg );
        lowerThresholdSet = true;
        std::cout << "  Lower threshold:                   " << lowerThreshold << std::endl;
        break;
      case 'u':
        upperThreshold = atoi( optarg );
        upperThresholdSet = true;
        std::cout << "  Upper threshold:                   " << upperThreshold << std::endl;
        break;
      case 'f':
        insideValue = atoi( optarg );
        std::cout << "  Foreground value:                  " << insideValue << std::endl;
        break;
      case 'b':
        outsideValue = atoi( optarg );
        std::cout << "  Background value:                  " << outsideValue << std::endl;
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
  
  if( inputImageFilename == NULL )
  {
    std::cerr << "No input image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputImageFilename == NULL )
  {
    std::cerr << "No output image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( !lowerThresholdSet && !upperThresholdSet)
  {
    std::cerr << "Lower AND upper threshold not set!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input image:
  
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
  
  // -------------------------------------------------------------
  // Thresholding image:
  
  if( !lowerThresholdSet )
  {
    lowerThreshold = itk::NumericTraits<ImageType::PixelType>::min();
  }
  if( !upperThresholdSet )
  {
    upperThreshold = itk::NumericTraits<ImageType::PixelType>::max();
  }
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Thresholding image ..." << std::endl;
  std::cout << "  Values < " << lowerThreshold << " --> " << outsideValue << std::endl;
  std::cout << "  Values > " << upperThreshold << " --> " << outsideValue << std::endl;
  std::cout << "  Inside values: " << insideValue << std::endl;
  
  BinaryThresholdFilterType::Pointer thresholdFilter = BinaryThresholdFilterType::New();
  thresholdFilter->SetOutsideValue( outsideValue );
  thresholdFilter->SetInsideValue( insideValue );
  thresholdFilter->SetLowerThreshold( lowerThreshold );
  thresholdFilter->SetUpperThreshold( upperThreshold );
  thresholdFilter->SetInput( inputImageReader->GetOutput() );
  try
  {
    thresholdFilter->Update();
  }
  catch( itk::ExceptionObject &err )
  {
    std::cerr << "Cannot execute itkBinaryThresholdFilter! Exception error: " << err << std::endl;
    return EXIT_FAILURE;
  }

  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetInput( thresholdFilter->GetOutput() );
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
  std::cout << "OK." << std::endl;
  
  return EXIT_SUCCESS;
}
