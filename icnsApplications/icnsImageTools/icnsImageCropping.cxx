/** \file icnsImageCropping.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2016 Department of Computational Neuroscience,
 *     University Medical Center Hamburg-Eppendorf
 *
 ****************************************************************************/

#include <iostream>
extern "C"
{
#include "getopt.h"
}

#include <itkExtractImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkChangeInformationImageFilter.h>

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 3>                          ImageType;
typedef itk::ImageFileReader<ImageType>               ImageReaderType;
typedef itk::ImageFileWriter<ImageType>               ImageWriterType;
typedef itk::ExtractImageFilter<ImageType, ImageType> ExtractImageFilterType;
typedef itk::ChangeInformationImageFilter<ImageType>  ChangeImageInformationFilterType;


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
    std::cout << "\n";
    std::cout << "Usage:\n";
    std::cout << "icnsImageCropping -I <input image> -O <output image> "
                 "-s <x-size> <y-size> <z-size> [-i <x-index> <y-index> <z-index>]\n\n";
    
    std::cout << "-I <input image>               Filename of input image.\n";
    std::cout << "-O <output image>              Filename of output image.\n";
    std::cout << "-i <x> <y> <z>                 Starting index of cropping region.\n";
    std::cout << "-s <x> <y> <z>                 Size of cropping region.\n";
  
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
  
  std::cout << "========================================" << std::endl;
  std::cout << "icnsImageCropping" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Reading parameters ... " << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputImageFilename  = NULL;
  char* outputImageFilename = NULL;
  
  ImageType::IndexType croppingRegionStartingIndex;
  croppingRegionStartingIndex.Fill( 0 );
  ImageType::SizeType croppingRegionSize;
  croppingRegionSize.Fill( 0 );
  
  // Reading parameters: 
  while( (c = getopt( argc, argv, "I:O:i:::s:::h?" )) != -1 )
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
      case 'i':
#ifdef __APPLE__
        std::cout << "  -- getopt with APPLE system --" << std::endl;
        croppingRegionStartingIndex[0] = atoi(argv[optind-1]);
        croppingRegionStartingIndex[1] = atoi(argv[optind]);
        croppingRegionStartingIndex[2] = atoi(argv[optind+1]);
        optind = optind+2;
#elif __unix
        std::cout << "  -- getopt with unix system --" << std::endl;
        croppingRegionStartingIndex[0] = atoi(argv[optind]);
        croppingRegionStartingIndex[1] = atoi(argv[optind+1]);
        croppingRegionStartingIndex[2] = atoi(argv[optind+2]);
        optind = optind+3;
#endif
        std::cout << "  Cropping region starting index:    " << croppingRegionStartingIndex << std::endl;
        break;
      case 's':
#ifdef __APPLE__
        std::cout << "  -- getopt with APPLE system --" << std::endl;
        croppingRegionSize[0] = atoi(argv[optind-1]);
        croppingRegionSize[1] = atoi(argv[optind]);
        croppingRegionSize[2] = atoi(argv[optind+1]);
        optind = optind+2;
#elif __unix
        std::cout << "  -- getopt with unix system --" << std::endl;
        croppingRegionSize[0] = atof(argv[optind]);
        croppingRegionSize[1] = atof(argv[optind+1]);
        croppingRegionSize[2] = atof(argv[optind+2]);
        optind = optind+3;
#endif
        std::cout << "  Cropping region size:              " << croppingRegionSize << std::endl;
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
  // Cropping image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Cropping image:" << std::endl;
  
  ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
  
  ImageType::RegionType desiredRegion;
  desiredRegion.SetSize( croppingRegionSize );
  desiredRegion.SetIndex( croppingRegionStartingIndex );
  
  std::cout << "  Extracted region: " << desiredRegion;
  
  extractImageFilter->SetInput( inputImageReader->GetOutput() );
  extractImageFilter->SetExtractionRegion( desiredRegion );
  
  try
  {
    extractImageFilter->Update();
  }
  catch( itk::ExceptionObject &err )
  {
    std::cerr << "ERROR: Cannot extract sub-image! Exception error: " << err << std::endl;
    return EXIT_FAILURE;
  }
  
  // Setting offset and starting index to zero.
  // Note: this is just necessary due to the ITK image file writer implementation,
  // which computes a new origin if the starting index of an image is not zero.
  
  ChangeImageInformationFilterType::Pointer changeImageInformationFilter = ChangeImageInformationFilterType::New();
  changeImageInformationFilter->SetInput( extractImageFilter->GetOutput() );
  
  ImageType::SpacingType inputSpacing = inputImageReader->GetOutput()->GetSpacing();
  ImageType::PointType outputOrigin;
  outputOrigin.Fill( 0 );
  
  outputOrigin[0] = (-1) * croppingRegionStartingIndex[0] * inputSpacing[0];
  outputOrigin[1] = (-1) * croppingRegionStartingIndex[1] * inputSpacing[1];
  outputOrigin[2] = (-1) * croppingRegionStartingIndex[2] * inputSpacing[2];
  
  changeImageInformationFilter->SetOutputOrigin( outputOrigin );
  changeImageInformationFilter->ChangeOriginOn();
  changeImageInformationFilter->Update();
  
  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetInput( changeImageInformationFilter->GetOutput() );
  imageWriter->SetFileName( outputImageFilename );
  try
  {
    imageWriter->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "ERROR while writing warped target image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;
  
  return EXIT_SUCCESS;
}
