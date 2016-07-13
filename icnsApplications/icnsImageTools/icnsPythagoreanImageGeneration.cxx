/** \file icnsPythagoreanImageGeneration
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2016 Department of Computational Neuroscience,
 *     University Medical Center Hamburg-Eppendorf
 *
 ****************************************************************************/

// System includes:
#include <string>
#include <fstream>
#include <vector>
extern "C"
{
#include "getopt.h"
}

// ITK includes:
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkAddImageFilter.h>
#include <itkShiftScaleImageFilter.h>

// Project includes: NONE so far

// Global typedefs:
typedef itk::Image<short, 3>                            ImageType;
typedef itk::ImageFileReader<ImageType>                 ImageReaderType;
typedef itk::ImageFileWriter<ImageType>                 ImageWriterType;
typedef itk::AddImageFilter<ImageType>                  AddImageFilterType;
typedef itk::ShiftScaleImageFilter<ImageType,ImageType> ScaleImageFilterType;

// Additional method declarations:
ImageType::Pointer ArithmeticMeanImageComputation( const std::vector<ImageType::Pointer> inputImages );
ImageType::Pointer GeometricMeanImageComputation( const std::vector<ImageType::Pointer> inputImages ){ return NULL; };
ImageType::Pointer HarmonicMeanImageComputation( const std::vector<ImageType::Pointer> inputImages ){ return NULL; };;


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "icnsPhythagoreanMeanImageGeneration -I <input file> <input file> [<input file> ... <input file>] <output file> [-o <option>]\n";
    
  std::cout << "-I <... input images ...>      Filenames of input images to be averaged.\n";
  std::cout << "-O <output image>              Filename of output image.\n";
  std::cout << "-o [0|1|2]                     Averaging option:\n";
  std::cout << "                                 0: arithmetic mean (default).\n";
  std::cout << "                                 1: geometric mean (NYI).\n";
  std::cout << "                                 2: harmonic mean (NYI).\n";
  
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
  
  std::cout << "========================================" << std::endl;
  std::cout << "icnsPhythagoreanMeanImageGeneration" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Reading parameters ... " << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  std::vector<std::string> inputImageFilenames;
  std::vector<ImageType::Pointer> inputImages;
  
  int paramCounter                   = 0;
  char* outputImageFilename          = NULL;
  unsigned int meanComputationOption = 0;
  
  // Actually reading parameters:
  
  while( (c = getopt( argc, argv, "I:O:o:h?" )) != -1 )
  {
    switch( c )
    {
      // Reading and interpreting params as filenames -- as long as the next sign isn't '-':
      case 'I':
        optind--;
        paramCounter = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, paramCounter++ )
        {
          inputImageFilenames.push_back( argv[optind] );
          std::cout << "  Input image [" << paramCounter << "] filename:  " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << "  Output image filename:     " << outputImageFilename << std::endl;
        break;
      case 'o':
        meanComputationOption = atoi( optarg );
        if( meanComputationOption == 0 )      std::cout << "  Mean image approach:       Arithmetic mean." << std::endl;
        else if( meanComputationOption == 1 ) std::cout << "  Mean image approach:       Geometric mean." << std::endl;
        else if( meanComputationOption == 2 ) std::cout << "  Mean image approach:       Harmonic mean." << std::endl;
        else
        {
          std::cerr << "  Mean image approach:       Unknown approach. Aborting computation." << std::endl;
          return EXIT_FAILURE;
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
  
  if( inputImageFilenames.size() < 2 )
  {
    std::cerr << "ERROR: Less than 2 image filenames specified!" << std::endl;
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
  std::cout << "Loading input images ... " << std::endl;
  
  ImageReaderType::Pointer inputImageReader;
  for( unsigned int i = 0; i < inputImageFilenames.size(); i++ )
  {
    ImageType::Pointer inputImage;
    inputImageReader = ImageReaderType::New();
    inputImageReader->SetFileName( inputImageFilenames[i].c_str() );
    inputImage = inputImageReader->GetOutput();
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
    inputImage->DisconnectPipeline();
    inputImages.push_back( inputImage );
    
    std::cout << "  Image [" << i << "] read." << std::endl;
  }
 
  // -------------------------------------------------------------
  // Computing mean image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Computing mean image " << std::flush;
  
  ImageType::Pointer outputImage;
  
  if( meanComputationOption == 0 )
  {
    std::cout << "[arithmetic mean] ... " << std::flush;
    outputImage = ArithmeticMeanImageComputation( inputImages );
  }
  else if( meanComputationOption == 1 )
  {
    std::cout << "[geometric mean] ... " << std::flush;
    outputImage = GeometricMeanImageComputation( inputImages );
  }
  else if( meanComputationOption == 2 )
  {
    std::cout << "[harmonic mean] ... " << std::flush;
    outputImage = HarmonicMeanImageComputation( inputImages );
  }
  else
  {
    std::cerr << "ERROR: computation option UNKNOWN. Aborting computing." << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetInput( outputImage );
  imageWriter->SetFileName( outputImageFilename );
  try
  {
    imageWriter->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "ERROR while writing output image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;

  std::cout << "==========================================" << std::endl;
  return EXIT_SUCCESS;
}

// ---------------------------------------------------------------
// Definition of the additional methods:
// ---------------------------------------------------------------

/** Computing image average; option: arithmetic mean. */

ImageType::Pointer ArithmeticMeanImageComputation( const std::vector<ImageType::Pointer> inputImages )
{
  // Adding all images (TODO: check for overflow!):
  
  for( unsigned int i = 1; i < inputImages.size(); i++ )
  {
    if( inputImages[i]->GetRequestedRegion().GetSize() != inputImages[0]->GetRequestedRegion().GetSize() )
    {
      std::cerr << "ERROR: Images have different size!" << std::endl;
      return NULL;
    }
    if( inputImages[i]->GetSpacing() != inputImages[0]->GetSpacing() )
    {
      std::cerr << "ERROR: Images have different spacing!" << std::endl;
      return NULL;
    }
 
    std::cout << "  Adding two Images ... " << std::endl;

    AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
    addFilter->InPlaceOn();
    addFilter->SetInput1( inputImages[0] );
    addFilter->SetInput2( inputImages[i] );
    addFilter->GraftOutput( inputImages[0] );
    try
    {
      addFilter->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << "ERROR: Processing failed wit exception:" << excep << std::endl;
      return NULL;
    }
    inputImages[0]->Graft( addFilter->GetOutput() );
    
    std::cout << "OK." << std::endl;
  }

  // Scaling output image:
  
  double scaleFactor = 1.0 / inputImages.size();

  std::cout << "  Scaling image by " << scaleFactor << " ... " << std::flush;
  ScaleImageFilterType::Pointer scaleFilter = ScaleImageFilterType::New();
  scaleFilter->SetInput( inputImages[0] );
  scaleFilter->SetScale( scaleFactor );
  try
  {
    scaleFilter->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "ERROR: Processing failed with exception:" << excep << std::endl;
    return NULL;
  }
  
  return scaleFilter->GetOutput();
}
