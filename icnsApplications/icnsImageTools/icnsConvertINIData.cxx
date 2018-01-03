/** \file icnsConvertINIData.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2017 Department of Computational Neuroscience,
 *     University Medical Center Hamburg-Eppendorf
 *
 ****************************************************************************/

// System includes:
#include <string>
#include <iostream>
extern "C"
{
  #include "getopt.h"
}

// ITK includes:
#include <itkImageFileReader.h>
#include <itkRawImageIO.h>
#include <itkImageFileWriter.h>
#include <itkTileImageFilter.h>

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 2>                 ImageType2D;
typedef itk::ImageFileReader<ImageType2D>    ImageReaderType2D;
typedef itk::ImageFileWriter<ImageType2D>    ImageWriterType2D;

typedef itk::Image<short, 3>                 ImageType;
typedef itk::ImageFileReader<ImageType>      ImageReaderType;
typedef itk::ImageFileWriter<ImageType>      ImageWriterType;

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsConvertINIData -D <input dir> -l <number of slices> -O <output image> [...]" << std::endl;
  std::cout << std::endl;
  std::cout << "-D <input dir>                 Input directory to read from." << std::endl;
  std::cout << "-l <number of slices>          Number of slices in input dir." << std::endl;
  std::cout << "-O <output image>              Filename of the 3D image (= output image)." << std::endl;
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
  std::cout << "icnsConvertINIData                        " << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputDirectory = NULL;
  char* outputImageFilename = NULL;
  
  std::string inputImagePrefix = "I.0";
  unsigned int nSlices         = 60;
  
  // Reading parameters: 
  while( (c = getopt( argc, argv, "D:l:O:h?" )) != -1 )
  {
    switch( c )
    {
      case 'D':
        inputDirectory = optarg;
        std::cout << "  Input directory:                   " << inputDirectory << std::endl;
        break;
      case 'l':
        nSlices = atoi( optarg );
        std::cout << "  number of slices:                  " << nSlices << std::endl;
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << "  Output image filename:             " << outputImageFilename << std::endl;
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
  
  if( inputDirectory == NULL )
  {
    std::cerr << "ERROR: No input directory specified!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputImageFilename == NULL )
  {
    std::cerr << "ERROR: No output image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input images:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading images in input directory: " << std::endl;
  
  typedef itk::TileImageFilter< ImageType2D, ImageType > TilerType;
  TilerType::Pointer filter = TilerType::New();
  
  itk::FixedArray< short, 3 > layout;
    layout[0] = 1;
    layout[1] = 1;
    layout[2] = 0;
  
  filter->SetLayout( layout );
  
  for( int i = 0; i < nSlices; ++i )
  {
    // READING DATA:
  
    std::string inputFilename = inputDirectory + inputImagePrefix + std::to_string(35) + ".raw";
    /*
    std::string inputFilename = inputDirectory + inputImagePrefix + std::to_string(i+1);
    if( (i+1) < 10 )
    {
      inputFilename = inputDirectory + inputImagePrefix + "0" + std::to_string(i+1);
      //inputFilename = inputDirectory + inputImagePrefix + "0" + std::to_string(9);
    }
     */
    std::cout << inputFilename << " ... " << std::flush;
    
    itk::RawImageIO<unsigned short, 2>::Pointer io =
    itk::RawImageIO<unsigned short, 2>::New();
    
    io->SetFileName( inputFilename );
    
    unsigned int dim[2] = {256,256};
    double spacing[2]   = {1.0,1.0};
    double origin[2]    = {0.0,0.0};
    
    for( unsigned int i=0; i<2; i++ )
    {
      io->SetDimensions( i, dim[i] );
      io->SetSpacing( i, spacing[i] );
      io->SetOrigin( i, origin[i] );
    }
    io->SetHeaderSize( 7300 );
    io->SetByteOrderToLittleEndian();
    io->SetPixelType( itk::ImageIOBase::SCALAR );
    io->SetNumberOfComponents( 1 );
    
    ImageReaderType2D::Pointer inputImageReader = ImageReaderType2D::New();
    inputImageReader->SetFileName( inputFilename );
    
    try
    {
      inputImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading input reference image." << std::endl;
      std::cerr << excp.GetDescription() << std::endl;
      return EXIT_FAILURE;
    }
    
    // STEAKING DATA:
    
    ImageType2D::Pointer inputImage = inputImageReader->GetOutput();
    inputImage->DisconnectPipeline();
    
    filter->SetInput( i, inputImage );
    
    std::cout << " OK." << i << std::endl;
  }
  
  short filler = 0;
  filter->SetDefaultPixelValue( filler );
  filter->Update();

  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  //ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  //imageWriter->SetInput( histogramMatchingFilter->GetOutput() );
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetInput( filter->GetOutput() );
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
