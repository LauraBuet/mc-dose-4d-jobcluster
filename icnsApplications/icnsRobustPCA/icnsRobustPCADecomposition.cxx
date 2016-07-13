/** \file imiImageRPCADecomposition.cxx
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2015 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// System includes:
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
extern "C"
{
#include "getopt.h"
}

// ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// Project includes:
#include "icnsMatrixFunctions.h"
#include "icnsRobustPCA.h"

using namespace imi;

// Global typedefs:

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Info :\n";
  std::cout << "RPCA decomposition of matrix D (cols = vectorized images) into a low-rank matrix L and a sparse matrix S using the inexact augmented lagrangian multiplier method: "<<std::endl;
  std::cout << " min_{L,S} ||L||_* + lambda ||S||_1, s.t. ||D-L-S||_{fro} < eps"<<std::endl;
  std::cout << "\n\n";
  std::cout << "Usage :\n";
  std::cout << "imiImageRPCADecomposition -I <list of input images> -L <list of low-rank output images> -S <list of sparse output images> ...\n\n";
  std::cout << "Input: " << std::endl;
  std::cout << "-I  <list of input images>           Filenames of input images.\n";
  std::cout << "-M  <mask image>                     Filenames of an optional binary mask image.\n";
  std::cout << std::endl;
  std::cout << "Output: " << std::endl;
  std::cout << "-L  <list of input images>           Filenames of low rank output images.\n";
  std::cout << "-S  <list of input images>           Filenames of sparse output images.\n";
  std::cout << std::endl;
  std::cout << "Parameters: " << std::endl;
  std::cout << "-l  <lambda>                         Lambda parameter used for term balancing. (default = 1 /sqrt(max( rows(D), cols(D) )/)\n ";
  std::cout << "-i  <num iterations>                 Set maximum number of iterations. (default = 1000)\n ";
  std::cout << "-x  <dbg-level>                      Set debug level\n ";
  std::cout << "-h                                   Print this help.\n\n";
}

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
    PrintHelp();
    return 1;
  }

  std::cout << "================================================" << std::endl;
  std::cout << "icnsImageRPCADecomposition" << std::endl;
  std::cout << "================================================" << std::endl;
  std::cout << "Reading parameters...\n" << std::endl;

  // Typedefs.
  typedef short ImagePixelType;
  typedef itk::Image<ImagePixelType, 3> ImageType;
  typedef ImageType::Pointer ImagePointerType;
  
  typedef unsigned char BinarySegmentationPixelType;
  typedef itk::Image<BinarySegmentationPixelType, 3> BinarySegmentationType;
  typedef BinarySegmentationType::Pointer BinarySegmentationPointerType;
  
  typedef itk::ImageFileReader<ImageType> ImageFileReaderType;
  typedef itk::ImageFileReader<BinarySegmentationType> BinaryImageFileReaderType;
  typedef itk::ImageFileWriter<ImageType> ImageFileWriterType;
  
  // Initialize parameters with default values.
  std::vector<std::string> inputImageFilenames;
  std::vector<std::string> outputLRFilenames;
  std::vector<std::string> outputSparseFilenames;
  std::string maskFilename;

  std::vector<ImagePointerType> inputImages;
  BinarySegmentationPointerType maskImage;

  double lambda = 0;
  int maxIters = 1000;

  // Reading parameters:
  
  char c;
  int cnt = 0;

  while( (c = getopt( argc, argv, "I:L:S:M:l:i:x:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        //reads the input filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputImageFilenames.push_back( argv[optind] );
          std::cout << "  image [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'L':
        //reads the output LR filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputLRFilenames.push_back( argv[optind] );
          std::cout << "  LR [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'S':
        //reads the output sparse filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputSparseFilenames.push_back( argv[optind] );
          std::cout << "  sparse [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'M':
        optind--;
        maskFilename = argv[optind];
        std::cout << "  mask = " << maskFilename << std::endl;
        break;
      case 'l':
        optind--;
        lambda = atof( argv[optind] );
        std::cout << "  lambda = " << lambda << std::endl;
        break;
      case 'i':
        optind--;
        maxIters = atoi( argv[optind] );
        std::cout << "  max. Iterations = " << maxIters << std::endl;
        break;
      case 'x':
        optind--;
        //intVal = atoi( argv[optind] );
        //std::cout << "  Debug-Level:                     " << intVal << std::endl;
        //imiObject::SetGlobalDebugLevel( intVal );
        break;
      case 'h':
      case '?':
        PrintHelp();
        exit( 1 );
        break;
      default:
        std::cerr << "Argument "<<(char)c<<" not processed !\n" << std::endl;
        break;
    }
  }
  // Initialize time for log-output
  //imiTime::imiInitTime();

  if( inputImageFilenames.size() < 2 )
  {
    std::cerr << "At least two input images are needed!\n" << std::endl;
    return EXIT_FAILURE;
  }

  if( inputImageFilenames.size() != outputLRFilenames.size() || inputImageFilenames.size() != outputSparseFilenames.size() )
  {
    std::cerr << "The same number of filenames need be provided for input and output images!\n" << std::endl;
    return EXIT_FAILURE;
  }

  // Reading input images:
  
  ImageFileReaderType::Pointer inputImageReader;
  
  std::cout << std::endl;
  std::cout << "  :::Loading image data..." << std::endl;
  
  for( unsigned int i = 0; i < inputImageFilenames.size(); i++ )
  {
    ImagePointerType inputImage;
    inputImageReader = ImageFileReaderType::New();
    inputImageReader->SetFileName( inputImageFilenames[i].c_str() );
    inputImage = inputImageReader->GetOutput();
    inputImage->Update();
    
    inputImage->DisconnectPipeline();
    inputImages.push_back( inputImage );
  }

  if( !maskFilename.empty() )
  {
    BinaryImageFileReaderType::Pointer maskImageReader = BinaryImageFileReaderType::New();
    maskImageReader->SetFileName( maskFilename.c_str() );
    maskImage = maskImageReader->GetOutput();
    maskImage->Update();
  }

  //convert images to columns of the data matrix D
  imiRPCA::ArmaMatrixType D( inputImages[0]->GetLargestPossibleRegion().GetNumberOfPixels(), inputImages.size(), arma::fill::zeros );
  for( unsigned int i = 0; i < inputImages.size(); i++ )
  {
    imiMatrixFunctions::ArmaVectorType column( D.n_rows );
    if( !imiMatrixFunctions::ConvertImageToColumnVector( inputImages[i], maskImage, column ) )
    {
      std::cerr << "Can not convert images to matrix!" << std::endl;
      return EXIT_FAILURE;
    }
    D.col( i ) = column;
  }

  //set lambda to default value
  if( !std::fabs( lambda ) > 0 )
  {
    //see Candes et al., "Robust Principal Component Analysis?", Journal of the ACM, 2011 for an explanation
    lambda = 1 / sqrt( std::max( D.n_rows, D.n_cols ) );
  }

  imiRPCA *RPCAFilter = imiRPCA::New();
  RPCAFilter->SetD( D );
  RPCAFilter->SetLambda( lambda );
  RPCAFilter->SetMaxIterations( maxIters );

  if( !RPCAFilter->Execute() )
  {
    std::cerr << "RPCA decomposition failed!" << std::endl;
    return EXIT_FAILURE;
  }

  imiRPCA::ArmaMatrixType L = RPCAFilter->GetL();
  imiRPCA::ArmaMatrixType S = RPCAFilter->GetS();

  //write resulting images to disk
  for( unsigned int i = 0; i < outputLRFilenames.size(); i++ )
  {
    ImagePointerType outputImage;
    imiMatrixFunctions::ArmaVectorType column = L.col(i);
    if( !imiMatrixFunctions::ConvertColumnVectorToImage( inputImages[0],outputImage, maskImage, column ) )
    {
      std::cerr << "Can not convert matrix matrix to images!" << std::endl;
      return EXIT_FAILURE;
    }

    std::cout << "Saving LR image..." << std::endl;
    ImageFileWriterType::Pointer outputImageWriter = ImageFileWriterType::New();
    outputImageWriter->SetInput( outputImage );
    outputImageWriter->SetFileName( outputLRFilenames[i].c_str() );
    outputImageWriter->Write();
  }

  for( unsigned int i = 0; i < outputSparseFilenames.size(); i++ )
  {
    ImagePointerType outputImage;
    imiMatrixFunctions::ArmaVectorType column = S.col(i);
    if( !imiMatrixFunctions::ConvertColumnVectorToImage( inputImages[0],outputImage, maskImage, column ) )
    {
      std::cerr << "Can not convert matrix matrix to images!" << std::endl;
      return EXIT_FAILURE;
    }

    std::cout << "Saving sparse image..." << std::endl;
    ImageFileWriterType::Pointer outputImageWriter = ImageFileWriterType::New();
    outputImageWriter->SetInput( outputImage );
    outputImageWriter->SetFileName( outputSparseFilenames[i].c_str() );
    outputImageWriter->Write();
  }

  //RPCAFilter->Delete();

  std::cout << "RPCA execution finished." << std::endl;
  //imiINFO( " Process time needed:  " << imi::imiTime::getProcessTime() );
  //imiINFO( " Overall time needed:  " << imiTime::getMonotonicTime() <<"\n" );

  return EXIT_SUCCESS;
}
