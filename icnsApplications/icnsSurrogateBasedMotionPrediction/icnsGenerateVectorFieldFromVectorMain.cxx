/** \file imiGenerateVectorFieldFromVectorMain.cxx
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2014 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// System includes:
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
//#include "float.h"
extern "C"
{
#include "getopt.h"
}

// ITK includes:
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "icnsSurrogateBasedMotionPredictionTypeDefinitions.h"

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "imiGenerateVectorFieldFromVector ..." << std::endl;
  std::cout << std::endl;
  std::cout << "Generates a vector field from a MATLAB matrix/vector (.mat)." << std::endl;
  std::cout << "-I            Filenames of the input vector fields." << std::endl;
  std::cout << "-O            Filename of the output Matlab matrix." << std::endl;
  std::cout << "-R            Filename of the reference field." << std::endl;
  std::cout << "-M            Filename of the mask." << std::endl;
  std::cout << "-h            Print this help." << std::endl;
  std::cout << std::endl;
}

int main( int argc, char *argv[] )
{

  std::cout << "===========================================" << std::endl;
  std::cout << "===  imiGenerateVectorFieldFromVector  ====" << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << "Reading parameters ... " << std::endl;

  std::string outputFilename;
  std::string maskFilename;
  std::string refFieldFilename;
  std::string inputFilename;

  char c;
  int cnt = 0;


  // Reading parameters
  while( (c = getopt( argc, argv, "O:M:R:I:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        optind--;
        inputFilename = argv[optind];
        std::cout << "  input matrix = " << inputFilename << std::endl;
        break;
      case 'O':
        optind--;
        outputFilename = argv[optind];
        std::cout << "  output field = " << outputFilename << std::endl;
        break;
      case 'M':
        optind--;
        maskFilename = argv[optind];
        std::cout << "  mask = " << maskFilename << std::endl;
        break;
      case 'R':
        optind--;
        refFieldFilename = argv[optind];
        std::cout << "  reference field = " << refFieldFilename << std::endl;
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_SUCCESS;
      default:
        std::cerr << "Argument "<<(char)c<<" not processed!" << std::endl;
        return EXIT_FAILURE;
    }
  }

  // Plausibility checks? TODO
  // Loading input data:
  
  VnlMatrixType inputMatrix;
  DisplacementFieldType::Pointer outputField;

  icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileRealType(inputMatrix,inputFilename);

  // std::cout<<"inputMatrix ("<<inputMatrix.rows()<<","<<inputMatrix.cols()<<"): "<<inputMatrix<<std::endl;

  icnsSurrogateBasedMotionPredictionHelpers::GenerateVectorFieldFromVnlVector(inputMatrix,outputField,refFieldFilename,maskFilename);

  typedef itk::ImageFileWriter< DisplacementFieldType > FieldWriterType;
  FieldWriterType::Pointer displacementFieldWriter = FieldWriterType::New();

  displacementFieldWriter->SetInput( outputField);
  displacementFieldWriter->SetFileName( outputFilename );
  try
  {
    displacementFieldWriter->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << " Save failed with exception:" << excep << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "imiGenerateVectorFieldFromVector FINISHED." << std::endl;
  std::cout << "==========================================\n" << std::endl;
  
  return EXIT_SUCCESS;
}
