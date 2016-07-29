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
#include "float.h"
extern "C"
{
#include "getopt.h"
}

// ITK includes:
#include "imiImageReader.h"
#include "imiImageWriter.h"

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenerateVectorFieldFromVector ... \n\n";
  std::cout << "Generates a vector field from a MATLAB matrix/vector (.mat). \n";
  std::cout << "-I            Filenames of the input vector fields.\n";
  std::cout << "-O            Filename of the output Matlab matrix. \n";
  std::cout << "-R            Filename of the reference field. \n";
  std::cout << "-M            Filename of the mask. \n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{

  // Initialize time for log-output
  imiTime::imiInitTime();

  imiINFO( "===========================================" );
  imiINFO( "===  imiGenerateVectorFieldFromVector  ====" );
  imiINFO( "===========================================" );
  imiINFO( "Reading parameters ...\n" );

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
        exit( 1 );
        break;
      default:
        imiERROR( "Argument "<<(char)c<<" not processed !\n" );
        break;
    }
  }

  /*
   * Loading input data
   */
  VnlMatrixType inputMatrix;

  DisplacementFieldType::Pointer outputField;

  imiSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileRealType(inputMatrix,inputFilename);

//  std::cout<<"inputMatrix ("<<inputMatrix.rows()<<","<<inputMatrix.cols()<<"): "<<inputMatrix<<std::endl;

  imiSurrogateBasedMotionPredictionHelpers::GenerateVectorFieldFromVnlVector(inputMatrix,outputField,refFieldFilename,maskFilename);

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
      imiERROR(" Save failed with exception:" << excep << std::endl);
      return -1;
    }




  imiINFO( "imiGenerateVectorFieldFromVector FINISHED." );
  imiINFO( "==========================================\n" );
  return 0;
}

