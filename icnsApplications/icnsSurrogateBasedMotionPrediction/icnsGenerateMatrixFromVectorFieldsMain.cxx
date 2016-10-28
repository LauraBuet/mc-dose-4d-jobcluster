/** \file icnsGenerateMatrixFromVectorFieldsMain.cxx
 *
 *  Original author:
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2014 Institute of Medical Informatics,
 *     University of Luebeck
 *
 *  Modified by Rene Werner (2016, ICNS)
 *
 ****************************************************************************/

// System includes:
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
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
  std::cout << "icnsGenerateMatrixFromVectorFields ... " << std::cout;
  std::cout << "Generates a MATLAB matrix (.mat) from a set of vector fields." << std::cout;
  std::cout << "Each (masked) vector field is serialized (x-comp, y-comp, z-comp, x-comp,...)" << std::cout;
  std::cout << "and stored as a column of the output matrix." << std::cout;
  std::cout << std::endl;
  std::cout << "-I            Filenames of the input vector fields." << std::endl;
  std::cout << "-O            Filename of the output Matlab matrix." << std::endl;
  std::cout << "-M            Filename of the mask." << std::endl;
  std::cout << "-h            Print this help." << std::endl;
  std::cout << std::endl;
}

int main( int argc, char *argv[] )
{

  // Initialize time for log-output

  std::cout << "============================================" << std::endl;
  std::cout << "===  imiGenerateMatrixFromVectorFields  ====" << std::endl;
  std::cout << "============================================" << std::endl;
  std::cout << "Reading parameters ...\n" << std::endl;

  std::string outputFilename;
  std::string maskFilename;
  std::vector<std::string> inputFilenames;

  char c;
  int cnt = 0;

  // Reading parameters:
  
  while( (c = getopt( argc, argv, "O:M:I:h?" )) != -1 )
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
          inputFilenames.push_back( argv[optind] );
          std::cout << "  vector field [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        optind--;
        outputFilename = argv[optind];
        std::cout << "  output matrix = " << outputFilename << std::endl;
        break;
      case 'M':
        optind--;
        maskFilename = argv[optind];
        std::cout << "  mask = " << maskFilename << std::endl;
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_SUCCESS;
      default:
        std::cerr <<  "Argument "<<(char)c<<" not processed !\n" << std::endl;
        return EXIT_FAILURE;
    }
  }

  // Loading input data:
  
  VnlMatrixType outputMatrix;
  icnsSurrogateBasedMotionPredictionHelpers::GenerateMatrixFromVectorFields( outputMatrix, inputFilenames, maskFilename );

  std::cout << outputMatrix.cols() << " " << outputMatrix.rows() << std::endl;
  
  vnl_matlab_filewrite outFileWriter( outputFilename.c_str() );
  outFileWriter.write( outputMatrix, "fields" );

  std::cout << "imiGenerateMatrixFromVectorFields FINISHED." << std::endl;
  std::cout << "===========================================" << std::endl;
  std::cout << std::endl;
  
  return EXIT_SUCCESS;
}

