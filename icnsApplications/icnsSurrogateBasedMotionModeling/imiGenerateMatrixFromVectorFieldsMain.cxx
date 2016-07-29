/** \file imiGenerateMatrixFromVectorFieldsMain.cxx
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
  std::cout << "imiGenerateMatrixFromVectorFields ... \n\n";
  std::cout << "Generates a MATLAB matrix (.mat) from a set of vector fields. \n";
  std::cout << "Each (masked) vector field is serialized (x-comp, y-comp, z-comp, x-comp,...) and stored as a column of the output matrix. \n\n";
  std::cout << "-I            Filenames of the input vector fields.\n";
  std::cout << "-O            Filename of the output Matlab matrix. \n";
  std::cout << "-M            Filename of the mask. \n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{

  // Initialize time for log-output
  imiTime::imiInitTime();

  imiINFO( "============================================" );
  imiINFO( "===  imiGenerateMatrixFromVectorFields  ====" );
  imiINFO( "============================================" );
  imiINFO( "Reading parameters ...\n" );

  std::string outputFilename;
  std::string maskFilename;
  std::vector<std::string> inputFilenames;

  char c;
  int cnt = 0;


  // Reading parameters
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
  VnlMatrixType outputMatrix;

  imiSurrogateBasedMotionPredictionHelpers::GenerateMatrixFromVectorFields(outputMatrix,inputFilenames,maskFilename);

  vnl_matlab_filewrite outFileWriter( outputFilename.c_str() );

  outFileWriter.write(outputMatrix,"fields");

  imiINFO( "imiGenerateMatrixFromVectorFields FINISHED." );
  imiINFO( "===========================================\n" );
  return 0;
}

