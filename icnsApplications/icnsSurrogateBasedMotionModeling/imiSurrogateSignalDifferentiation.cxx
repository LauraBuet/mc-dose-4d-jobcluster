/** \file imiSurrogateSignalDifferentiation.cxx
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
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

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiObject.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiSurrogateSignalDifferentiation ... \n\n";

  std::cout << "-I            List of input signal observations (.mat) Order matters! | \n";
  std::cout << "-O            List of regressor observations (.mat|.mha)\n";
  std::cout << "-c            Cyclic differentiation (e.g. for 4D-CT data sets)\n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{
  imiINFO( "==================================================" );
  imiINFO( "====    imiSurrogateSignalDifferentiation     ====" );
  imiINFO( "==================================================" );
  imiINFO( "Reading parameters ...\n" );

  if( argc < 2 )
  {
    PrintHelp();
    return 1;
  }

  imiObject::SetGlobalDebugLevel( 8 );

  imiTime::imiInitTime();

  // VARIABLES AND CALL PARAMS:

  std::vector<std::string> inputFilenames;
  std::vector<std::string> outputFilenames;

  bool cyclicFlag = false;

  // Running through cmd line params:

  char c;
  int cnt = 0;

  // Reading parameters
  while( (c = getopt( argc, argv, "O:I:c?" )) != -1 )
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
          std::cout << "  input [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        //reads the output filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputFilenames.push_back( argv[optind] );
          std::cout << "  output [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'c':
        cyclicFlag = true;
        std::cout << "     Cyclic differentiation: ON" << std::endl;
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

  //Read input signal
  std::vector<VnlMatrixType> inputSignal;
  std::vector<VnlMatrixType> outputSignal;

  VnlMatrixType tempInputMatrix;
  for( unsigned int i = 0; i < inputFilenames.size(); i++ )
  {
    imiINFO( "  Reading signal " << inputFilenames[i] << " ..." );
    vcl_ifstream fid;
    fid.open( inputFilenames[i].c_str() );
    vnl_matlab_read_or_die( fid, tempInputMatrix );
    inputSignal.push_back( tempInputMatrix );
  }

  VnlMatrixType currVector;

  currVector.set_size( inputSignal[0].rows() * 2, 1 );

  if( cyclicFlag )
  {
    //copy old content
    for( unsigned int i = 0; i < inputSignal[0].rows(); i++ )
    {
      currVector.put( i, 0, inputSignal[0].get( i, 0 ) );
    }
    //differentiation
    for( unsigned int i = inputSignal[0].rows(); i < currVector.rows(); i++ )
    {
      currVector.put( i, 0, inputSignal[0].get( i - inputSignal[0].rows(), 0 ) - inputSignal[inputSignal.size() - 1].get( i - inputSignal[0].rows(), 0 ) );
    }

    VnlMatrixType firstVector = inputSignal[0] - inputSignal[inputSignal.size() - 1];
  }
  else
  {
    //copy old content
    for( unsigned int i = 0; i < inputSignal[0].rows(); i++ )
    {
      currVector.put( i, 0, inputSignal[0].get( i, 0 ) );
    }
    //differentiation
    for( unsigned int i = inputSignal[0].rows(); i < currVector.rows(); i++ )
    {
      currVector.put( i, 0, 0 );
    }

  }

  outputSignal.push_back( currVector );

  for( unsigned int signal = 1; signal < inputSignal.size(); signal++ )
  {
    //copy old content
    for( unsigned int i = 0; i < inputSignal[0].rows(); i++ )
    {
      currVector.put( i, 0, inputSignal[signal].get( i, 0 ) );
    }
    //differentiation
    for( unsigned int i = inputSignal[0].rows(); i < currVector.rows(); i++ )
    {
      currVector.put( i, 0, inputSignal[signal].get( i - inputSignal[0].rows(), 0 ) - inputSignal[signal - 1].get( i - inputSignal[0].rows(), 0 ) );
    }

    outputSignal.push_back( currVector );

  }

  for( int unsigned index = 0; index < outputFilenames.size(); index++ )
  {
    vnl_matlab_filewrite outVectorFile( outputFilenames[index].c_str() );
    VnlMatrixType writeMatrix = outputSignal[index];
    outVectorFile.write( writeMatrix, "vector" );
  }

  imiINFO( "\n------------------------------------------" );
  imiINFO( "imiSurrogateSignalDifferentiation finished." );
  imiINFO( "===========================================\n" );

  return 0;
} // end of main

