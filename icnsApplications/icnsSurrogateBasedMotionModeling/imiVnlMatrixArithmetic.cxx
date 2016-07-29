/** \file imiVnlMatrixArithmetic.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// System includes:
#include <string>
#include <iostream>
extern "C"
{
#include "getopt.h"
}

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"




using namespace imi;

enum TypeOfProcess
{
  operationAdd,
  operationSub,
  operationMul,
  operationDiv,
  operationMin,
  operationMax,
  operationAnd,
  operationOr
};

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage (imiVnlVectorArithmetic):\n";
  std::cout << "imiVnlVectorArithmetic -F <first matlab file> -S <second matlab file> -O <output matlab file> [-p <arithmetic process>]\n\n";
  std::cout << "-F <image>                 Filename of the first matrix.\n";
  std::cout << "-S <image>                 Filename of the second matrix.\n";
  std::cout << "-O <image>                 Filename of the output matrix.\n";
  std::cout << "-p <process>               Tag for the arithmetric operation:\n";
  std::cout << "                             'a' for Addition (NYI)\n";
  std::cout << "                             's' for Subtraction\n";
  std::cout << "                             'm' for Multiplication (NYI)\n";
  std::cout << "                             'd' for Division (NYI)\n";
  std::cout << "-h                         Print this help.\n";
  std::cout << "\n";
}


int main( int argc, char *argv[] )
{
  if( argc < 5 )
  {
    PrintHelp();
    return 1;
  }

  char* firstFilename = NULL;
  char* secondFilename = NULL;
  char* outputFilename = NULL;
  std::string process = "";
  std::string processName = "addition";
  int c;
  int operation = operationAdd;

  imiINFO( "==============================================" );
  imiINFO( "======      imiVnlVectorArithmetic      ======" );
  imiINFO( "==============================================" );
  imiINFO( "READING parameters...\n" );

  while( (c = getopt( argc, argv, "F:S:O:p:h" )) != -1 )
  {
    switch( c )
    {
      case 'F':
        firstFilename = optarg;
        std::cout << " First filename:                " << firstFilename << std::endl;
        break;
      case 'S':
        secondFilename = optarg;
        std::cout << " Second filename:               " << secondFilename << std::endl;
        break;
      case 'O':
        outputFilename = optarg;
        std::cout << " Output filename:               " << outputFilename << std::endl;
        break;
      case 'p':
        process = optarg;
        if( process == "a" )
        {
          std::cout << " Process:                             Add" << std::endl;
          operation = operationAdd;
          processName = "addition";
        }
        else if( process == "s" )
        {
          std::cout << " Process:                             Substract" << std::endl;
          operation = operationSub;
          processName = "subtraction";
        }
        else if( process == "m" )
        {
          std::cout << " Process:                             Multiply" << std::endl;
          operation = operationMul;
          processName = "multiplication";
        }
        else if( process == "d" )
        {
          std::cout << " Process:                             Divide" << std::endl;
          operation = operationDiv;
          processName = "division";
        }
        else if( process == "i" )
        {
          std::cout << " Process:                             Min" << std::endl;
          operation = operationMin;
          processName = "minimum";
        }
        else if( process == "x" )
        {
          std::cout << " Process:                             Max" << std::endl;
          operation = operationMax;
          processName = "maximum";
        }
        else if( process == "n" )
        {
          std::cout << " Process:                             And" << std::endl;
          operation = operationAnd;
          processName = "and";
        }
        else if( process == "o" )
        {
          std::cout << " Process:                             Or" << std::endl;
          operation = operationOr;
          processName = "or";
        }
        else
        {
          imiERROR( "Unsupported arithmetic process!" );
          return -1;
        }
        imiINFO( "  arithmetic process:                 " << processName );
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


  // check valid arguments.
  if( firstFilename == NULL )
  {
    imiERROR( "No input 1!" );
    return -1;
  }
  if( secondFilename == NULL )
  {
    imiERROR( "No input 2!" );
    return -1;
  }
  if( outputFilename == NULL )
  {
    imiERROR( "No output!" );
    return -1;
  }

  //
  // Load input   //
  imiINFO( "Loading input 1..." );
  VnlMatrixType inputMatrix1;

  vcl_ifstream fid1;
  fid1.open( firstFilename );
  vnl_matlab_read_or_die( fid1, inputMatrix1 );

  imiINFO( "Loading input 2..." );
  VnlMatrixType inputMatrix2;

  vcl_ifstream fid2;
  fid2.open( secondFilename );
  vnl_matlab_read_or_die( fid2, inputMatrix2 );

  VnlMatrixType outputMatrix;


  if (operation == operationSub)
  {
    outputMatrix=inputMatrix1-inputMatrix2;
  }

  imiINFO( "Writing output..." );
  vnl_matlab_filewrite outFile( outputFilename );
  outFile.write( outputMatrix, "matrix" );


  imiINFO( "Finished. " );
  imiINFO( "==========================================\n" );
  return 0;
}
