/** \file imiMATLABMatrixCompare.cxx
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#include <string>
#include <fstream>
extern "C"
{
#include "getopt.h"
}

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace std;
using namespace imi;

/** print help function. */

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Compares two MATLAB matrices (.mat files) . Return value is 0 (SUCCESS), if matrices are equal, and -1 otherwise.\n";
  std::cout << "Usage :\n";
  std::cout << "imiMATLABMatrixCompare -F <input MATLAB file 1> -S <input MATLAB file 2> [-t <thresh> -d <0|1>]\n\n";
  std::cout << "-F <input MATLAB file 1>    Filename of the first input vector field.\n";
  std::cout << "-S <input MATLAB file 2>    Filename of the second input vector field.\n";
  std::cout << "-t <thresh>                 Threshold for matrix elements to be different (default=0).\n";
  std::cout << "-d <0|1>                    Data type of the matrix elements (0=float=default,1=double).\n";
  std::cout << "-h                          Print this help.\n";
  std::cout << "\n";
}

template<class TValueType> bool MatrixCompare( string inputMatrix1Filename, string inputMatrix2Filename, double diffThreshold)
{

  typedef vnl_matrix<TValueType> MatrixType;

  MatrixType matrix1;
  MatrixType matrix2;


  ////////////////////////////////////////////////////////////////
  //
  // load matrices
  //


  string fileTypePattern1, fileTypePattern2;


  fileTypePattern1 = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( inputMatrix1Filename );
  fileTypePattern2 = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( inputMatrix2Filename );

  if( strcmp( fileTypePattern1.c_str(), "mat" ) == 0 )
  {
    vcl_ifstream fid1;
    fid1.open( inputMatrix1Filename.c_str() );
    if(!vnl_matlab_read_or_die( fid1, matrix1 ))
    {
      imiERROR( "Could not read file: "<< inputMatrix1Filename);
      return false;
    }
    fid1.close();
  }
  else
  {
    imiERROR( fileTypePattern1 << " files are not supported! Use .mat files instead." );
    return false;
  }

  if( strcmp( fileTypePattern2.c_str(), "mat" ) == 0 )
  {
    vcl_ifstream fid2;
    fid2.open( inputMatrix2Filename.c_str() );
    if(!vnl_matlab_read_or_die( fid2, matrix2 ))
    {
      imiERROR( "Could not read file: "<< inputMatrix2Filename);
      return false;
    }
    fid2.close();
  }
  else
  {
    imiERROR( fileTypePattern2 << " files are not supported! Use .mat files instead." );
    return false;
  }

  ////////////////////////////////////////////////////////////////
  //
  // check sizes
  //
  if (matrix1.rows()!=matrix2.rows() || matrix1.cols()!=matrix2.cols() )
  {
    imiERROR("Matrices are of different size!");
    return false;
  }

  ////////////////////////////////////////////////////////////////
  //
  // check element-wise differences
  //

  //imiINFO("matrix1: "<<matrix1.rows()<<" x "<<matrix1.cols());
  //imiINFO("matrix2: "<<matrix2.rows()<<" x "<<matrix2.cols());

  double totalDifference=0;
  unsigned int elemCount=matrix1.rows()*matrix1.cols();
  unsigned int numberOfDiffElems=0;
  for (unsigned int i=0; i< matrix1.rows();i++)
  {
    for (unsigned int j=0; j< matrix1.cols();j++)
    {
      //actually a bad idea to use fabs(a-b) because of catastrophic cancellation etc.
      const double diff = fabs(matrix1.get(i,j)-matrix2.get(i,j));

      //imiINFO("diff: "<<diff);

      if (diff > diffThreshold)
      {
        totalDifference+=diff;
        numberOfDiffElems++;
      }
    }
  }

  const double epsilon = 0.0000001;
  if(fabs(totalDifference) > epsilon)
  {

    imiINFO( "Matrices are different!\n");
    imiINFO( "Mean Difference: " << totalDifference/static_cast<double>(elemCount)<< "\n");
    imiINFO( "Total Difference: " << totalDifference << "\n");
    imiINFO( "Number of elements with differences: " << numberOfDiffElems << "\n");
    return false;
  }
  else
  {
    imiINFO( "Matrices are equal!" );
  }

  return true;
}

int main( int argc, char *argv[] )
{
  if( argc < 4 )
  {
    PrintHelp();
    return 1;
  }

  imiINFO( "========================================" );
  imiINFO( "=====    imiMATLABMatrixCompare    =====" );
  imiINFO( "========================================" );
  imiINFO( "READING parameters...\n" );

  // Initialize parameters with default values.
  int c;

  double diffThreshold = 0.0;
  string inputMatrix1Filename;
  string inputMatrix2Filename;
  bool useDouble = false;
  unsigned int  retVal=0;

  // Reading parameters
  while( (c = getopt( argc, argv, "S:F:t:d:h" )) != -1 )
  {
    switch( c )
    {
      case 'F':
        inputMatrix1Filename = optarg;
        std::cout << " Input matrix 1 filename:         " << inputMatrix1Filename << std::endl;
        break;
      case 'S':
        inputMatrix2Filename = optarg;
        std::cout << " Input matrix 2 filename:         " << inputMatrix2Filename << std::endl;
        break;
      case 't':
        diffThreshold = atof( optarg );
        std::cout << "Difference threshold:                 " << diffThreshold << std::endl;
        break;
      case 'd':
        if( atoi( optarg ) == 0 )
        {
          useDouble = false;
          std::cout << "Float data type is used to load/store the MATLAB matrices!" << std::endl;
        }
        else if( atoi( optarg ) == 1 )
        {
          useDouble = true;
          std::cout << "Double data type is used to load/store the MATLAB matrices!" << std::endl;
        }
        break;
      case 'h':
      case '?':
        PrintHelp();
        exit( 1 );
        break;
      default:
        imiERROR( " Argument "<<(char)c<<" not processed !\n" );
    }
  }

  // Read image data.
  if( inputMatrix1Filename.empty() )
  {
    imiERROR( "No input matrix 1!" );
    return 1;
  }
  if( inputMatrix2Filename.empty() )
  {
    imiERROR( "No input matrix 2!" );
    return 1;
  }

  if(useDouble)
  {
    if(!MatrixCompare<double>( inputMatrix1Filename, inputMatrix2Filename, diffThreshold))
    {
      retVal=1;
    }
  } else
  {
    if(!MatrixCompare<float>( inputMatrix1Filename, inputMatrix2Filename, diffThreshold))
    {
      retVal=1;
    }
  }

  imiINFO( "Finished." );
  imiINFO( "==========================================\n" );
  return retVal;
}
