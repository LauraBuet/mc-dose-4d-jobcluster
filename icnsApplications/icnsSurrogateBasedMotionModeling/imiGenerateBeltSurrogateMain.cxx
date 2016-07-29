/** \file imiGenerateBeltSurrogateMain.cxx
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
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
#include "imiBeltSurrogate.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenereateBeltSurrogate ... \n\n";
  std::cout << "-I            Filenames of the input body segmentation data.\n";
  std::cout << "-O            Filename of the output Matlab matrix. \n";
  std::cout << "-V            Filenames of the output Matlab 'vectors'. (one file per breathing phase) \n";
  std::cout << "-s            Starting slice for belt simulation. \n";
  std::cout << "-n            Number of slices used for belt simulation. \n";
  std::cout << "-d            Save surrogate signal and first derivative. (state augmentation using finite differences) \n";
  std::cout << "-m            'Normalize' first derivative. (1: increasing, -1: decreasing signal) \n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{

  // Initialize time for log-output
  imiTime::imiInitTime();

  imiINFO( "=========================================" );
  imiINFO( "===  imiGenerateBeltSurrogate  ====" );
  imiINFO( "=========================================" );
  imiINFO( "Reading parameters ...\n" );

  std::vector<char*> outputVectorFilenames;
  char* outputWholeMatrixFilename = NULL;
  std::vector<char*> inputSegmentationFilenames;

  //chest belt to measure the patients breathing volume for a specific slice
  int beltLine = 0;
  int beltSize = 0;

  char c;
  int cnt = 0;

  unsigned int roiXLow = 0;
  unsigned int roiXHigh = 0;
  unsigned int roiYLow = 0;
  unsigned int roiYHigh = 0;

  // Reading parameters
  while( (c = getopt( argc, argv, "O:I:V:r::::n:s:h?" )) != -1 )
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
          inputSegmentationFilenames.push_back( argv[optind] );
          std::cout << "  body segmentation [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        optind--;
        outputWholeMatrixFilename = argv[optind];
        std::cout << "  output matrix = " << outputWholeMatrixFilename << std::endl;
        break;
      case 'V':
        //reads in the output filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputVectorFilenames.push_back( argv[optind] );
          std::cout << "  output vector [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'r':
        roiXLow = atoi( argv[optind] );
        roiXHigh = atoi( argv[++optind] );
        roiYLow = atoi( argv[++optind] );
        roiYHigh = atoi( argv[++optind] );
        optind++;
        std::cout << " ROI:                 x= [" << roiXLow << " , " << roiXHigh << "]" << " z= [" << roiYLow << " , " << roiYHigh << "]" << std::endl;
        break;
      case 's':
        optind--;
        std::cout << std::endl;
        beltLine = atoi( argv[optind] );
        std::cout << "  Starting slice = " << beltLine << std::endl;
        break;
      case 'n':
        optind--;
        std::cout << std::endl;
        beltSize = atoi( argv[optind] );
        std::cout << "  Number of slices = " << beltSize << std::endl;
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

  //variables
  std::vector<BinarySegmentationPointerType> inputSegmentationVector;

  if( inputSegmentationFilenames.size() > 0 )
  {
    std::cout << std::endl;
    std::cout << "  :::Loading segmentation data (*.dcm)..." << std::endl;
    int inputCount = 0;
    int SEGcount = inputSegmentationFilenames.size();
    imiImageReader *inputImageReader = imiImageReader::New();
    BinarySegmentationPointerType inputImage;
    for( int z = 0; z < SEGcount; z++ )
    {
      if( !inputImageReader->LoadImageData( inputSegmentationFilenames[z] ) )
      {
        imiERROR( "Cannot load input image(s)!" );
        return -1;
      }
      inputImage = inputImageReader->GetITKBinarySegmentation();
      inputSegmentationVector.push_back( inputImage );
      inputCount++;
    }
  }

  roiXHigh=inputSegmentationVector[0]->GetLargestPossibleRegion().GetSize(0)-1;
  roiYHigh=inputSegmentationVector[0]->GetLargestPossibleRegion().GetSize(1)-1;

  VnlMatrixType observationMatrix;
  VnlMatrixType regressorMatrix;

  //computing the regressor matrix Z
  imiBeltSurrogate* belt = imiBeltSurrogate::New();
  belt->SetInputImages( inputSegmentationVector );
  belt->SetBeltPosition( beltLine );
  belt->SetBeltSize( beltSize );
  belt->SetROI(roiXLow,roiXHigh,roiYLow,roiYHigh);
  belt->ComputeMeasurements();
  belt->GetMeasurementMatrix( regressorMatrix );

  VnlMatrixType pureBelt = VnlMatrixType( 1, regressorMatrix.cols() );
  pureBelt.set_row( 0, regressorMatrix.get_row( 0 ) );

  std::string fileEnding = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( outputWholeMatrixFilename );
  imiINFO("fileEnding: "<<fileEnding<<std::endl);

  if( fileEnding == "mat" )
  {
    imiINFO("Writing matlab file...\n");
    vnl_matlab_filewrite outSystemMatrixFileX( outputWholeMatrixFilename );
    outSystemMatrixFileX.write( pureBelt, "beltMatrix" );
  }
  else if( fileEnding == "txt" )
  {
    imiINFO("Writing txt file...\n");
    std::ofstream file;
    file.open( outputWholeMatrixFilename, std::ios::out );

    for( unsigned int index = 0; index < pureBelt.cols(); index++ )
    {
      file << pureBelt.get( 0, index ) << "\t";
    }
    file << std::endl;
    file.close();
  }
  else
  {
    imiERROR( "Unsupported matrix file ending "<<fileEnding<<"!" );
    return -1;
  }

  for( int unsigned index = 0; index < outputVectorFilenames.size(); index++ )
  {
    vnl_matlab_filewrite outBeltVectorFile( outputVectorFilenames[index] );
    VnlMatrixType writeMatrix = VnlMatrixType( 1, 1 );
    writeMatrix.set_column( 0, pureBelt.get_column( index ) );
    outBeltVectorFile.write( writeMatrix, "beltVector" );
  }

  imiINFO( "imiGenerateBeltSurrogate FINISHED." );
  imiINFO( "==========================================\n" );
  return 0;
}

