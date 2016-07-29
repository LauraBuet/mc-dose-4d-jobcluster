/** \file imiGenerateDiaphragmSurrogateMain.cxx
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

// ITK includes:
#include "imiImageReader.h"
#include "imiImageWriter.h"

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiDiaphragmSurrogate.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenerateDiaphragmSurrogate ... \n\n";
  std::cout << "-I            Filenames of the input image data.\n";
  std::cout << "-M            Filenames of the mask files.\n";
  std::cout << "-O            Filename of the output Matlab matrix. \n";
  std::cout << "-V            Filenames of the output Matlab 'vectors'. (one file per breathing phase) \n";
  std::cout << "-s            Coronal slice (y-coord) for diaphragm tracking. \n";
  std::cout << "-p            Starting point (x,z) for diaphragm tracking. (multiple starting points are allowed!) \n";
  std::cout << "-n            Number of diaphragm position samplings (in x direction). \n";
  std::cout << "-f            Head first data (z=0 is anatomically located above the dome of the diaphragm). \n";
  std::cout << "-t            Lung tissue threshold (default=650)\n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{

  // Initialize time for log-output
  imiTime::imiInitTime();


  imiINFO( "=========================================" );
  imiINFO( "===  imiGenerateDiaphragmSurrogate  ====" );
  imiINFO( "=========================================" );
  imiINFO( "Reading parameters ...\n" );

  if( argc < 4 )
  {
    PrintHelp();
    return 1;
  }

  std::vector<char*> outputVectorFilenames;
  char* outputWholeMatrixFilename = NULL;
  std::vector<char*> inputImageFilenames;
  std::vector<char*> inputMaskFilenames;

  unsigned int startSlice=0;
  std::vector<unsigned int> startX;
  std::vector <unsigned int> startZ;
  unsigned int numberOfSamplings=0;
  unsigned int threshold=650;
  bool headFirst=false;

  char c;
  int cnt = 0;

  // Reading parameters
  while( (c = getopt( argc, argv, "O:I:V:M:p:n:s:t:h?df" )) != -1 )
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
          std::cout << "  input image [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'M':
        //reads the mask filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputMaskFilenames.push_back( argv[optind] );
          std::cout << "  input mask [" << cnt << "] = " << argv[optind] << std::endl;
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
      case 's':
        optind--;
        std::cout << std::endl;
        startSlice = atoi( argv[optind] );
        std::cout << "  Starting slice = " << startSlice << std::endl;
        break;
      case 'n':
        optind--;
        std::cout << std::endl;
        numberOfSamplings = atoi( argv[optind] );
        std::cout << "  Number of samplings = " << numberOfSamplings << std::endl;
        break;
      case 'p':
        optind--;
        std::cout << std::endl;
        startX.push_back(atoi( argv[optind] ));
        optind++;
        startZ.push_back(atoi( argv[optind] ));
        std::cout << "  Starting point = ("<< startX[startX.size()-1]<<","<<startZ[startZ.size()-1]<<")"<< std::endl;
        break;
      case 't':
        optind--;
        std::cout << std::endl;
        threshold = atoi( argv[optind] );
        std::cout << "  Lung tissue threshold = " << threshold << std::endl;
        break;
      case 'f':
        std::cout << std::endl;
         headFirst = true;
        std::cout << "  Head first data " << std::endl;
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
  std::vector<ImageType::Pointer> inputImageVector;
  std::vector<ImageType::Pointer> inputMaskVector;

  if( inputImageFilenames.size() > 0 )
  {
    std::cout << std::endl;
    std::cout << "  :::Loading image data (*.dcm)..." << std::endl;
    int inputCount = 0;
    int IMGcount = inputImageFilenames.size();
    imiImageReader *inputImageReader = imiImageReader::New();
    ImageType::Pointer inputImage;
    for( int z = 0; z < IMGcount; z++ )
    {
      if( !inputImageReader->LoadImageData( inputImageFilenames[z] ) )
      {
        imiERROR( "Cannot load input image(s)!" );
        return -1;
      }
      inputImage = inputImageReader->GetITKImage();
      inputImageVector.push_back( inputImage );
      inputCount++;
    }
  }


  if( inputMaskFilenames.size() > 0 )
  {
    std::cout << std::endl;
    std::cout << "  :::Loading mask data (*.dcm)..." << std::endl;
    int inputCount = 0;
    int IMGcount = inputMaskFilenames.size();
    imiImageReader *inputImageReader = imiImageReader::New();
    ImageType::Pointer inputImage;
    for( int z = 0; z < IMGcount; z++ )
    {
      if( !inputImageReader->LoadImageData( inputImageFilenames[z] ) )
      {
        imiERROR( "Cannot load mask image(s)!" );
        return -1;
      }
      inputImage = inputImageReader->GetITKImage();
      inputMaskVector.push_back( inputImage );
      inputCount++;
    }
  }

  VnlMatrixType observationMatrix;
  VnlMatrixType regressorMatrix;

  //computing the regressor matrix Z
  imiDiaphragmSurrogate* diaphragm = imiDiaphragmSurrogate::New();
  diaphragm->SetMaskImages( inputMaskVector );
  diaphragm->SetInputImages( inputImageVector );
  diaphragm->SetCoronalSlicePosition(startSlice);
  diaphragm->SetStartingPoints(startX,startZ);
  diaphragm->SetXSize(numberOfSamplings);
  diaphragm->SetLungThreshold(threshold);
  diaphragm->SetHeadFirst(headFirst);
  diaphragm->ComputeMeasurements();
  diaphragm->GetMeasurementMatrix( regressorMatrix );


    //without differentiation
    VnlMatrixType diaphragmMatrix = VnlMatrixType( 1, regressorMatrix.cols() );
    diaphragmMatrix.set_row( 0, regressorMatrix.get_row( 0 ) );

    vnl_matlab_filewrite outSystemMatrixFileX( outputWholeMatrixFilename );
    outSystemMatrixFileX.write( diaphragmMatrix, "diaphragmMatrix" );

    for( int unsigned index = 0; index < outputVectorFilenames.size(); index++ )
    {
      vnl_matlab_filewrite outBeltVectorFile( outputVectorFilenames[index] );
      VnlMatrixType writeMatrix = VnlMatrixType( 1, 1 );
      writeMatrix.set_column( 0, diaphragmMatrix.get_column( index ) );
      outBeltVectorFile.write( writeMatrix, "diaphragmVector" );
    }

  imiINFO( "imiGenerateDiaphragmSurrogate FINISHED." );
  imiINFO( "==========================================\n" );
  return 0;
}

