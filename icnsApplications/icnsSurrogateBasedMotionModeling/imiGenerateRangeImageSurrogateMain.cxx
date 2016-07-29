/** \file imiGenerateRangeImageSurrogateMain.cxx
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
#include "imiRangeImageSurrogate.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenerateRangeImageSurrogate ... \n\n";
  std::cout << "-I            Filenames of the input images.\n";
  std::cout << "-O            Filenames of the output range images (Matlab matrices - one file per breathing phase). \n";
  std::cout << "-V            Filenames of the output vectors (Matlab 'vectors' - one file per breathing phase) \n";
  std::cout
      << "-C            Base filenames (one file per breathing phase; the suffix _col_$COLINDEX.mat will be automatically added) used to generate Matlab 'vectors' for each column of each range image.\n";
  std::cout
      << "-R            Base filenames (one file per breathing phase; the suffix _row_$ROWINDEX.mat will be automatically added) used to generate Matlab 'vectors' for each row of each range image.\n";
  std::cout
      << "-D            Base filenames (one file per breathing phase; the suffix _diag_$DIAGINDEX.mat will be automatically added) used to generate Matlab 'vectors' for both diagonals of each range image (only valid for quadratic range images!).\n";
  std::cout
      << "-P            Base filenames (one file per breathing phase; the suffix _pt_$POINTINDEX.mat will be automatically added) used to generate Matlab 'vectors' for each point of each range image.\n";
  std::cout << "-W            Filename of a single output file containing all vectors (.mat or .txt)\n";
  std::cout << "-r            ROI for surface sampling <x low> <x high> <z low> <z high>. \n";
  std::cout << "-s            Number of sampling points in <x-direction> and <z-direction>. \n";
  std::cout << "-t            Skin surface threshold. (default=500)\n";
  std::cout << "-d            Distance between camera and body in mm (default=1000)\n";
  std::cout << "-k            Kinect noise model (adds noise to the distance values)\n";
  std::cout << "-c            Center sampling points inside the ROI. (adds a margin between first/last point and the ROI border) \n";
  std::cout << "-b            Bilateral filtering of the input images for noise removal/reduction (<sigma spatial> in spacing units and <sigma intensities> in units of intensity) \n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{

  // Initialize time for log-output
  imiTime::imiInitTime();

  imiINFO( "==========================================" );
  imiINFO( "===  imiGenerateRangeImageSurrogate  ====" );
  imiINFO( "==========================================" );
  imiINFO( "Reading parameters ...\n" );

  std::vector<char*> outputVectorFilenames;
  std::vector<char*> outputMatrixFilenames;
  std::vector<char*> inputFilenames;
  std::vector<char*> outputRowVectorFilenames;
  std::vector<char*> outputColVectorFilenames;
  std::vector<char*> outputDiagVectorFilenames;
  std::vector<char*> outputPointVectorFilenames;
  char* outputWholeMatrixFilename = NULL;

  unsigned int samplingX = 0;
  unsigned int samplingZ = 0;
  float roiXLow = 0;
  float roiXHigh = 0;
  float roiZLow = 0;
  float roiZHigh = 0;
  int skinThreshold = 500;
  bool centerSamplingPoints = false;
  float sigmaSpatial = 1.0;
  float sigmaIntens = 200;
  bool useBilateralFilter = false;
  bool useKinectNoise = false;
  float cameraDistance = -1;

  char c;
  int cnt = 0;

  // Reading parameters
  while( (c = getopt( argc, argv, "O:I:V:C:D:P:R:W:r::::s::b::t:d:kch?" )) != -1 )
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
          std::cout << "  Input image [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        //reads in the output filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputMatrixFilenames.push_back( argv[optind] );
          std::cout << "  output matrix [" << cnt << "] = " << argv[optind] << std::endl;
        }
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
      case 'C':
        //reads in the output filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputColVectorFilenames.push_back( argv[optind] );
          std::cout << "  output column vector [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'R':
        //reads in the output filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputRowVectorFilenames.push_back( argv[optind] );
          std::cout << "  output row vector [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'D':
        //reads in the output filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputDiagVectorFilenames.push_back( argv[optind] );
          std::cout << "  output diagonal vector [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'P':
        //reads in the output filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputPointVectorFilenames.push_back( argv[optind] );
          std::cout << "  output point vector [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'W':
        optind--;
        outputWholeMatrixFilename = argv[optind];
        std::cout << "  output matrix (single file) = " << outputWholeMatrixFilename << std::endl;
        break;
      case 'r':
        roiXLow = atof( argv[optind] );
        roiXHigh = atof( argv[++optind] );
        roiZLow = atof( argv[++optind] );
        roiZHigh = atof( argv[++optind] );
        optind++;
        std::cout << " ROI:                 x= [" << roiXLow << " , " << roiXHigh << "]" << " z= [" << roiZLow << " , " << roiZHigh << "]" << std::endl;
        break;
      case 's':
        samplingX = atof( argv[optind] );
        samplingZ = atof( argv[++optind] );
        optind++;
        std::cout << " Sampling:                 x= " << samplingX << " points; z= " << samplingZ << " points" << std::endl;
        break;
      case 'b':
        sigmaSpatial = atof( argv[optind] );
        sigmaIntens = atof( argv[++optind] );
        useBilateralFilter = true;
        optind++;
        std::cout << " Bilateral filtering:                 sigma (spatial)= " << sigmaSpatial << "; sigma (intensities)= " << sigmaIntens << std::endl;
        break;
      case 'd':
        cameraDistance = atof( optarg );
        std::cout << " Camera distance:                 d= " << cameraDistance << " mm" << std::endl;
        break;
      case 't':
        skinThreshold = atoi( optarg );
        std::cout << " Skin surface threshold:       " << skinThreshold << std::endl;
        break;
      case 'k':
        useKinectNoise = true;
        std::cout << " Use Kinect noise model:       true" << std::endl;
        break;
      case 'c':
        centerSamplingPoints = true;
        std::cout << " Center sampling points:       true" << std::endl;
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

  if( inputFilenames.size() > 0 )
  {
    std::cout << std::endl;
    std::cout << "  :::Loading image data ..." << std::endl;
    int inputCount = 0;
    int IMGcount = inputFilenames.size();
    imiImageReader *inputImageReader = imiImageReader::New();
    ImageType::Pointer inputImage;
    for( int z = 0; z < IMGcount; z++ )
    {
      if( !inputImageReader->LoadImageData( inputFilenames[z] ) )
      {
        imiERROR( "Cannot load input image(s)!" );
        return -1;
      }
      inputImage = inputImageReader->GetITKImage();
      inputImageVector.push_back( inputImage );
      inputCount++;
    }
  }

  VnlMatrixType observationMatrix;
  VnlMatrixType regressorMatrix;

  imiRangeImageSurrogate* surrogate = imiRangeImageSurrogate::New();
  surrogate->SetImageListForChestComputing( inputImageVector );

  //surrogate->SetROI(0,inputImageVector[0]->GetLargestPossibleRegion().GetSize(0)-1,0,inputImageVector[0]->GetLargestPossibleRegion().GetSize(2)-1);
  surrogate->SetROI( roiXLow, roiXHigh, roiZLow, roiZHigh );
  surrogate->SetSampling( samplingX, samplingZ );
  surrogate->SetSkinThreshold( skinThreshold );
  surrogate->SetCenterSamplingPoints( centerSamplingPoints );
  surrogate->SetUseBilateralFilter( useBilateralFilter );
  surrogate->SetSigmaIntens( sigmaIntens );
  surrogate->SetSigmaSpatial( sigmaSpatial );
  surrogate->SetAddKinectNoise( useKinectNoise );
  surrogate->SetCameraDistance( cameraDistance );

  surrogate->PerformSurfaceTracking();

  std::vector<VnlMatrixType> rangeImages = surrogate->GetRangeImages();

  for( int unsigned index = 0; index < outputMatrixFilenames.size(); index++ )
  {
    vnl_matlab_filewrite outRangeImageFile( outputMatrixFilenames[index] );
    VnlMatrixType writeMatrix = rangeImages[index];
    outRangeImageFile.write( writeMatrix, "rangeImage" );
  }

  std::vector<VnlMatrixType> rangeVectors;

  for( int unsigned index = 0; index < rangeImages.size(); index++ )
  {
    VnlMatrixType writeMatrix = rangeImages[index];
    //convert matrix to vector
    VnlMatrixType writeVector = VnlMatrixType( writeMatrix.rows() * writeMatrix.cols(), 1 );
    for( unsigned int col = 0; col < writeMatrix.cols(); col++ )
    {
      for( unsigned int row = 0; row < writeMatrix.rows(); row++ )
      {
        writeVector.put( col * writeMatrix.rows() + row, 0, writeMatrix.get( row, col ) );
      }
    }
    rangeVectors.push_back( writeVector );
  }

  for( int unsigned index = 0; index < outputVectorFilenames.size(); index++ )
  {
    vnl_matlab_filewrite outRangeImageFile( outputVectorFilenames[index] );
    outRangeImageFile.write( rangeVectors[index], "rangeImageVector" );
  }

  for( int unsigned index = 0; index < outputColVectorFilenames.size(); index++ )
  {
    for( unsigned int col = 0; col < rangeImages[index].cols(); col++ )
    {
      std::stringstream filename;
      filename << outputColVectorFilenames[index];
      filename << "_col_" << col << ".mat";
      vnl_matlab_filewrite outRangeImageFile( filename.str().c_str() );
      outRangeImageFile.write( rangeImages[index].get_n_columns( col, 1 ), "rangeImageColumn" );
    }
  }

  for( int unsigned index = 0; index < outputRowVectorFilenames.size(); index++ )
  {
    for( unsigned int row = 0; row < rangeImages[index].rows(); row++ )
    {
      std::stringstream filename;
      filename << outputRowVectorFilenames[index];
      filename << "_row_" << row << ".mat";
      vnl_matlab_filewrite outRangeImageFile( filename.str().c_str() );
      outRangeImageFile.write( rangeImages[index].get_n_rows( row, 1 ).transpose(), "rangeImageRow" );
    }
  }

  if( outputDiagVectorFilenames.size() > 0 )
  {
    if( samplingX != samplingZ )
    {
      imiERROR( "Cannot determine diagonals of non-quadratic matrices!" );
    }
    else
    {
      for( int unsigned index = 0; index < outputDiagVectorFilenames.size(); index++ )
      {
        VnlMatrixType diagonal1 = VnlMatrixType( samplingX, 1 );
        VnlMatrixType diagonal2 = VnlMatrixType( samplingX, 1 );

        //primary diagonal top left->bottom right
        for( unsigned int i = 0; i < samplingX; i++ )
        {
          diagonal1.put( i, 0, rangeImages[index].get( i, i ) );
        }

        //secondary diagonal bottom left->top right
        for( unsigned int r = samplingX-1, c = 0; c < samplingX; c++, r-- )
        {
          diagonal2.put( c, 0, rangeImages[index].get( r, c ) );
        }

        std::stringstream filename1;
        filename1 << outputDiagVectorFilenames[index];
        filename1 << "_diag_0.mat";
        vnl_matlab_filewrite outRangeImageFile1( filename1.str().c_str() );
        outRangeImageFile1.write( diagonal1, "rangeImageDiag" );

        std::stringstream filename2;
        filename2 << outputDiagVectorFilenames[index];
        filename2 << "_diag_1.mat";
        vnl_matlab_filewrite outRangeImageFile2( filename2.str().c_str() );
        outRangeImageFile2.write( diagonal2, "rangeImageDiag" );
      }
    }
  }

  //points
  for( int unsigned index = 0; index < outputPointVectorFilenames.size(); index++ )
  {
    for( unsigned int row = 0; row < rangeImages[index].rows(); row++ )
    {
      for( unsigned int col = 0; col < rangeImages[index].cols(); col++ )
      {
        VnlMatrixType point=VnlMatrixType(1,1);
        point.put(0,0,rangeImages[index].get(row,col));

        std::stringstream filename;
        filename << outputPointVectorFilenames[index];
        filename << "_pt_" << col * rangeImages[index].rows() + row << ".mat";
        vnl_matlab_filewrite outRangeImageFile( filename.str().c_str() );
        outRangeImageFile.write( point, "rangeImagePoint" );
      }
    }
  }

  if( outputWholeMatrixFilename != NULL )
  {

    std::string fileEnding = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( outputWholeMatrixFilename );

    if( fileEnding == "mat" )
    {
      VnlMatrixType writeMatrix = VnlMatrixType( rangeVectors[0].rows(), rangeVectors.size() );
      for( unsigned int col = 0; col < rangeVectors.size(); col++ )
      {
        writeMatrix.set_column( col, rangeVectors[col].get_column( 0 ) );
      }

      vnl_matlab_filewrite outRangeImageFile( outputWholeMatrixFilename );
      outRangeImageFile.write( writeMatrix, "rangeImageMatrix" );
    }
    else if( fileEnding == "txt" )
    {
      std::ofstream file;
      file.open( outputWholeMatrixFilename, std::ios::out );

      for( unsigned int row = 0; row < rangeVectors[0].rows(); row++ )
      {
        for( unsigned int index = 0; index < rangeVectors.size(); index++ )
        {
          file << rangeVectors[index].get( row, 0 ) << "\t";
        }
        file << std::endl;
      }
      file.close();
    }
  }

  imiINFO( "imiGenerateRangeImageSurrogate FINISHED." );
  imiINFO( "==========================================\n" );
  return 0;
}
