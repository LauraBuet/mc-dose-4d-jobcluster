/** \file imiGenerateAirVolumeSurrogate.cxx
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
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiImageReader.h"
#include "imiImageWriter.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenereateAirVolumeSurrogate ... \n\n";
  std::cout << "-I            Filenames of the input gray images.\n";
  std::cout << "-M            Filenames of the mask images.\n";
  std::cout << "-O            Filename of the output Matlab matrix. \n";
  std::cout << "-V            Filenames of the output Matlab 'vectors'. (one file per breathing phase) \n";
  std::cout << "-o            CT offset value. \n";
  std::cout << "-a            Air HU value. \n";
  std::cout << "-t            Soft tissue HU value \n";
  std::cout << "-h            Print this help.\n";
}

bool CalculateAirVolume(
    ImageShort3DType::Pointer inputSegmImage,
    ImageShort3DType::Pointer inputCTImage,
    int thresh,
    double SoftTissueHounsfield,
    double AirHounsfield,
    double &airVolume,
    int &voxelWithAirCount,
    int &voxelCount )
{
  ////////////////////////////////////////////////////////////////
  //
  // setup warper
  //

  typedef itk::ImageRegionConstIterator< ImageShort3DType > ConstIteratorType;
  ConstIteratorType inputSegmIt(   inputSegmImage, inputSegmImage->GetRequestedRegion()  );
  ConstIteratorType inputCTIt(   inputCTImage, inputCTImage->GetRequestedRegion()  );


  ImageShort3DType::PixelType pixelValue;
  const double DifferenceTissueAndAirHounsfield = SoftTissueHounsfield - AirHounsfield;

  /*  unsigned int voxelCount = 0;
    unsigned int voxelWithAirCount = 0;
    double airVolume = 0;
  */
  voxelCount = 0;
  voxelWithAirCount = 0;
  double airContent = 0;
  airVolume = 0;

  for ( inputSegmIt.GoToBegin(), inputCTIt.GoToBegin(); !inputSegmIt.IsAtEnd(); ++inputSegmIt, ++inputCTIt)
  {
    if(inputSegmIt.Get() > thresh)
    {
      voxelCount++;
      pixelValue = inputCTIt.Get();
      if(pixelValue < SoftTissueHounsfield)
      {
        airContent = 1.0 - ((pixelValue - AirHounsfield)/DifferenceTissueAndAirHounsfield);
        airVolume += airContent;
        voxelWithAirCount++;
      }

    }
  }

  ImageShort3DType::SpacingType spacing;
  spacing = inputCTImage->GetSpacing();

  double voxelsize=spacing[0]*spacing[1]*spacing[2];
  airVolume *= voxelsize;

  return true;
}

double CalculateAirVolume(
    ImageShort3DType::Pointer inputSegmImage,
    ImageShort3DType::Pointer inputCTImage,
    int thresh,
    double SoftTissueHounsfield,
    double AirHounsfield )
{
//  imiINFO("Calculate Air volume ... ");

  ImageShort3DType::SpacingType spacing = inputCTImage->GetSpacing();

  int voxelCount = 0;
  int voxelWithAirCount = 0;
  double airVolume = 0;

  if(!CalculateAirVolume( inputSegmImage, inputCTImage, thresh, SoftTissueHounsfield, AirHounsfield,
                          airVolume, voxelWithAirCount, voxelCount))
  {
    imiERROR("Calculating Air volume failed!");
    return -1;
  }

//  imiINFO("Air volume: "<<airVolume<<" , segmented voxel volume: "<<voxelVolume<<" , number of voxels: "<<voxelCount<<" , number of voxels with air: "<<voxelWithAirCount);

  return airVolume;
}






int main( int argc, char *argv[] )
{

  // Initialize time for log-output
  imiTime::imiInitTime();

  imiINFO( "=========================================" );
  imiINFO( "===  imiGenerateAirVolumeSurrogate  ====" );
  imiINFO( "=========================================" );
  imiINFO( "Reading parameters ...\n" );

  std::vector<char*> outputVectorFilenames;
  char* outputWholeMatrixFilename = NULL;
  std::vector<char*> inputSegmentationFilenames;
  std::vector<char*> inputImageFilenames;

  int CTOffsetValue = 1024;
  int AirHUValue = -1000;
  int SoftTissueHUValue = 40;

  char c;
  int cnt = 0;

  // Reading parameters
  while( (c = getopt( argc, argv, "O:M:I:V:o:a:t:?" )) != -1 )
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
          std::cout << "  image [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'M':
        //reads the input filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputSegmentationFilenames.push_back( argv[optind] );
          std::cout << "  mask [" << cnt << "] = " << argv[optind] << std::endl;
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
      case 'o':
        optind--;
        CTOffsetValue = atoi(optarg);
        std::cout << "  CT Offset Value:        " << CTOffsetValue << std::endl;
        break;
      case 'a':
        optind--;
        AirHUValue = atoi(optarg);
        std::cout << "     Air HU Value:        " << AirHUValue << std::endl;
        break;
      case 't':
        optind--;
        SoftTissueHUValue = atoi(optarg);
        std::cout << " Soft Tissue HU Value:    " << SoftTissueHUValue << std::endl;
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

  std::vector<ImagePointerType> inputImageVector;

  if( inputSegmentationFilenames.size() > 0 )
  {
    std::cout << std::endl;
    std::cout << "  :::Loading image data (*.dcm)..." << std::endl;
    int IMGcount = inputImageFilenames.size();
    imiImageReader *inputImageReader = imiImageReader::New();
    ImagePointerType inputImage;
    for( int z = 0; z < IMGcount; z++ )
    {
      if( !inputImageReader->LoadImageData( inputImageFilenames[z] ) )
      {
        imiERROR( "Cannot load input image(s)!" );
        return -1;
      }
      inputImage = inputImageReader->GetITKImage();
      inputImageVector.push_back( inputImage );
    }
  }


  std::vector<ImagePointerType> inputSegmentationVector;

  if( inputSegmentationFilenames.size() > 0 )
  {
    std::cout << std::endl;
    std::cout << "  :::Loading segmentation data (*.dcm)..." << std::endl;
    int SEGcount = inputSegmentationFilenames.size();
    imiImageReader *inputImageReader = imiImageReader::New();
    ImagePointerType inputImage;
    for( int z = 0; z < SEGcount; z++ )
    {
      if( !inputImageReader->LoadImageData( inputSegmentationFilenames[z] ) )
      {
        imiERROR( "Cannot load input image(s)!" );
        return -1;
      }
      inputImage = inputImageReader->GetITKImage();
      inputSegmentationVector.push_back( inputImage );
    }
  }

  VnlMatrixType observationMatrix;

  observationMatrix.set_size(1,inputImageVector.size());


  for (unsigned int i=0; i<inputImageVector.size();i++)
  {
    observationMatrix.put(0,i,static_cast<MatrixValueType>(CalculateAirVolume(inputSegmentationVector[i],inputImageVector[i],0,SoftTissueHUValue + CTOffsetValue,AirHUValue + CTOffsetValue)));
  }

  /*std::cout<<"inputImageVector: "<<inputImageVector.size()<<std::endl;
  std::cout<<"observationMatrix: "<<observationMatrix.cols()<<std::endl;*/


    vnl_matlab_filewrite outSystemMatrixFileX( outputWholeMatrixFilename );
    outSystemMatrixFileX.write( observationMatrix, "airVolumeMatrix" );


    for( int unsigned index = 0; index < outputVectorFilenames.size(); index++ )
    {
      vnl_matlab_filewrite outVectorFile( outputVectorFilenames[index] );
      VnlMatrixType writeMatrix = VnlMatrixType( 1, 1 );
      writeMatrix.set_column( 0, observationMatrix.get_column( index ) );
      outVectorFile.write( writeMatrix, "airVolumeVector" );
    }

  imiINFO( "imiGenerateAirVolumeSurrogate FINISHED." );
  imiINFO( "==========================================\n" );
  return 0;
}

