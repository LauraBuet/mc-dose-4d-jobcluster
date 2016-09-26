/** \file imiValidateLandmarks.cpp
 *
 *  \b Initial \b Author: Alexander Schmidt-Richberg\n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 *  modified: RW, 2016-09
 *
 ****************************************************************************/

// System includes:
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>

extern "C"
{
#include "getopt.h"
}

// ITK includes
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include <itkCastImageFilter.h>
#include <itkWarpImageFilter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkVectorLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkMultiplyImageFilter.h>
#include <itkVotingBinaryIterativeHoleFillingImageFilter.h>
#include <itkSignedMaurerDistanceMapImageFilter.h>

// Project includes:

// Global typedefs:

typedef itk::Image<short, 3>                         ImageType;
typedef itk::ImageFileReader<ImageType>              ImageReaderType;
typedef itk::ImageFileWriter<ImageType>              ImageWriterType;
typedef itk::Image<short, 3>                         BinarySegmentationType;
typedef itk::ImageFileReader<BinarySegmentationType> BinarySegmentationReaderType;
typedef itk::Image<itk::Vector<float, 3>, 3>         DisplacementFieldType;
typedef itk::ImageFileReader<DisplacementFieldType>  DisplacementFieldReaderType;

typedef DisplacementFieldType::PointType             PointType;
typedef std::vector<PointType>                       PointsType;

bool bVerbose = false;
bool bUseShiftLandmarks = false;
float shiftLandmarks[3] = { 0, 0, 0 };
bool bUseCorrdinateTransformation = false;
float spacing[3] = { 0, 0, 0 };

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage :\n";
  std::cout << "icnsValidateLandmarks -R <reference landmarks> -T <target landmarks> -D <displ field>\n\n";

  std::cout << "-R <reference landmarks>         Filename of the reference landmark file.\n";
  std::cout << "-T <target landmarks>            Filename of the target landmark file.\n";
  std::cout << "-D <displ field>                 Filename of the displacement field.\n";
  std::cout << "-M <mask image>                  Filename of reference mask image (optional).\n";
  std::cout << "-L <logfile>                     Filename of the logfile.\n";
  std::cout << "-O <output filename>             Write output image files with plotted landmarks.\n";
  std::cout << "-s <x-shift> <y-shift> <z-shift> Shift all Landmarks with given values (world coordinates).\n";
  std::cout << "-w <x-spac.> <y-spac.> <z-spac.> Spacing for transformation in world coordinates.\n";
  std::cout << "-f <scale factor>                Scale displacement field with factor.\n";
  std::cout << "-d                               Plot surface distances (mask required).\n";
  std::cout << "-c                               Snap warped LM position to voxel during evaluation.\n";
  std::cout << "-v                               Verbose checking (plot point coordinates and image index).\n";
  std::cout << "-h                               Print this help.\n";
  std::cout << "\n";
}

///////////////////////////////////////////
//
//  ReadLandmarks
//
///////////////////////////////////////////
PointsType ReadLandmarks( const char* filename )
{
  PointsType newPoints;

  // Define filepointer.
  std::ifstream f;
  f.open( filename );

  std::string line;
  std::istringstream iss;
  iss >> std::skipws;

  while( std::getline( f, line ) )
  {
    PointType newPoint;
    iss.clear();
    iss.str( line );
    iss >> newPoint[0] >> newPoint[1] >> newPoint[2];

    if( bUseCorrdinateTransformation )
    {
      newPoint[0] *= spacing[0];
      newPoint[1] *= spacing[1];
      newPoint[2] *= spacing[2];
    }

    if( bUseShiftLandmarks )
    {
      newPoint[0] += shiftLandmarks[0];
      newPoint[1] += shiftLandmarks[1];
      newPoint[2] += shiftLandmarks[2];
    }

    newPoints.push_back( newPoint );
  }

  f.close();

  // Return pointlist.
  return newPoints;
}

///////////////////////////////////////////
//
//  WarpLandmarks
//
///////////////////////////////////////////
PointsType WarpLandmarks( const PointsType points, const DisplacementFieldType::Pointer displField, const bool snapToVoxel = false )
{
  PointsType warpedPoints;

  //  Define interpolator and landmark.
  typedef itk::VectorLinearInterpolateImageFunction<
      DisplacementFieldType,
      DisplacementFieldType::PointType::ValueType> InterpolatorType;

  InterpolatorType::OutputType interpolatedVector;

  // Setup interpolator.
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  interpolator->SetInputImage( displField );

  // Interpolate all landmarks.
  for( unsigned int i = 0; i < points.size(); i++ )
  {
    PointType newPoint;

    // Evaluate interpolated vector.
    interpolatedVector = interpolator->Evaluate( points[i] );

    // Calculate deformed landmark.
    for( unsigned int dim = 0; dim < 3; dim++ )
    {
      newPoint[dim] = points[i][dim] + interpolatedVector[dim];
      if( snapToVoxel )
      {
        InterpolatorType::IndexType discrIndex;
        displField->TransformPhysicalPointToIndex( newPoint, discrIndex );
        displField->TransformIndexToPhysicalPoint( discrIndex, newPoint );
      }
    }
    if( bVerbose )
    {
      InterpolatorType::ContinuousIndexType contIndex;
      displField->TransformPhysicalPointToContinuousIndex( points[i], contIndex );
      std::cout << " " << i << ":\t" << points[i] <<
        "\t transl. " << interpolatedVector << "\t index " << contIndex << std::endl;
    }
    // Store result.
    warpedPoints.push_back( newPoint );
  }

  return warpedPoints;
}

///////////////////////////////////////////
//
//  CheckPointsWithImage
//
///////////////////////////////////////////
bool CheckPointsWithImage( const PointsType points, ImageType::Pointer inputImage )
{
  if( inputImage.IsNull() )
  {
      std::cerr << "Image to check landmarks is NULL!" << std::endl;
    return false;
  }

  bool bIsOK = true;

  if( bVerbose )
  {
    std::cout << " Image Information:" << std::endl;
    std::cout << "    Origin: " << inputImage->GetOrigin() << std::endl;
    std::cout << "   Spacing: " << inputImage->GetSpacing() << std::endl;
    std::cout << "      Size: " << inputImage->GetLargestPossibleRegion().GetSize() << std::endl;
  }
  ImageType::SizeType size = inputImage->GetLargestPossibleRegion().GetSize();
  ImageType::IndexType index;
  bool bIndexOutside = false;
  for( unsigned int i = 0; i < points.size(); i++ )
  {
    inputImage->TransformPhysicalPointToIndex( points[i], index );
    if( bVerbose )
      std::cout << "  " << i << ":\t" << points[i] << "\t --> " << index << std::endl;

    bIndexOutside = false;
    for( int k = 0; k < 3; k++ )
    {
      if( index[k] < 0 || static_cast<unsigned int> ( index[k] ) >= size[k] )
      {
        if( bVerbose )
          std::cout << "Index outside image region!"<< std::endl;
        bIndexOutside = true;
        bIsOK = false;
      }
    }
    if( !bIndexOutside )
    {
      inputImage->SetPixel( index, 255 );
    }
  }

  return bIsOK;
}

///////////////////////////////////////////
//
//  main
//
///////////////////////////////////////////
int main( int argc, char *argv[] )
{
  if( argc < 4 )
  {
    PrintHelp();
    return 1;
  }

  std::cout << "==========================================" << std::endl;
  std::cout << "=====     icnsValidateLandmarks      =====" << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "READING parameters...\n" << std::endl;

  // Initalising parameters with default values
  int c;
  char* logFilename = NULL;
  char* targLandmarksFilename = NULL;
  char* refLandmarksFilename = NULL;
  char* displFieldFilename = NULL;
  char* maskImageFilename = NULL;
  char* outputImageFilename = NULL;
  bool bSurfDist = false;
  bool bMaskedEvaluation = false;
  bool bScaleDisplacementField = false;
  bool bSnapToVoxel = false;
  double displFieldScaleFactor = 1.0;

  // Reading parameters
  while( (c = getopt( argc, argv, "R:T:D:M:L:O:s:::w:::f:dcvh" )) != -1 )
  {
    switch( c )
    {
      case 'T':
        targLandmarksFilename = optarg;
        std::cout << "  Target landmark file:        " << targLandmarksFilename << std::endl;
        break;
      case 'R':
        refLandmarksFilename = optarg;
        std::cout << "  Reference landmark file:     " << refLandmarksFilename << std::endl;
        break;
      case 'D':
        displFieldFilename = optarg;
        std::cout << "  Displ field filename:        " << displFieldFilename << std::endl;
        break;
      case 'M':
        bMaskedEvaluation = true;
        maskImageFilename = optarg;
        std::cout << "  Mask image filename:         " << maskImageFilename << std::endl;
        break;
      case 'L':
        logFilename = optarg;
        std::cout << "  Log filename:                " << logFilename << std::endl;
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << "  Output landmark image:       " << outputImageFilename << std::endl;
        break;
      case 's':
        bUseShiftLandmarks = true;
        shiftLandmarks[0] = atof( argv[optind] );
        shiftLandmarks[1] = atof( argv[++optind] );
        shiftLandmarks[2] = atof( argv[++optind] );
        optind++;
        std::cout << "  Shift landmarks:             [ " << shiftLandmarks[0] << " , " << shiftLandmarks[1] << " , " << shiftLandmarks[2] << " ]" << std::endl;
        break;
      case 'w':
        bUseCorrdinateTransformation = true;
        spacing[0] = atof( argv[optind] );
        spacing[1] = atof( argv[++optind] );
        spacing[2] = atof( argv[++optind] );
        optind++;
        std::cout << "  Transformation with spacing: [ " << spacing[0] << " , " << spacing[1] << " , " << spacing[2] << " ]" << std::endl;
        break;
      case 'f':
        bScaleDisplacementField = true;
        displFieldScaleFactor = atof( optarg );
        std::cout << "  Scale displacements with:    " << displFieldScaleFactor << std::endl;
        break;
      case 'd':
        bSurfDist = true;
        std::cout << "  Plot surface distances:      On" << std::endl;
        break;
      case 'c':
        bSnapToVoxel = true;
        std::cout << "  Snap warped LM to vx center: On" << std::endl;
        break;
      case 'v':
        bVerbose = true;
        std::cout << "  Verbose checking:            On" << std::endl;
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_FAILURE;
      default:
        std::cout << "ERROR: Argument " << (char) c << " not processed !\n" << std::endl;
        break;
    }
  }

  if( targLandmarksFilename == NULL )
  {
    std::cerr << "No target landmark file!" << std::endl;
    return EXIT_FAILURE;
  }
  if( refLandmarksFilename == NULL )
  {
    std::cerr << "No reference landmark file!" << std::endl;
    return EXIT_FAILURE;
  }
  if( displFieldFilename == NULL )
  {
    std::cerr << "No displacement field!" << std::endl;
    return EXIT_FAILURE;
  }

  ////////////////////////////////////////////////////////////////
  //
  // Load landmarks.
  //
  std::cout << "Reading target landmarks... " << std::endl;
  PointsType targetPoints = ReadLandmarks( targLandmarksFilename );
  std::cout << "Reading reference landmarks... " << std::endl;
  PointsType refPoints = ReadLandmarks( refLandmarksFilename );

  if( targetPoints.size() == 0 )
  {
    std::cerr << "No valid target landmarks!" << std::endl;
    return EXIT_FAILURE;
  }
  if( refPoints.size() == 0 )
  {
    std::cerr << "No valid target landmarks!" << std::endl;
    return EXIT_FAILURE;
  }
  if( refPoints.size() != targetPoints.size() )
  {
    std::cout << "RefPointsSize: " << refPoints.size();
    std::cout << "TargetPointsSize: " << targetPoints.size();
    std::cerr << "Different number of reference and target landmarks!" << std::endl;
    return EXIT_FAILURE;
  }

  // Loading displacement field:
    
  std::cout << "Loading displacement field..." << std::endl;
  typedef itk::ImageFileReader<DisplacementFieldType> ImageReaderType;
  ImageReaderType::Pointer imageReader = ImageReaderType::New();

  imageReader->SetFileName( displFieldFilename );
  try
  {
    imageReader->Update();
  }
  catch( itk::ExceptionObject &err )
  {
    std::cerr << "Cannot read displacement field! Exception error: " << err << std::endl;
    return EXIT_FAILURE;
  }

  DisplacementFieldType::Pointer displField = imageReader->GetOutput();

  // Loading mask image:
  
  BinarySegmentationType::Pointer maskImage = NULL;
  if( bMaskedEvaluation )
  {
    std::cout << "  Loading mask image ..." << std::flush;
    BinarySegmentationReaderType::Pointer maskImageReader = BinarySegmentationReaderType::New();
    maskImageReader->SetFileName( maskImageFilename );
    try
    {
      maskImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading reference image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    maskImage = maskImageReader->GetOutput();
    maskImage->Update();
    std::cout << "OK." << std::endl;
  }
  
  if( bSurfDist && maskImage.IsNull() )
  {
    std::cout << "No mask image given, surface distances cannot be calculated!" << std::endl;
    bSurfDist = false;
  }

  ////////////////////////////////////////////////////////////////
  //
  // Generate image to plot landmarks.
  //
  bool bSaveLandmarkFiles = false;
  bool bLandmarksAreValid = true;
  std::string origOutputName;
  std::string::size_type endingPos = 0;
  if( outputImageFilename != NULL )
  {
    bSaveLandmarkFiles = true;
    // to generate new output filenames
    origOutputName = static_cast<const char*> ( outputImageFilename );
    std::string pattern = ".";
    // Get position where ending starts.
    endingPos = origOutputName.rfind( pattern, origOutputName.length() );
    if( endingPos == std::string::npos )
      endingPos = origOutputName.length();
  }
  else
  {
    bSaveLandmarkFiles = false;
  }

  ////////////////////////////////////////////////////////////////
  //
  // Check and plot target landmarks.
  //
  std::cout << "Generating image to check landmarks..." << std::endl;
  ImageType::Pointer landmarkImage = ImageType::New();
  landmarkImage->SetRegions( displField->GetLargestPossibleRegion() );
  landmarkImage->SetSpacing( displField->GetSpacing() );
  landmarkImage->SetOrigin( displField->GetOrigin() );
  landmarkImage->CopyInformation( displField );
  landmarkImage->Allocate();

  std::cout << "Checking target landmarks..." << std::endl;
  landmarkImage->FillBuffer( 0 );
  if( !CheckPointsWithImage( targetPoints, landmarkImage ) )
  {
    std::cout << "Checking failed for Target Landmarks!" << std::endl;
    bLandmarksAreValid = false;
  }

  if( bSaveLandmarkFiles )
  {
    // new target landmark image filename
    std::string outputTargetName = origOutputName;
    outputTargetName.insert( endingPos, "_targ" );

    std::cout << "Saving target landmark image to " << outputTargetName << "..." << std::flush;
    
    ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    imageWriter->SetInput( landmarkImage );
    imageWriter->SetFileName( outputTargetName.c_str() );
    try
    {
      imageWriter->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while writing output image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "OK." << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  //
  // Check and plot reference landmarks.
  //
  std::cout << "Checking reference landmarks..." << std::endl;
  landmarkImage->FillBuffer( 0 );
  landmarkImage->Modified();
  if( !CheckPointsWithImage( refPoints, landmarkImage ) )
  {
    std::cout << "Checking failed for Reference Landmarks!" << std::endl;
    bLandmarksAreValid = false;
  }

  if( bSaveLandmarkFiles )
  {
    // new reference landmark image filename
    std::string outputReferenceName = origOutputName;
    outputReferenceName.insert( endingPos, "_ref" );

    std::cout << "Saving reference landmark image to " << outputReferenceName << "..." << std::flush;
    
    ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    imageWriter->SetInput( landmarkImage );
    imageWriter->SetFileName( outputReferenceName.c_str() );
    try
    {
      imageWriter->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while writing output image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "OK." << std::endl;
  }

  ////////////////////////////////////////////////////////////////
  //
  // Scale displacement field if desired.
  //
  if( bScaleDisplacementField )
  {
    std::cout << "Saving scaled displacement field (factor=" << displFieldScaleFactor << "..." << std::endl;

    typedef itk::MultiplyImageFilter<
        DisplacementFieldType,
        DisplacementFieldType,
        DisplacementFieldType> MultiplyByConstantType;
    MultiplyByConstantType::Pointer multiplier = MultiplyByConstantType::New();
    multiplier->SetConstant( displFieldScaleFactor );
    multiplier->SetInput( displField );
    try
    {
      multiplier->Update();
    }
    catch( itk::ExceptionObject &err )
    {
      std::cerr << "Cannot scale deformation field! Exception error: " << err << std::endl;
      return false;
    }
    displField = multiplier->GetOutput();
  }

  ////////////////////////////////////////////////////////////////
  //
  // Call warping function.
  //
  std::cout << "Warping target landmarks with displacement field..." << std::endl;
  PointsType warpedPoints = WarpLandmarks( refPoints, displField, bSnapToVoxel );

  std::cout << "Checking warped landmarks..." << std::endl;
  landmarkImage->FillBuffer( 0 );
  landmarkImage->Modified();
  if( !CheckPointsWithImage( warpedPoints, landmarkImage ) )
  {
    std::cout << "Checking failed for Warped Landmarks!" << std::endl;
    bLandmarksAreValid = false;
  }

  if( bSaveLandmarkFiles )
  {
    // new warped landmark image filename
    std::string outputWarpedName = origOutputName;
    outputWarpedName.insert( endingPos, "_warp" );

    std::cout << "Saving warped landmark image to " << outputWarpedName << "..." << std::flush;
    
    ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    imageWriter->SetInput( landmarkImage );
    imageWriter->SetFileName( outputWarpedName.c_str() );
    try
    {
      imageWriter->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while writing output image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "OK" << std::endl;
  }


  ////////////////////////////////////////////////////////////////
  //
  // Calculate surface distances.
  //
  std::vector<float> surfDistances;
  PointsType surfacePoints;

  if( bSurfDist )
  {
    std::cout << "Calculating surface distance..." << std::endl;

    BinarySegmentationType::SizeType size = maskImage->GetLargestPossibleRegion().GetSize();
    BinarySegmentationType::SpacingType spacing = maskImage->GetSpacing();
    BinarySegmentationType::IndexType index;
    PointType point;

    for( index[2] = 0; index[2] < static_cast<long>(size[2]); ++index[2] )
    {
      for( index[1] = 0; index[1] < static_cast<long>(size[1]); ++index[1] )
      {
        // Get first point from left.
        index[0] = 0;
        while( maskImage->GetPixel( index ) == 0 && index[0] < static_cast<long>(size[0]) - 1 )
        {
          index[0]++;
        }
        if( index[0] < static_cast<long>(size[0]) - 1 )
        {
          maskImage->TransformIndexToPhysicalPoint( index, point );
          point[0] -= 0.5 * spacing[0];
          surfacePoints.push_back( point );
        }

        // Get first point from right.
        index[0] = size[0] - 1;
        while( maskImage->GetPixel( index ) == 0 && index[0] > 0)
        {
          index[0]--;
        }
        if( index[0] > 0 )
        {
          maskImage->TransformIndexToPhysicalPoint( index, point );
          point[0] += 0.5 * spacing[0];
          surfacePoints.push_back( point );
        }
      }
    }

    for( unsigned int i = 0; i < refPoints.size(); ++i )
    {
      // Calculate distance to closest surface point.
      float distanceToSurfacePoint = itk::NumericTraits<float>::max();
      float minDistanceToSurface = itk::NumericTraits<float>::max();
      for( unsigned int j = 0; j < surfacePoints.size(); ++j )
      {
        distanceToSurfacePoint = refPoints[i].EuclideanDistanceTo( surfacePoints[j] );
        if( distanceToSurfacePoint < minDistanceToSurface )
        {
          minDistanceToSurface = distanceToSurfacePoint;
        }
      }
      surfDistances.push_back( minDistanceToSurface );
    }
  }


  ////////////////////////////////////////////////////////////////
  //
  // Calculate measures.
  //
  std::cout << "Calculating euclidean distances:\n"  << std::endl;

  std::vector<float> distances;
  std::vector<float> distancesOld;

  std::vector<float> distancesInsideMask;
  std::vector<float> distancesInsideMask_x;
  std::vector<float> distancesInsideMask_y;
  std::vector<float> distancesInsideMask_z;
  std::vector<float> distancesInsideMask_Old;
  std::vector<float> distancesOutsideMask;
  std::vector<float> distancesOutsideMask_Old;
  ImageType::IndexType index;

  float dist, distOld;
  float diff, diffOld;
  for( unsigned int i = 0; i < warpedPoints.size(); ++i )
  {
    dist = 0.0;
    distOld = 0.0;
    for( unsigned int d = 0; d < 3; ++d )
    {
      diff = warpedPoints[i][d] - targetPoints[i][d];
      dist += diff * diff;
      diffOld = targetPoints[i][d] - refPoints[i][d];
      distOld += diffOld * diffOld;
    }
    dist = vcl_sqrt( dist );
    distOld = vcl_sqrt( distOld );

    distances.push_back( dist );
    distancesOld.push_back( distOld );

    std::cout << "Distance Point " << i + 1 << ": " << dist <<
                 " [" << warpedPoints[i][0] - targetPoints[i][0] <<
                 ", " << warpedPoints[i][1] - targetPoints[i][1] <<
                 ", " << warpedPoints[i][2] - targetPoints[i][2] << "]; was " << distOld;

    if( bMaskedEvaluation && maskImage.IsNotNull() )
    {
      if( bSurfDist )
      {
        std::cout << " (surf distance " << surfDistances[i] << ")";
      }

      maskImage->TransformPhysicalPointToIndex( refPoints[i], index );
      if( maskImage->GetPixel( index ) != 0 )
      {
        distancesInsideMask.push_back( dist );
        distancesInsideMask_Old.push_back( distOld );
        distancesInsideMask_x.push_back( warpedPoints[i][0] - targetPoints[i][0] );
        distancesInsideMask_y.push_back( warpedPoints[i][1] - targetPoints[i][1] );
        distancesInsideMask_z.push_back( warpedPoints[i][2] - targetPoints[i][2] );
      }
      else
      {
        distancesOutsideMask.push_back( dist );
        distancesOutsideMask_Old.push_back( distOld );
        std::cout << "   WARNING: Point outside mask!";
      }
    }

    std::cout << std::endl;
  }

  float mean = 0.0;
  float meanOld = 0.0;
  float maximum = 0.0;
  float maximumOld = 0.0;

  float meanInsideMask = 0.0;
  float meanInsideMask_Old = 0.0;
  float maximumInsideMask = 0.0;
  float maximumInsideMask_Old = 0.0;

  float meanOutsideMask = 0.0;
  float meanOutsideMask_Old = 0.0;
  float maximumOutsideMask = 0.0;
  float maximumOutsideMask_Old = 0.0;

  //
  // Computation of mean values:
  for( unsigned int i = 0; i < distances.size(); ++i )
  {
    mean += distances[i];
    meanOld += distancesOld[i];
    if( maximum < distances[i] )
      maximum = distances[i];
    if( maximumOld < distancesOld[i] )
      maximumOld = distancesOld[i];
  }
  mean /= distances.size();
  meanOld /= distancesOld.size();
  std::cout << "Mean Distance: " << mean << " (was " << meanOld << ")" << std::endl;

  if( bMaskedEvaluation )
  {
    for( unsigned int i = 0; i < distancesInsideMask.size(); ++i )
    {
      meanInsideMask += distancesInsideMask[i];
      meanInsideMask_Old += distancesInsideMask_Old[i];
      if( maximumInsideMask < distancesInsideMask[i] )
        maximumInsideMask = distancesInsideMask[i];
      if( maximumInsideMask_Old < distancesInsideMask_Old[i] )
        maximumInsideMask_Old = distancesInsideMask_Old[i];
    }
    if( distancesInsideMask.size() != 0 )
    {
      meanInsideMask /= distancesInsideMask.size();
      meanInsideMask_Old /= distancesInsideMask_Old.size();
      std::cout << "Mean Distance inside Mask: " << meanInsideMask << " (was " << meanInsideMask_Old << "; no of Points: " << distancesInsideMask.size() << ")" << std::endl;
    }
    else
    {
      std::cout << "No landmarks inside Mask!" << std::endl;
    }

    for( unsigned int i = 0; i < distancesOutsideMask.size(); ++i )
    {
      meanOutsideMask += distancesOutsideMask[i];
      meanOutsideMask_Old += distancesOutsideMask_Old[i];
      if( maximumOutsideMask < distancesOutsideMask[i] )
        maximumOutsideMask = distancesOutsideMask[i];
      if( maximumOutsideMask_Old < distancesOutsideMask_Old[i] )
        maximumOutsideMask_Old = distancesOutsideMask_Old[i];
    }
    if( distancesOutsideMask.size() != 0 )
    {
      meanOutsideMask /= distancesOutsideMask.size();
      meanOutsideMask_Old /= distancesOutsideMask_Old.size();
      std::cout << "Mean Distance outside Mask: " << meanOutsideMask << " (was " << meanOutsideMask_Old << "; no of Points: " << distancesOutsideMask.size() << ")" << std::endl;
    }
    else
    {
      std::cout << "No landmarks outside mask!" << std::endl;
    }
  }

  //
  // Computation of variances:
  float variance = 0.0;
  float varianceOld = 0.0;

  float varianceInsideMask = 0.0;
  float varianceInsideMask_Old = 0.0;
  float varianceOutsideMask = 0.0;
  float varianceOutsideMask_Old = 0.0;

  for( unsigned int i = 0; i < distances.size(); ++i )
  {
    variance += (distances[i] - mean) * (distances[i] - mean);
    varianceOld += (distancesOld[i] - meanOld) * (distancesOld[i] - meanOld);
  }
  variance /= distances.size();
  varianceOld /= distancesOld.size();
  std::cout << "Variance: " << variance << " (was " << varianceOld << ")" << std::endl;
  std::cout << " Std-Abw: " << vcl_sqrt( variance ) << " (was " << vcl_sqrt( varianceOld ) << ")" << std::endl;
  std::cout << " Maximum: " << maximum << " (was " << maximumOld << ")" << std::endl;

  if( bMaskedEvaluation )
  {
    // Computation of LM-Variances for LMs inside mask:
    for( unsigned int i = 0; i < distancesInsideMask.size(); ++i )
    {
      varianceInsideMask += (distancesInsideMask[i] - meanInsideMask) * (distancesInsideMask[i] - meanInsideMask);
      varianceInsideMask_Old += (distancesInsideMask_Old[i] - meanInsideMask_Old) * (distancesInsideMask_Old[i] - meanInsideMask_Old);
    }
    if( distancesInsideMask.size() != 0 )
    {
      varianceInsideMask /= distancesInsideMask.size();
      varianceInsideMask_Old /= distancesInsideMask_Old.size();
    }

    // Computation of LM-Variances for LMs outside mask:
    for( unsigned int i = 0; i < distancesOutsideMask.size(); ++i )
    {
      varianceOutsideMask += (distancesOutsideMask[i] - meanOutsideMask) * (distancesOutsideMask[i] - meanOutsideMask);
      varianceOutsideMask_Old += (distancesOutsideMask_Old[i] - meanOutsideMask_Old) * (distancesOutsideMask_Old[i] - meanOutsideMask_Old);
    }
    if( distancesOutsideMask.size() != 0 )
    {
      varianceOutsideMask /= distancesOutsideMask.size();
      varianceOutsideMask_Old /= distancesOutsideMask_Old.size();
    }

    // Output:
    if( distancesInsideMask.size() != 0 )
    {
      std::cout << "Variance inside Mask: " << varianceInsideMask << " (was " << varianceInsideMask_Old << "; no of Points: " << distancesInsideMask.size() << ")" << std::endl;
    }
    if( distancesOutsideMask.size() != 0 )
    {
      std::cout << "Variance outside Mask: " << varianceOutsideMask << " (was " << varianceOutsideMask_Old << "; no of Points: " << distancesOutsideMask.size() << ")" << std::endl;
    }
    if( distancesInsideMask.size() != 0 )
    {
      std::cout << "Maximum inside Mask:  " << maximumInsideMask << " [" << *std::min_element( distancesInsideMask_x.begin(), distancesInsideMask_x.end() ) <<
                                                                    "/" << *std::max_element( distancesInsideMask_x.begin(), distancesInsideMask_x.end() ) << "; "
                                                                    " [" << *std::min_element( distancesInsideMask_y.begin(), distancesInsideMask_y.end() ) <<
                                                                    "/" << *std::max_element( distancesInsideMask_y.begin(), distancesInsideMask_y.end() ) << "; "
                                                                    " [" << *std::min_element( distancesInsideMask_z.begin(), distancesInsideMask_z.end() ) <<
                                                                    "/" << *std::max_element( distancesInsideMask_z.begin(), distancesInsideMask_z.end() ) << "]; was " << maximumInsideMask_Old << std::endl;
    }
    if( distancesOutsideMask.size() != 0 )
    {
      std::cout << "Maximum outside Mask:  " << maximumOutsideMask << " (was " << maximumOutsideMask_Old << ")" << std::endl;
    }
  }

  //
  // Calculate quantiles and outliers.
  std::vector<float> sortedDistances;
  sortedDistances.resize( distances.size());

  std::copy( distances.begin(), distances.end(), sortedDistances.begin() );
  std::sort( sortedDistances.begin(), sortedDistances.end() );

  unsigned int q1index = 0.25 * sortedDistances.size();
  unsigned int q3index = 0.75 * sortedDistances.size();

  float q1 = sortedDistances[q1index];
  float q3 = sortedDistances[q3index];

  const float K = 1.5;
  float lowerWedge = q1 - K * (q3 - q1);
  float upperWedge = q3 + K * (q3 - q1);

  std::cout << "Q1: " << q1 << " Q3: " << q3 << " LW: " << lowerWedge << " UW: " << upperWedge << std::endl;

  if( bSurfDist )
  {
    std::cout << "  " << refPoints.size() << "  " << distances.size() << "  " << surfDistances.size() <<  std::endl;

    // Calculate outlier statistics an image.
    landmarkImage->FillBuffer( 0 );
    landmarkImage->Modified();

    float meanSurfaceDistanceAll = 0.0;
    float meanSurfaceDistanceOutliers = 0.0;
    float meanSurfaceDistanceInliers = 0.0;
    unsigned int numberOfOutliers = 0;

    for( unsigned int i = 0; i < refPoints.size(); ++i )
    {
      meanSurfaceDistanceAll += surfDistances[i];

      // Check if point is an outlier.
      if( distances[i] > upperWedge )
      {
        numberOfOutliers++;
        meanSurfaceDistanceOutliers += surfDistances[i];

        ImageType::IndexType index;
        landmarkImage->TransformPhysicalPointToIndex( refPoints[i], index );
        landmarkImage->SetPixel( index, 255 );
      }
      else
      {
        meanSurfaceDistanceInliers += surfDistances[i];
      }

    }

    meanSurfaceDistanceAll /= refPoints.size();
    meanSurfaceDistanceInliers /= (refPoints.size() - numberOfOutliers);
    meanSurfaceDistanceOutliers /= numberOfOutliers;

    std::cout << "Mean surface distance all:      " << meanSurfaceDistanceAll << std::endl;
    std::cout << "Mean surface distance outliers: " << meanSurfaceDistanceOutliers << std::endl;
    std::cout << "Mean surface distance inliers:  " << meanSurfaceDistanceInliers << std::endl;
    std::cout << "Number of outliers:             " << numberOfOutliers << std::endl;
  }

  if( bSaveLandmarkFiles )
  {
    // new warped landmark image filename
    std::string outputWarpedName = origOutputName;
    outputWarpedName.insert( endingPos, "_refoutliers" );

    std::cout << "Saving warped landmark image to " << outputWarpedName << "..."  << std::flush;
    
    ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    imageWriter->SetInput( landmarkImage );
    imageWriter->SetFileName( outputWarpedName.c_str() );
    try
    {
      imageWriter->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while writing output image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "OK." << std::endl;
  }

  //
  // Writing Log file:
  if( logFilename != NULL )
  {
    std::ofstream logFile( logFilename );
    std::cout << "Saving log file..." << std::endl;

    logFile << "Distance \t Distance(old) \t warpedPoints - - \t referencePoints - - \t targetPoints - - \n";
    logFile << "Distance \t Distance(old) \t x y z \t x y z \t x y z \n";
    logFile << "\n";
    for( unsigned int i = 0; i < warpedPoints.size(); ++i )
    {
      logFile << distances[i] << "\t" << distancesOld[i] << "\t";
      logFile << warpedPoints[i][0] << " " << warpedPoints[i][1] << " " << warpedPoints[i][2] << "\t";
      logFile << refPoints[i][0] << " " << refPoints[i][1] << " " << refPoints[i][2] << "\t";
      logFile << targetPoints[i][0] << " " << targetPoints[i][1] << " " << targetPoints[i][2] << "\t";
      if( bSurfDist )
      {
        logFile << surfDistances[i] << "\t";
      }
      logFile << "\n";
    }
    logFile.close();
  }

  if( bLandmarksAreValid )
  {
    std::cout << "" << std::endl;
    std::cout << "All landmarks are valid." << std::endl;
  }
  else
  {
    std::cout << "" << std::endl;
    std::cerr << "Landmarks are not valid!" << std::endl;
  }

  std::cout << "icnsValidateLandmarks FINISHED." << std::endl;
  std::cout << "==========================================\n" << std::endl;
  
  return EXIT_SUCCESS;
}
