/** \file icnsSegmentationAccuracyEvaluation
 *
 *  \b Initial \b Author: Jan Ehrhardt\n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// System includes:
#include <string>
#include <iostream>
#include <iomanip>
extern "C"
{
#include "getopt.h"
}

// ITK includes
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkLabelOverlapMeasuresImageFilter.h>
//#include "itkSubtractImageFilter.h"
//#include "itkDifferenceImageFilter.h"
//#include "itkHausdorffDistanceImageFilter.h"
//#include "itkSignedMaurerDistanceMapImageFilter.h"
//#include "itkSignedDanielssonDistanceMapImageFilter.h"
//#include "itkNeighborhoodIterator.h"

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 3>                 ImageType;
typedef itk::ImageFileReader<ImageType>      ImageReaderType;
typedef itk::LabelOverlapMeasuresImageFilter<ImageType> LabelOverlapMeasuresFilterType;
//typedef itk::ImageFileWriter<ImageType>      ImageWriterType;
//typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> BinaryThresholdFilterType;
//typedef itk::IdentityTransform<double, ImageType::ImageDimension> TransformType;
//typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
//typedef itk::NearestNeighborInterpolateImageFunction<ImageType, double> NNInterpolatorType;
//typedef itk::BSplineInterpolateImageFunction<ImageType, double> BSplineInterpolatorType;
//typedef itk::LinearInterpolateImageFunction<ImageType, double> LinearInterpolatorType;

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsSegmentationAccuracyEvaluation -R <reference segmentation> -T <test segmentation> [...]" << std::endl;
  std::cout << std::endl;
  std::cout << "-R <reference segmentation>    Filename of reference segmentation." << std::endl;
  std::cout << "-T <test segmentation>         Filename of segmentation to be evaluated." << std::endl;
  std::cout << "-L <log file>                  Filename of log file." << std::endl;
  std::cout << "-d [0|1]                       Compute distance measures: Hausdorff distance (NYI)." << std::endl;
  std::cout << "-l                             Lower threshold of intensity interval considered as INSIDE segmentation values (NYI)." << std::endl;
  std::cout << "-h                             Upper threshold of intensity interval considered as INSIDE segmentation values (NYI)." << std::endl;
  std::cout << "-h                             Print this help." << std::endl;
  std::cout << std::endl;
}



// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 2 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  std::cout << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsEvaluateSegmentationAccuracy"           << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;
    
  // Initializing parameters with default values:
    
  int c;
  
  char* referenceSegmentationFilename = NULL;
  char* testSegmentationFilename      = NULL;
  char* logFilename                   = NULL;
  bool computeDistanceMeasures        = false;

  // Reading parameters: 
  
  while( (c = getopt( argc, argv, "R:T:L:dh:l:h?" )) != -1 )
  {
    switch( c )
    {
      case 'R':
        referenceSegmentationFilename = optarg;
        std::cout << "  Reference segmentation filename:   " << referenceSegmentationFilename << std::endl;
        break;
      case 'T':
        testSegmentationFilename = optarg;
        std::cout << "  Test segmentation filename:        " << testSegmentationFilename << std::endl;
        break;
      case 'L':
        logFilename = optarg;
        std::cout << "  Logfile filename:                  " << logFilename << " (FYI)." << std::endl;
        break;
      case 'd':
        computeDistanceMeasures = true;
        std::cout << "  Computing distance measures:       TRUE (NYI)." << std::endl;
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_SUCCESS;
      default:
        std::cout << "  Argument " << (char)c << " not processed!\n" << std::endl;
        return EXIT_FAILURE;
    }
  }

  // Check validity of arguments:
  
  if( referenceSegmentationFilename == NULL )
  {
    std::cerr << "ERROR: No reference segmentation filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( testSegmentationFilename == NULL )
  {
    std::cerr << "ERROR: No test segmentation filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input images:
  // (1/2) Reference segmentation:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading reference segmentation ... "     << std::flush;
  
  ImageReaderType::Pointer referenceSegmentationReader = ImageReaderType::New();
  referenceSegmentationReader->SetFileName( referenceSegmentationFilename );
  try
  {
    referenceSegmentationReader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading reference segmentation." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  
  // (2/2) Test segmentation:
  
  ImageReaderType::Pointer testSegmentationReader;
  std::cout << "Loading test segmentation ... " << std::flush;
    
  testSegmentationReader = ImageReaderType::New();
  testSegmentationReader->SetFileName( testSegmentationFilename );
  try
  {
    testSegmentationReader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading test segmentation." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
    
  std::cout << "OK." << std::endl;
  
  // Init pipeline input data:
  
  ImageType::Pointer referenceSegmentationImage;
  referenceSegmentationImage = referenceSegmentationReader->GetOutput();
  referenceSegmentationImage->Update();
  referenceSegmentationImage->DisconnectPipeline();
  
  ImageType::Pointer testSegmentationImage;
  testSegmentationImage = testSegmentationReader->GetOutput();
  testSegmentationImage->Update();
  testSegmentationImage->DisconnectPipeline();

  // -------------------------------------------------------------
  // Compute overlap measures: using itkLabelOverlapMeasuresImageFilter
  
  LabelOverlapMeasuresFilterType::Pointer labelOverlapMeasuresFilter = LabelOverlapMeasuresFilterType::New();
  labelOverlapMeasuresFilter->SetSourceImage( referenceSegmentationImage );
  labelOverlapMeasuresFilter->SetTargetImage( testSegmentationImage );
  labelOverlapMeasuresFilter->Update();
  
  std::cout << "                                       "
            << "************ Individual Labels *************" << std::endl;
  std::cout << std::setw( 10 ) << "Label"
            << std::setw( 17 ) << "Target"
            << std::setw( 17 ) << "Union (jaccard)"
            << std::setw( 17 ) << "Mean (dice)"
            << std::setw( 17 ) << "Volume sim."
            << std::setw( 17 ) << "False negative"
            << std::setw( 17 ) << "False positive" << std::endl;
  
  LabelOverlapMeasuresFilterType::MapType labelMap = labelOverlapMeasuresFilter->GetLabelSetMeasures();
  LabelOverlapMeasuresFilterType::MapType::const_iterator it;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
  {
    if( (*it).first == 0 ) continue;
    
    int label = (*it).first;
    
    std::cout << std::setw( 10 ) << label;
    std::cout << std::setw( 17 ) << labelOverlapMeasuresFilter->GetTargetOverlap( label );
    std::cout << std::setw( 17 ) << labelOverlapMeasuresFilter->GetUnionOverlap( label );
    std::cout << std::setw( 17 ) << labelOverlapMeasuresFilter->GetMeanOverlap( label );
    std::cout << std::setw( 17 ) << labelOverlapMeasuresFilter->GetVolumeSimilarity( label );
    std::cout << std::setw( 17 ) << labelOverlapMeasuresFilter->GetFalseNegativeError( label );
    std::cout << std::setw( 17 ) << labelOverlapMeasuresFilter->GetFalsePositiveError( label );
    std::cout << std::endl;
  }
  
/*
  imiINFO("--------------------------------------");
  imiINFO("RESULTS:");
  imiINFO("       DICE : "<<calculateOverlapMeasures->GetDiceCoefficient() );
  imiINFO("    JACCARD : "<<calculateOverlapMeasures->GetJaccardCoefficient() );
  imiINFO("   Center A : "<<calculateOverlapMeasures->GetCenter1() );
  imiINFO("   Center B : "<<calculateOverlapMeasures->GetCenter2() );
  imiINFO(" CenterDiff : "<<calculateOverlapMeasures->GetCenterDistance() );
  imiINFO("   Volume A : "<<calculateOverlapMeasures->GetNumberOfVoxels1() <<" ( "<<calculateOverlapMeasures->GetVolume1()<<" cm^3 )" );
  imiINFO("   Volume B : "<<calculateOverlapMeasures->GetNumberOfVoxels2() <<" ( "<<calculateOverlapMeasures->GetVolume2()<<" cm^3 )" );

  if(!bComputeHausdorffDistance)
  {
    // Open log file if necessary.
    if( logFilename != NULL )
    {
      std::ofstream logfile;

      logfile.open( logFilename, std::ios::app );

      logfile << "# Jaccard\tCenter Diff" << std::endl;
      logfile << calculateOverlapMeasures->GetJaccardCoefficient() << "\t" << calculateOverlapMeasures->GetCenterDistance()<< std::endl << std::endl;

      logfile.close();
    }
  }
  else
  {
    typedef itk::HausdorffDistanceImageFilter<ImageShort3DType,ImageShort3DType> HausdorffDistanceImageFilterType;

    HausdorffDistanceImageFilterType::Pointer hausdorffDistanceFilter = HausdorffDistanceImageFilterType::New();
    hausdorffDistanceFilter->SetInput1(inputImage1);
    hausdorffDistanceFilter->SetInput2(inputImage2);
    hausdorffDistanceFilter->Update();

    double hausdorffDistance = hausdorffDistanceFilter->GetHausdorffDistance();

    imiINFO("   Hausdorff: "<<hausdorffDistance);

    // Open log file if necessary.
    if( logFilename != NULL )
    {
      std::ofstream logfile;

      logfile.open( logFilename, std::ios::app );

      logfile << "# Jaccard\tCenter Diff\tHausdorff" << std::endl;
      logfile << calculateOverlapMeasures->GetJaccardCoefficient() << "\t" << calculateOverlapMeasures->GetCenterDistance()<< "\t" <<hausdorffDistance << std::endl;

      logfile.close();
    }
  }

  if( bComputeSurfaceDistance )
  {
    //typedef itk::SignedMaurerDistanceMapImageFilter<ImageShort3DType, ImageFloat3DType> DistanceMapFilterType;
    typedef itk::SignedDanielssonDistanceMapImageFilter<ImageShort3DType, ImageFloat3DType> DistanceMapFilterType;
    DistanceMapFilterType::Pointer distanceFilter1 = DistanceMapFilterType::New();
    DistanceMapFilterType::Pointer distanceFilter2 = DistanceMapFilterType::New();

    distanceFilter1->SetInput( inputImage1 );
    distanceFilter1->UseImageSpacingOff();
    distanceFilter1->SquaredDistanceOff();
    distanceFilter1->Update();

    ImageFloat3DTypePointer distanceImage1 = distanceFilter1->GetOutput();
    imiImageFunctions::WriteImage( distanceImage1, "dist1.mha" );

    distanceFilter2->SetInput( inputImage2 );
    distanceFilter2->UseImageSpacingOff();
    distanceFilter2->SquaredDistanceOff();
    distanceFilter2->Update();

    ImageFloat3DTypePointer distanceImage2 = distanceFilter2->GetOutput();
    imiImageFunctions::WriteImage( distanceImage2, "dist2.mha" );

    unsigned int ctr1 = 0;
    double meanDistance1 = 0.0;

    unsigned int ctr2 = 0;
    double meanDistance2 = 0.0;

    typedef itk::ImageRegionConstIterator<ImageFloat3DType> IteratorType;
    IteratorType it1( distanceImage1, distanceImage1->GetRequestedRegion() );

    for( it1.GoToBegin(); !it1.IsAtEnd(); ++it1 )
    {
      if( it1.Get() <= 0.5 )
      {
        ctr1++;
        meanDistance1 += vnl_math_abs( it1.Get() - distanceImage2->GetPixel( it1.GetIndex() ) );
      }
    }

    IteratorType it2( distanceImage2, distanceImage2->GetRequestedRegion() );

    for( it2.GoToBegin(); !it2.IsAtEnd(); ++it2 )
    {
      if( it2.Get() <= 0.5 )
      {
        ctr2++;
        meanDistance2 += vnl_math_abs( it2.Get() - distanceImage1->GetPixel( it2.GetIndex() ) );
      }
    }
    meanDistance1 /= ctr1;
    meanDistance2 /= ctr2;

    std::cout << "Surface distance: "<< meanDistance1 << " " << meanDistance2 << " " << (meanDistance1 + meanDistance2) /2.0 << "\n";
  }*/

  std::cout << "==========================================" << std::endl;
  return EXIT_SUCCESS;
}
