/** \file icnsAffineImageTransformation.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2016 Department of Computational Neuroscience,
 *     University Medical Center Hamburg-Eppendorf
 *
 ****************************************************************************/

// System includes
#include <iostream>
extern "C"
{
#include "getopt.h"
}

#include <itkAffineTransform.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkResampleImageFilter.h>
#include <itkTransformFileReader.h>

// Project includes: NONE so far.

// Global typedefs:
typedef itk::AffineTransform<double, 3>      AffineTransformType;
typedef itk::Image<short, 3>                 ImageType;
typedef itk::ImageFileReader<ImageType>      ImageReaderType;
typedef itk::ImageFileWriter<ImageType>      ImageWriterType;
typedef itk::NearestNeighborInterpolateImageFunction<ImageType> NNInterpolatorType;
typedef itk::LinearInterpolateImageFunction<ImageType>          LinearInterpolatorType;
typedef itk::BSplineInterpolateImageFunction<ImageType>         BSplineInterpolatorType;


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
    std::cout << "\n";
    std::cout << "Usage :\n";
    std::cout << "icnsAffineImageRegistration -I <input image> -T <transformation> -O  <transformed image> [-i <interpolation approach]\n\n";
    
    std::cout << "-I <input image>               Filename of the input image.\n";
    std::cout << "-T <transformation>            Transformation filename.\n";
    std::cout << "-O <output image>              Filename of the output (=transformed) image.\n";
    std::cout << "-i [0|1|2]                     Interpolation approach.\n";
    std::cout << "                                 0: nearest neighbor (choose for binary data!).\n";
    std::cout << "                                 1: linear (default).\n";
    std::cout << "                                 2: cubic B-splines.\n";
  
    std::cout << "-h                             Print this help.\n";
    std::cout << "\n";
}


// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  std::cout << "========================================" << std::endl;
  std::cout << "icnsAffineRegistration         " << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputImageFilename       = NULL;
  char* transformationFilename   = NULL;
  char* outputImageFilename      = NULL;
  
  short backgroundIntensityValue = 0;
  short interpolationApproach    = 1; // 0: NN; 1: linear; 2: cubic B-Splines
  
  // Reading parameters: 
  while( (c = getopt( argc, argv, "I:T:O:i:?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputImageFilename = optarg;
        std::cout << "  Input image filename:              " << inputImageFilename << std::endl;
        break;
      case 'T':
        transformationFilename = optarg;
        std::cout << "  Transformation filename:           " << transformationFilename << std::endl;
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << "  Output image filename:             " << outputImageFilename << std::endl;
        break;
      case 'i':
        interpolationApproach = atoi( optarg );
        if( interpolationApproach == 0 )
        {
          std::cout << "  Interpolation approach:            NEAREST NEIGHBOR" << std::endl;
        }
        else if( interpolationApproach == 2 )
        {
          std::cout << "  Interpolation approach:            CUBIC B-SPLINES" << std::endl;
        }
        else
        {
          std::cout << "  Interpolation approach:            LINEAR" << std::endl;
        }
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_FAILURE;
        break;
      default:
        std::cout << "  Argument " << (char)c << " not processed!\n" << std::endl;
        return EXIT_FAILURE;
    }
  }
  
  // Check validity of arguments:
  
  if( inputImageFilename == NULL )
  {
    std::cerr << "ERROR: No input image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( transformationFilename == NULL )
  {
    std::cerr << "ERROR: No output image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputImageFilename == NULL )
  {
    std::cerr << "ERROR: No output image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input data:"                      << std::endl;
  
  // Loading input image:
  
  std::cout << "  Input image ... " << std::flush;
  
  ImageReaderType::Pointer inputImageReader = ImageReaderType::New();
  inputImageReader->SetFileName( inputImageFilename );
  try
  {
    inputImageReader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading input image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  
  // Loading transform:
  
  AffineTransformType::Pointer affineTransform;
  
  std::cout << "  Transform file ... "  << std::flush;
  
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileReaderTemplate<double>::Pointer transformReader = itk::TransformFileReaderTemplate<double>::New();
#else
  itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
#endif
  transformReader->SetFileName( transformationFilename );
  try
  {
    transformReader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading transformation: " << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << " OK." << std::endl;
  
  // Check if read transform is an affine transform.
  // If so, assign it to the affineTransform pointer.
  
  typedef itk::TransformFileReaderTemplate<double> TransformReaderType;
  typedef TransformReaderType::TransformListType* TransformListType;
  
  TransformListType transforms = transformReader->GetTransformList();
  std::cout << "   Number of transforms = " << transforms->size() << std::endl;
  
  TransformReaderType::TransformListType::const_iterator it = transforms->begin();
  if( !strcmp((*it)->GetNameOfClass(), "AffineTransform") )
  {
    affineTransform = static_cast< AffineTransformType* >( (*it).GetPointer() );
    affineTransform->Print(std::cout);
  }
  else
  {
    std::cerr << "  ERROR: transformation type is not AffineTransform!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Applying transformation to input image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Applying transformation ... " << std::flush;

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  
  typename ResampleFilterType::Pointer resample = ResampleFilterType::New();
  resample->SetTransform( affineTransform );
  resample->SetInput( inputImageReader->GetOutput() );
  resample->SetReferenceImage( inputImageReader->GetOutput() );
  resample->UseReferenceImageOn();
  resample->SetDefaultPixelValue( backgroundIntensityValue );

  if( interpolationApproach == 0 )
  {
    // generate a nearest neighbor interpolator
    NNInterpolatorType::Pointer interpolator = NNInterpolatorType::New();
    resample->SetInterpolator( interpolator );
  }
  else if( interpolationApproach == 2 )
  {
    // generate a B-spline interpolator
    BSplineInterpolatorType::Pointer interpolator = BSplineInterpolatorType::New();
    resample->SetInterpolator( interpolator );
  }
  else
  {
    // generate a linear interpolator
    LinearInterpolatorType::Pointer interpolator = LinearInterpolatorType::New();
    resample->SetInterpolator( interpolator );
  }

  resample->Update();
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetInput( resample->GetOutput() );
  imageWriter->SetFileName( outputImageFilename );
  try
  {
    imageWriter->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "   ERROR while writing warped target image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;
  
  return EXIT_SUCCESS;
}
