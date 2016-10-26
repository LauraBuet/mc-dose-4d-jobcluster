/** \file imiGenericMLRPredictionMain.cxx
 *
 *  Original authors:
 *
 *  \b Initial \b Authors: Matthias Wilms, Rene Werner \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 *  Modified by: Rene Werner (ICNS, 2016)
 *
 ****************************************************************************/

// System includes:
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
extern "C"
{
#include "getopt.h"
}

// ITK includes:
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkExponentialDisplacementFieldImageFilter.h>
#include <itkWarpImageFilter.h>

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "icnsMLRMotionPrediction.h"
#include "icnsSurrogateBasedMotionPredictionTypeDefinitions.h"

// Global typedefs:
typedef itk::Image<short, 3>                                                 ImageType;
typedef itk::ImageFileReader<ImageType>                                      ImageReaderType;
typedef itk::ImageFileWriter<ImageType>                                      ImageWriterType;
typedef itk::LinearInterpolateImageFunction< ImageType, double >             InterpolatorType;
typedef itk::WarpImageFilter< ImageType, ImageType, DisplacementFieldType >  WarpImageFilterType;

using namespace imi;

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "icnsGenericMLRPrediction ... \n";

  std::cout << "-V_mean       IN: Mean regressand in Matlab format (.mat as generated during training).\n";
  std::cout << "-Z_mean       IN: Mean regressor in Matlab format (.mat as generated during training.\n";
  std::cout << "-B            IN: Filename of system matrix (.mat).\n";
  std::cout << "-Z_hat        IN: Input measurement (.mat in -v4 format).\n";
  std::cout << "-M            IN (optional): Mask file to restrict the computations.\n";
  std::cout << "              only usable with -v_hat and *.mha.";
  
  std::cout << "-V_hat        OUT: Prediction output in matlab format (.mat).\n";
  std::cout << "-v_hat        OUT: Prediction output field (.mha).\n";
  
  std::cout << "-h            Print this help.\n";
}

// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  std::cout << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsGenericMLRPrediction                  " << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;
  
  if( argc < 2 )
  {
    PrintHelp();
    return EXIT_SUCCESS;
  }

  // Initializing parameters with default values:

  std::string meanRegressorFilename;
  std::string meanRegressandFilename;
  std::string inSystemMatrixFilename;
  std::string inputMatrixFilename;
  std::string predictionMatrixFilename;
  std::string maskFilename;
  std::string templateFilename;
  std::string inputImageFilename;

  bool inverseField=false;

  // Running through cmd line params:

  std::cout << std::endl;
  for( int i = 1; i < argc; i++ )
  {
    if( strcmp( argv[i], "-h" ) == 0 )
    {
      PrintHelp();
      return EXIT_SUCCESS;
    }

    // Regressor related params:
    if( strcmp( argv[i], "-Z" ) == 0 )
    {
      i++;
      meanRegressorFilename = argv[i];
      std::cout << "Mean regressor (Z_mean) filename:        " << meanRegressorFilename << std::endl;
      continue;
    }

    // Regressand related params:
    if( strcmp( argv[i], "-V_mean" ) == 0 )
    {
      i++;
      meanRegressandFilename = argv[i];
      std::cout << "Mean regressand (V_mean) filename:        " << meanRegressandFilename << std::endl;
      continue;
    }

    // Sytem matrix related params:
    if( strcmp( argv[i], "-B" ) == 0 )
    {
      i++;
      inSystemMatrixFilename = argv[i];
      std::cout << "System matrix (B) filename: " << inSystemMatrixFilename << std::endl;
      continue;
    }

    // Input matrix related params:
    if( strcmp( argv[i], "-I" ) == 0 )
    {
      i++;
      inputMatrixFilename = argv[i];
      std::cout << "Input measurement (Z_hat) filename: " << inputMatrixFilename << std::endl;
      continue;
    }

    // Prediction matrix related params:
    if( strcmp( argv[i], "-O" ) == 0 )
    {
      i++;
      predictionMatrixFilename = argv[i];
      std::cout << "Prediction (V_hat) filename: " << predictionMatrixFilename << std::endl;
      continue;
    }

    //Mask
    if( strcmp( argv[i], "-M" ) == 0 )
    {
      i++;
      maskFilename = argv[i];
      std::cout << "Mask filename:        " << maskFilename << std::endl;
      continue;
    }

    //Template
    if( strcmp( argv[i], "-T" ) == 0 )
    {
      i++;
      templateFilename = argv[i];
      std::cout << "Template filename:        " << templateFilename << std::endl;
      continue;
    }

    //Input Image
    if( strcmp( argv[i], "-W" ) == 0 )
    {
      i++;
      inputImageFilename = argv[i];
      std::cout << "Input image filename:        " << inputImageFilename << std::endl;
      continue;
    }

    //Inverse Field?
    if( strcmp( argv[i], "-o" ) == 0 )
    {
      inverseField=true;
      std::cout << "Inverse field:        " << "true" << std::endl;
      continue;
    }

    // Anything not accounted for?
    std::cout << "Unknown param at " << i << ": " << argv[i] << std::endl;
  }

  icnsMLRMotionPrediction* predictionFilter = icnsMLRMotionPrediction::New();

  // -------------------------------------------------------------
  // Prediction step
  // Assuming regressor observation given as .mat.
  // Prediction will be saved unter specified name.
  // -------------------------------------------------------------

  vcl_ifstream fid;

  VnlMatrixType meanRegressor;
  VnlMatrixType meanRegressand;
  VnlMatrixType systemMatrix;
  VnlMatrixType measurement;
  VnlMatrixType prediction;

  if( meanRegressorFilename.empty() )
  {
    std::cerr << " No mean regressor specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "  Reading mean regressor (Z_mean) ..." << std::endl;
    fid.open( meanRegressorFilename.c_str() );
    vnl_matlab_read_or_die( fid, meanRegressor );
    fid.close();
  }

  if( meanRegressandFilename.empty() )
  {
    std::cerr << " No mean regressand specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "  Reading mean regressand (V_mean) ..." << std::endl;
    fid.open( meanRegressandFilename.c_str() );
    vnl_matlab_read_or_die( fid, meanRegressand );
    fid.close();
  }

  if( inSystemMatrixFilename.empty() )
  {
    std::cerr << " No system matrix specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "  Reading system matrix (B) ..." << std::endl;
    fid.open( inSystemMatrixFilename.c_str() );
    vnl_matlab_read_or_die( fid, systemMatrix );
    fid.close();
  }

  if( inSystemMatrixFilename.empty() )
  {
    std::cerr << " No measurement specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "  Reading measurement (Z_hat) ..." << std::endl;
    fid.open( inputMatrixFilename.c_str() );
    vnl_matlab_readhdr matlabHeader( fid );
    
    // Depending on type, use the appropriate reading function:
    if( matlabHeader.is_single() )
    {
      icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileRealType( measurement, inputMatrixFilename );
    }
    else
    {
      icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble( measurement, inputMatrixFilename );
    }
    fid.close();
  }

  // -------------------------------------------------------------
  // Actual prediction:
  // -------------------------------------------------------------
  
  predictionFilter->SetMeanRegressor( meanRegressor );
  predictionFilter->SetMeanRegressand( meanRegressand );
  predictionFilter->SetTrainedEstimator( systemMatrix );
  predictionFilter->PredictOutput( measurement, prediction );
  
  // -------------------------------------------------------------
  // Writing results:
  // -------------------------------------------------------------
  
  if( !predictionMatrixFilename.empty() && inputImageFilename.empty() )
  {
    std::cout << "Saving prediction output ... " << std::endl;
    icnsSurrogateBasedMotionPredictionHelpers::SaveVnlVectorAsVectorField( prediction, predictionMatrixFilename, templateFilename, maskFilename );

  }
  else if (!inputImageFilename.empty())
  {
    DisplacementFieldPointerType outputField;
    icnsSurrogateBasedMotionPredictionHelpers::GenerateVectorFieldFromVnlVector( prediction, outputField, templateFilename, maskFilename );

    std::cout << "Allocate displacement field ... " << std::endl;

    DisplacementFieldPointerType displacementField;

    displacementField = DisplacementFieldType::New();
    displacementField->CopyInformation( outputField );
    displacementField->SetBufferedRegion( outputField->GetBufferedRegion() );
    displacementField->SetRequestedRegion( outputField->GetRequestedRegion() );
    displacementField->Allocate();

    typedef itk::ExponentialDisplacementFieldImageFilter<DisplacementFieldType, DisplacementFieldType>  FieldExponentiatorType;
    
    FieldExponentiatorType::Pointer velocityFieldExponentiator = FieldExponentiatorType::New();
    velocityFieldExponentiator->SetInput( outputField );
    velocityFieldExponentiator->AutomaticNumberOfIterationsOff();
    velocityFieldExponentiator->SetMaximumNumberOfIterations( 4 );

    if( inverseField )
    {
      std::cout << "Calculate inverse displacement field ... " << std::flush;
      velocityFieldExponentiator->ComputeInverseOn();
      displacementField = velocityFieldExponentiator->GetOutput();
      displacementField->Update();
    }
    else
    {
      std::cout << "Calculate displacement field ... " << std::flush;
      displacementField = velocityFieldExponentiator->GetOutput();
      displacementField->Update();
    }
    std::cout << "OK." << std::endl;

    // If an image to load is specified, load it:
    
    std::cout << "Loading input image..." << std::flush;
    
    ImageReaderType::Pointer inputImageReader = ImageReaderType::New();
    inputImageReader->SetFileName( inputImageFilename.c_str() );
    try
    {
      inputImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "ERROR while loading input image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    ImageType::Pointer inputImage = inputImageReader->GetOutput();
    inputImage->Update();
    
    std::cout << "OK." << std::endl;

    std::cout << "Warping image ... " << std::endl;
    ImageType::Pointer warpedImage;
    
    WarpImageFilterType::Pointer warpImageFilter = WarpImageFilterType::New();
    InterpolatorType::Pointer interpolator = InterpolatorType::New();
    warpImageFilter->SetInterpolator( interpolator );
    warpImageFilter->SetOutputSpacing( displacementField->GetSpacing() );
    warpImageFilter->SetOutputOrigin( displacementField->GetOrigin() );
    warpImageFilter->SetOutputDirection( displacementField->GetDirection() );
    warpImageFilter->SetDisplacementField( displacementField );
    warpImageFilter->SetInput( inputImage );
    
    warpedImage = warpImageFilter->GetOutput();
    warpedImage->Update();
    
	  std::cout << "Saving output image..." << std::endl;
    ImageWriterType::Pointer outputImageWriter = ImageWriterType::New();
	  outputImageWriter->SetInput( warpedImage );
    outputImageWriter->SetFileName( predictionMatrixFilename.c_str() );

  }

  // -------------------------------------------------------------

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "imiGenericMLRPrediction finished."          << std::endl;
  std::cout << "==========================================" << std::endl;

  return EXIT_SUCCESS;
  
} // end of main
