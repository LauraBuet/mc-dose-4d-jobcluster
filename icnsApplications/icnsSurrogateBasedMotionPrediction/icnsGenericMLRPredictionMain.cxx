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
  std::cout << "-T            IN: Template field to allow for output conversion from matlab to mha.\n";
  
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
  std::string systemMatrixFilename;
  std::string inputMatrixFilename;
  
  std::string maskFilename;
  std::string templateFieldFilename;
  
  std::string predictionMatrixFilename;
  std::string predictionFieldFilename;
  
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
    if( strcmp( argv[i], "-Z_mean" ) == 0 )
    {
      i++;
      meanRegressorFilename = argv[i];
      std::cout << "Mean regressor (Z_mean) filename:    " << meanRegressorFilename << std::endl;
      continue;
    }

    // Regressand related params:
    if( strcmp( argv[i], "-V_mean" ) == 0 )
    {
      i++;
      meanRegressandFilename = argv[i];
      std::cout << "Mean regressand (V_mean) filename:   " << meanRegressandFilename << std::endl;
      continue;
    }

    // Sytem matrix related params:
    if( strcmp( argv[i], "-B" ) == 0 )
    {
      i++;
      systemMatrixFilename = argv[i];
      std::cout << "System matrix (B) filename:          " << systemMatrixFilename << std::endl;
      continue;
    }

    // Input matrix related params:
    if( strcmp( argv[i], "-Z_hat" ) == 0 )
    {
      i++;
      inputMatrixFilename = argv[i];
      std::cout << "Input measurement (Z_hat) filename:  " << inputMatrixFilename << std::endl;
      continue;
    }

    // Prediction matrix related params:
    if( strcmp( argv[i], "-V_hat" ) == 0 )
    {
      i++;
      predictionMatrixFilename = argv[i];
      std::cout << "Prediction (V_hat) filename:         " << predictionMatrixFilename << std::endl;
      continue;
    }
    
    // Prediction matrix related params:
    if( strcmp( argv[i], "-v_hat" ) == 0 )
    {
      i++;
      predictionFieldFilename = argv[i];
      std::cout << "Prediction (v_hat) filename:         " << predictionFieldFilename << std::endl;
      continue;
    }

    // Mask for writing appropriate field data:
    if( strcmp( argv[i], "-M" ) == 0 )
    {
      i++;
      maskFilename = argv[i];
      std::cout << "Mask filename:                       " << maskFilename << std::endl;
      continue;
    }

    // Template field for conversion to field.
    if( strcmp( argv[i], "-T" ) == 0 )
    {
      i++;
      templateFieldFilename = argv[i];
      std::cout << "Template field filename:             " << templateFieldFilename << std::endl;
      continue;
    }

    // Anything not accounted for?
    std::cout << "Unknown param at " << i << ": " << argv[i] << std::endl;
  }

  // -------------------------------------------------------------
  // Loading input data:
  // -------------------------------------------------------------

  icnsMLRMotionPrediction* predictionFilter = icnsMLRMotionPrediction::New();
  
  vcl_ifstream fid;

  VnlMatrixType meanRegressor;   // Z_mean
  VnlMatrixType meanRegressand;  // V_mean
  VnlMatrixType systemMatrix;    // B
  VnlMatrixType measurement;     // Z_hat
  VnlMatrixType prediction;      // V_hat

  // Reading mean regressor:
  if( meanRegressorFilename.empty() )
  {
    std::cerr << "  No mean regressor specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "  Reading mean regressor (Z_mean) ..." << std::endl;
    fid.open( meanRegressorFilename.c_str() );
    vnl_matlab_readhdr matlabHeader( fid );
    
    // Depending on type, use the appropriate reading function:
    if( matlabHeader.is_single() )
    {
      icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileRealType( meanRegressor, meanRegressorFilename.c_str() );
    }
    else
    {
      icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble( meanRegressor, meanRegressorFilename.c_str() );
    }
    
    fid.close();
  }
  
  // Reading mean regressor:
  if( meanRegressandFilename.empty() )
  {
    std::cerr << "  No mean regressand specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "  Reading mean regressand (V_mean) ..." << std::endl;
    fid.open( meanRegressandFilename.c_str() );
    vnl_matlab_readhdr matlabHeader( fid );
    
    // Depending on type, use the appropriate reading function:
    if( matlabHeader.is_single() )
    {
      icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileRealType( meanRegressand, meanRegressandFilename.c_str() );
    }
    else
    {
      icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble( meanRegressand, meanRegressandFilename.c_str() );
    }
    
    fid.close();
  }

  // Reading system matrix:
  if( systemMatrixFilename.empty() )
  {
    std::cerr << "  No system matrix specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    std::cout << "  Reading system matrix (B) ..." << std::endl;
    fid.open( systemMatrixFilename.c_str() );
    vnl_matlab_readhdr matlabHeader( fid );
    
    // Depending on type, use the appropriate reading function:
    if( matlabHeader.is_single() )
    {
      icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileRealType( systemMatrix, systemMatrixFilename.c_str() );
    }
    else
    {
      icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble( systemMatrix, systemMatrixFilename.c_str() );
    }
    
    fid.close();
  }

  // Reading input measurement:
  if( inputMatrixFilename.empty() )
  {
    std::cerr << "  No measurement specified. Aborting computation." << std::endl;
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
  
  // If V_hat has been set as input param, write prediction as mat-file:
  if( !predictionMatrixFilename.empty() )
  {
    std::cout << "  Saving prediction output as mat-file ... " << std::flush;
    
    vnl_matlab_filewrite outPredictionFile( predictionMatrixFilename.c_str() );
    outPredictionFile.write( prediction, "Prediction" );
    
    std::cout << "OK." << std::endl;
  }
  
  // If v_hat has been set as input param, write prediction as mha-field.
  // This, however, requires a templatefield to be given.
  if( !predictionFieldFilename.empty() && !templateFieldFilename.empty() )
  {
    std::cout << "  Saving prediction output as field data (.mha) ... " << std::flush;
    icnsSurrogateBasedMotionPredictionHelpers::SaveVnlVectorAsVectorField( prediction, predictionFieldFilename, templateFieldFilename, maskFilename );
    std::cout << "OK." << std::endl;
  }
 
  // -------------------------------------------------------------

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "icnsGenericMLRPrediction finished."         << std::endl;
  std::cout << "==========================================" << std::endl;

  return EXIT_SUCCESS;
  
} // end of main
