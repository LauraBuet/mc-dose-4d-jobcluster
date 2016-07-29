/** \file imiGenericMLRPredictionMain.cxx
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
#include "itkImageFileReader.h"

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiMLRMotionPrediction.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"
#include "imiImageFunctions.h"
#include "imiImageWriter.h"
#include "imiImageReader.h"

#include "imiVelocityFieldExponential.h"
#include "imiVelocityFieldScalingAndSquaringExponential.h"
#include "imiVelocityFieldEulerStepExponential.h"

using namespace imi;


void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenericMLRPredictionMain ... \n\n";

  std::cout << "-V            Mean regressand observation in standard Matlab format (.mat) | List of regressand observations (.mat|.mha)\n";
  std::cout << "-Z            Mean regressor observation in standard Matlab format (.mat) | List of regressor observations (.mat|.mha)\n";
  std::cout << "-B            Filename of system matrix (.mat). If not specified, a new system matrix will be estimated.\n";
  std::cout << "-I            Input measurement (.mat).\n";
  std::cout << "-O            Prediction output (.mat)|(.mha).\n";
  std::cout << "-M            Mask file to restrict the computations but, only usable with -O *.mha.";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{
  imiINFO( "==========================================" );
  imiINFO( "====    imiGenericMLRPredictionMain     ====" );
  imiINFO( "==========================================" );
  imiINFO( "Reading parameters ...\n" );

  if( argc < 2 )
  {
    PrintHelp();
    return 1;
  }

  imiObject::SetGlobalDebugLevel(8);


  // VARIABLES AND CALL PARAMS:

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
      return -1;
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
    if( strcmp( argv[i], "-V" ) == 0 )
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

  imiMLRMotionPrediction* predictionFilter = imiMLRMotionPrediction::New();

  /* ******************************************************************
   * Prediction step
   * Assuming regressor observation given as .mat.
   * Prediction will be saved unter specified name.
   * ******************************************************************/

  vcl_ifstream fid;

  VnlMatrixType meanRegressor;
  VnlMatrixType meanRegressand;
  VnlMatrixType systemMatrix;
  VnlMatrixType measurement;
  VnlMatrixType prediction;

  if( meanRegressorFilename.empty() )
  {
    imiERROR( " No mean regressor specified. Aborting computation." );
    return -1;
  }
  else
  {
    imiINFO( "  Reading mean regressor (Z_mean) ..." );
    fid.open( meanRegressorFilename.c_str() );
    vnl_matlab_read_or_die( fid, meanRegressor );
    fid.close();
  }

  if( meanRegressandFilename.empty() )
  {
    imiERROR( " No mean regressand specified. Aborting computation." );
    return -1;
  }
  else
  {
    imiINFO( "  Reading mean regressand (V_mean) ..." );
    fid.open( meanRegressandFilename.c_str() );
    vnl_matlab_read_or_die( fid, meanRegressand );
    fid.close();
  }

  if( inSystemMatrixFilename.empty() )
  {
    imiERROR( " No system matrix specified. Aborting computation." );
    return -1;
  }
  else
  {
    imiINFO( "  Reading system matrix (B) ..." );
    fid.open( inSystemMatrixFilename.c_str() );
    vnl_matlab_read_or_die( fid, systemMatrix );
    fid.close();
  }

  if( inSystemMatrixFilename.empty() )
  {
    imiERROR( " No measurement specified. Aborting computation." );
    return -1;
  }
  else
  {
    imiINFO( "  Reading measurement (Z_hat) ..." );
    fid.open( inputMatrixFilename.c_str() );
    imiSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble( measurement, inputMatrixFilename);
    fid.close();
  }

  predictionFilter->SetMeanRegressor(meanRegressor);
  predictionFilter->SetMeanRegressand(meanRegressand);
  predictionFilter->SetTrainedEstimator(systemMatrix);
  predictionFilter->PredictOutput(measurement,prediction);

  if( !predictionMatrixFilename.empty() && inputImageFilename.empty() )
  {


    imiINFO( "Saving prediction output ... " );
    imiSurrogateBasedMotionPredictionHelpers::SaveVnlVectorAsVectorField( prediction, predictionMatrixFilename, templateFilename, maskFilename );

  } else if (!inputImageFilename.empty())
  {
    DisplacementFieldPointerType outputField;
    imiSurrogateBasedMotionPredictionHelpers::GenerateVectorFieldFromVnlVector( prediction, outputField, templateFilename, maskFilename );

    imiINFO("Allocate displacement field ... ");

    DisplacementFieldPointerType displacementField;

    displacementField = DisplacementFieldType::New();
    displacementField->CopyInformation( outputField );
    displacementField->SetBufferedRegion( outputField->GetBufferedRegion() );
    displacementField->SetRequestedRegion( outputField->GetRequestedRegion() );
    displacementField->Allocate();

    imiVelocityFieldExponential *velocityFieldExponentFilter;

    velocityFieldExponentFilter = imiVelocityFieldScalingAndSquaringExponential::New();
    velocityFieldExponentFilter->SetAccuracy( 4 );

    bool ret=false;

    if (inverseField)
    {
      imiINFO("Calculate inverse displacement field ...");
      ret = velocityFieldExponentFilter->GetInverseDisplacementFromVelocityField(outputField, displacementField);
    } else {
      imiINFO("Calculate displacement field ...");
      ret = velocityFieldExponentFilter->GetDisplacementFromVelocityField(outputField, displacementField);
    }

    if(!ret)
    {
      imiERROR("Calculating exponential failed!");
      return 1;
    }

    imiINFO("                             Done." );

  imiINFO( "Loading input image..." );
  imiImageReader *inputImageReader = imiImageReader::New();
  if( !inputImageReader->LoadImageData( inputImageFilename.c_str() ) )
  {
    imiERROR( "Cannot load input image!" );
    return -1;
  }


    ImageType::Pointer inputImage=inputImageReader->GetITKImage();
    ImageType::Pointer warpedImage;
    imiINFO( "Warping image ... " );

    if( !imiImageFunctions::WarpImageWithFactor(
        outputField, inputImage, warpedImage, 1, 0, 1 ) ) //0=NN 1=LIN
    {
      imiERROR( "Image warping failed!" );
      return false;
    }

	  imiINFO( "Saving output image..." );
	  imiImageWriter *outputImageWriter = imiImageWriter::New();
	  outputImageWriter->SetITKImageData( warpedImage );
	  if( !outputImageWriter->SaveImageData( predictionMatrixFilename.c_str() ) )
	  {
	    imiERROR( "Cannot save output image!" );
	    return false;
	  }

  }



  /* *****************************************************************/

  predictionFilter->Delete();

  imiINFO( "\n------------------------------------------" );
  imiINFO( "imiGenericMLRPrediction finished." );
  imiINFO( "==========================================\n" );

  return 0;
} // end of main
