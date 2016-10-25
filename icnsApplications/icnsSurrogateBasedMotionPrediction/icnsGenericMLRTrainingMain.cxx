/** \file icnsGenericMLRTrainingMain.cxx
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

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "icnsMLRMotionPrediction.h"
#include "icnsSurrogateBasedMotionPredictionTypeDefinitions.h"

// Global typedefs:
typedef float                       MatrixValueType;
typedef vnl_vector<MatrixValueType> VnlVectorType;
typedef vnl_matrix<MatrixValueType> VnlMatrixType;

using namespace imi;

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "icnsGenericMLRTrainingMain ... \n";

  std::cout << "-V            Regressand observations in standard Matlab format (.mat as float and -v4 ? | mha)\n";
  std::cout << "-V_mean       Mean regressand observation output (.mat)\n";
  std::cout << "-Z            Regressor observations in standard Matlab format (.mat as double and -v4 )\n";
  std::cout << "-Z_mean       Mean regressor observation output (.mat)\n";
  std::cout << "-B            Filename of system matrix for saving (.mat).\n";
  std::cout << "-h            Print this help.\n";
}

// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  std::cout << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsGenericMLRTraining                    " << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;
  
  if( argc < 2 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }

  // Initializing parameters with default values:
  
  std::vector<std::string> regressorObservationsFilenames;
  bool regressorObservations_flag         = false;
  unsigned int noOfRegressorObservations  = 0;
  
  std::vector<std::string> regressandObservationsFilenames;
  bool regressandObservations_flag        = false;
  unsigned int nRegressandObservations = 0;
  
  std::string meanRegressorFilename;
  std::string meanRegressandFilename;
  std::string maskFilename;
  
  std::string outSystemMatrixFilename;

  bool majorFlagChange = false;

  // Running through cmd line params:

  std::cout << std::endl;
  for( int i = 1; i < argc; i++ )
  {
    if( strcmp( argv[i], "-h" ) == 0 )
    {
      PrintHelp();
      return -1;
    }

    if( majorFlagChange )
    {
      regressorObservations_flag = false;
      regressandObservations_flag = false;
      majorFlagChange = false;
    }

    // Regressor related params:
    if( strcmp( argv[i], "-Z" ) == 0 )
    {
      regressorObservations_flag = true;
      regressandObservations_flag = false;
      std::cout << "Reading filenames of regressor observations:" << std::endl;
      continue;
    }
    if( strcmp( argv[i], "-Z_mean" ) == 0 )
    {
      i++;
      meanRegressorFilename = argv[i];
      std::cout << "Mean regressor (Z_mean) filename:        " << meanRegressorFilename << std::endl;
      majorFlagChange = true;
      continue;
    }

    // Regressand related params:
    if( strcmp( argv[i], "-V" ) == 0 )
    {
      regressorObservations_flag = false;
      regressandObservations_flag = true;
      std::cout << "Reading filenames of regressand observations:" << std::endl;
      continue;
    }
    if( strcmp( argv[i], "-V_mean" ) == 0 )
    {
      i++;
      meanRegressandFilename = argv[i];
      std::cout << "Mean regressand (V_mean) filename:        " << meanRegressandFilename << std::endl;
      majorFlagChange = true;
      continue;
    }

    // Sytem matrix related params:
    if( strcmp( argv[i], "-B" ) == 0 )
    {
      i++;
      outSystemMatrixFilename = argv[i];
      std::cout << "System matrix will be stored at (B_out): " << outSystemMatrixFilename << std::endl;
      majorFlagChange = true;
      continue;
    }

    //Mask
    if( strcmp( argv[i], "-M" ) == 0 )
    {
      i++;
      maskFilename = argv[i];
      std::cout << "Mask filename:        " << maskFilename << std::endl;
      majorFlagChange = true;
      continue;
    }

    // Now reading filenames corresponding to variable length params:
    if( regressorObservations_flag )
    {
      regressorObservationsFilenames.push_back( argv[i] );
      noOfRegressorObservations++;
      std::cout << "  Regressor observation " << noOfRegressorObservations << ":  " << regressorObservationsFilenames[noOfRegressorObservations - 1] << std::endl;
      continue;
    }
    if( regressandObservations_flag )
    {
      regressandObservationsFilenames.push_back( argv[i] );
      nRegressandObservations++;
      std::cout << "  Regressand observation " << nRegressandObservations << ": " << regressandObservationsFilenames[nRegressandObservations - 1] << std::endl;
      continue;
    }

    // Anything not accounted for?
    std::cout << "Unknown param at " << i << ": " << argv[i] << std::endl;
  }

  // -------------------------------------------------------------
  // Plausibility checks:
  
  if( nRegressandObservations != noOfRegressorObservations )
  {
    std::cerr << "Numbers of regressand and regressor observations are not the same!" << std::endl;
    return EXIT_FAILURE;
  }

  icnsMLRMotionPrediction* predictionFilter = icnsMLRMotionPrediction::New();

  // -------------------------------------------------------------
  // REGRESSOR PART:
  // -------------------------------------------------------------
  // Current implementation: Centering inside prediction filter. Thus:
  // First step: Generating regressor matrix for training purposes:
  //
  // If -S is set and a series of filenames is given, they are stacked
  // together into a single regressor matrix.
  // -------------------------------------------------------------

  VnlMatrixType regressorTrainingMatrix;
  std::string fileTypePattern;

  if( noOfRegressorObservations == 0 )
  {
    std::cerr << " No regressor observations specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    fileTypePattern = icnsSurrogateBasedMotionPredictionHelpers::GetFileEnding( regressorObservationsFilenames[0] );

    if( (noOfRegressorObservations == 1) && (strcmp( fileTypePattern.c_str(), "mat" ) == 0) )
    {
      std::cerr << "Reading regressor training matrix not yet supported. Aborting computation." << std::endl;
      return EXIT_FAILURE;
      //ReadMatrixFromMatlabFile( regressorTrainingMatrix, regressorObservationsFilenames[0] );
    }
    else if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      std::cout << "Generating regressor training matrix from individual samples ... " << std::endl;

      VnlMatrixType currentRegressorObservationMatrix;
      for( unsigned int obsCounter = 0; obsCounter < noOfRegressorObservations; obsCounter++ )
      {
        std::cout << "  Reading observation " << regressorObservationsFilenames[obsCounter] << " ..." << std::endl;
        vcl_ifstream fid;
        fid.open( regressorObservationsFilenames[obsCounter].c_str() );
        vnl_matlab_readhdr matlabHeader( fid );
        
        // Depending on type, use the appropriate reading function:
        if( matlabHeader.is_single() )
        {
          icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileRealType(currentRegressorObservationMatrix,regressorObservationsFilenames[obsCounter]);
        }
        else
        {
          icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble(currentRegressorObservationMatrix,regressorObservationsFilenames[obsCounter]);
        }
        
        if( obsCounter == 0 )
        {
          regressorTrainingMatrix.set_size( currentRegressorObservationMatrix.rows(), regressorObservationsFilenames.size() );
          regressorTrainingMatrix.fill( 0.0 );
        }
        regressorTrainingMatrix.set_column( obsCounter, currentRegressorObservationMatrix.get_column( 0 ) );
      }
    }
    else
    {
      std::cout << "MLR not implemented for regressor file ending: " << fileTypePattern << std::endl;
      return EXIT_FAILURE;
    }
  }

  // -------------------------------------------------------------
  // REGRESSAND PART:
  // -------------------------------------------------------------
  // STANDARD USE CASE: Regressand observations are given as .mat
  // files.
  // -------------------------------------------------------------

  std::cout << "Generating regressand training matrix ... " << std::endl;
  VnlMatrixType regressandTrainingMatrix;

  if( nRegressandObservations == 0 )
  {
    std::cout << " No regressand observations specified. Aborting computation." << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    fileTypePattern = icnsSurrogateBasedMotionPredictionHelpers::GetFileEnding( regressandObservationsFilenames[0] );
    if( strcmp( fileTypePattern.c_str(), "mha" ) == 0 )
    {
      std::cout << "Generating float regressand training matrix from individual motion fields ... " << std::endl;
      icnsSurrogateBasedMotionPredictionHelpers::GenerateMatrixFromVectorFields( regressandTrainingMatrix, regressandObservationsFilenames, maskFilename );
    }
    
    else if( (nRegressandObservations == 1) && (strcmp( fileTypePattern.c_str(), "mat" ) == 0) )
    {
      std::cout << "Reading regressand training matrix not yet supported. Aborting computation." << std::endl;
      return EXIT_FAILURE;
    }
    
    else if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      std::cout << "Generating regressand training matrix from individual mat samples ... " << std::endl;

      VnlMatrixType currentRegressandObservationMatrix;
      for( unsigned int obsCounter = 0; obsCounter < nRegressandObservations; obsCounter++ )
      {
        std::cout << "  Reading observation " << regressandObservationsFilenames[obsCounter] << " ... " << std::endl;
        vcl_ifstream fid;
        fid.open( regressandObservationsFilenames[obsCounter].c_str() );
        vnl_matlab_readhdr matlabHeader( fid );
        
        // Depending on type, use the appropriate reading function:
        if( matlabHeader.is_single() )
        {
          icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileRealType( currentRegressandObservationMatrix, regressandObservationsFilenames[obsCounter]);
        }
        else
        {
          icnsSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble( currentRegressandObservationMatrix, regressandObservationsFilenames[obsCounter]);
        }
        
        if( obsCounter == 0 )
        {
          regressandTrainingMatrix.set_size( currentRegressandObservationMatrix.rows(), regressandObservationsFilenames.size() );
          regressandTrainingMatrix.fill( 0.0 );
        }
        regressandTrainingMatrix.set_column( obsCounter, currentRegressandObservationMatrix.get_column( 0 ) );
      }
    }
    else
    {
      std::cout << "MLR not implemented for regressor file ending: " << fileTypePattern << std::endl;
      return EXIT_FAILURE;
    }
  }

  // -------------------------------------------------------------
  // MLR TRAINING PART:
  // -------------------------------------------------------------

  std::cout << "  Training OLS system matrix ... " << std::endl;

  predictionFilter->SetRegressorTrainingMatrix( regressorTrainingMatrix );
  predictionFilter->SetRegressandTrainingMatrix( regressandTrainingMatrix );
  //predictionFilter->SetMulticollinearityCheck( false );
  predictionFilter->SetMulticollinearityCheck( true );
  predictionFilter->TrainLSEstimator();

  // -------------------------------------------------------------
  // If B is set, then write the system matrix to file.
  // Assumption: Matlab format.
  // -------------------------------------------------------------

  if( !outSystemMatrixFilename.empty() )
  {
    // Converting imiRealType-matrix to double (= matlab standard format):
    std::cout << "  Writing system matrix to file ... " << std::endl;

    VnlMatrixType lsEstimator;
    predictionFilter->GetTrainedEstimator( lsEstimator );

    vnl_matlab_filewrite outSystemMatrixFile( outSystemMatrixFilename.c_str() );
    outSystemMatrixFile.write( lsEstimator, "SystemMatrix" );
  }

  // -------------------------------------------------------------
  // If Z_mean is set, then write Z_mean to file.
  // Assumption: Matlab format.
  // -------------------------------------------------------------

  if( !meanRegressorFilename.empty() )
  {
    std::cout << "  Writing mean regressor to file ... " << std::endl;

    VnlMatrixType meanRegressor;
    predictionFilter->GetMeanRegressor( meanRegressor );

    vnl_matlab_filewrite outMeanRegressorFile( meanRegressorFilename.c_str() );
    outMeanRegressorFile.write( meanRegressor, "MeanRegressor" );
  }

  // -------------------------------------------------------------
  // If V_mean is set, then write V_mean to file.
  // Assumption: Matlab format.
  // -------------------------------------------------------------

  if( !meanRegressandFilename.empty() )
  {
    std::cout << "  Writing mean regressand to file ... "  << std::endl;

    VnlMatrixType meanRegressand;
    predictionFilter->GetMeanRegressand( meanRegressand );

    vnl_matlab_filewrite outMeanRegressandFile( meanRegressandFilename.c_str() );
    outMeanRegressandFile.write( meanRegressand, "MeanRegressand" );
  }

  //predictionFilter->Delete(); // TODO

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "icnsGenericMLRTraining finished."           << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
  
} // end of main
