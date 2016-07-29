/** \file imiGenericMLRTrainingMain.cxx
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck; distributed to IPMI at ICNS at UKE
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

//Typedefs

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenericMLRTrainingMain ... \n\n";

  std::cout << "-V            Regressand observations in standard Matlab format (.mat)\n";
  std::cout << "-V_mean       Mean regressand observation output (.mat)\n";
  std::cout << "-Z            Regressor observations in standard Matlab format (.mat)\n";
  std::cout << "-Z_mean       Mean regressor observation output (.mat)\n";
  std::cout << "-B            Filename of system matrix for saving (.mat).\n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{
  imiINFO( "==========================================" );
  imiINFO( "====    imiGenericMLRTrainingMain     ====" );
  imiINFO( "==========================================" );
  imiINFO( "Reading parameters ...\n" );

  if( argc < 2 )
  {
    PrintHelp();
    return 1;
  }

  // TYPEDEFS:

  typedef float MatrixValueType;
  typedef vnl_vector<MatrixValueType> VnlVectorType;
  typedef vnl_matrix<MatrixValueType> VnlMatrixType;

  // VARIABLES AND CALL PARAMS:

  std::vector<std::string> regressorObservationsFilenames;
  bool regressorObservations_flag = false;
  unsigned int noOfRegressorObservations = 0;
  std::string meanRegressorFilename;

  std::string maskFilename;

  std::vector<std::string> regressandObservationsFilenames;
  bool regressandObservations_flag = false;
  unsigned int noOfRegressandObservations = 0;
  std::string meanRegressandFilename;

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
      noOfRegressandObservations++;
      std::cout << "  Regressand observation " << noOfRegressandObservations << ": " << regressandObservationsFilenames[noOfRegressandObservations - 1] << std::endl;
      continue;
    }

    // Anything not accounted for?
    std::cout << "Unknown param at " << i << ": " << argv[i] << std::endl;
  }

  if( noOfRegressandObservations != noOfRegressorObservations )
  {
    imiERROR( "Numbers of regressand and regressor observations are not the same!" );
    return -1;
  }

  imiMLRMotionPrediction* predictionFilter = imiMLRMotionPrediction::New();

  /* ***************************************************************
   *    REGRESSOR PART:
   *****************************************************************
   * Current implementation: Centering inside prediction filter. Thus:
   * First step: Generating regressor matrix for training purposes:
   *
   * If -S is set and a series of filenames is given, they are stacked
   * together into a single regressor matrix.
   *****************************************************************/

  VnlMatrixType regressorTrainingMatrix;
  std::string fileTypePattern;

  if( noOfRegressorObservations == 0 )
  {
    imiERROR( " No regressor observations specified. Aborting computation." );
    return -1;
  }
  else
  {
    fileTypePattern = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( regressorObservationsFilenames[0] );

    if( (noOfRegressorObservations == 1) && (strcmp( fileTypePattern.c_str(), "mat" ) == 0) )
    {
      imiINFO( "Reading regressor training matrix not yet supported. Aborting computation." );
      return -1;
      //ReadMatrixFromMatlabFile( regressorTrainingMatrix, regressorObservationsFilenames[0] );
    }
    else if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      imiINFO( "Generating regressor training matrix from individual samples ... " );

      VnlMatrixType currentRegressorObservationMatrix;
      for( unsigned int obsCounter = 0; obsCounter < noOfRegressorObservations; obsCounter++ )
      {
        imiINFO( "  Reading observation " << regressorObservationsFilenames[obsCounter] << " ..." );
        vcl_ifstream fid;
        fid.open( regressorObservationsFilenames[obsCounter].c_str() );
        //vnl_matlab_read_or_die( fid, currentRegressorObservationMatrix );
        imiSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble(currentRegressorObservationMatrix,regressorObservationsFilenames[obsCounter]);
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
      imiERROR( "MLR not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }
  }

  /* ***************************************************************
   *    REGRESSAND PART:
   *****************************************************************
   * STANDARD USE CASE: Regressand observations are given as .mat
   * files.
   * ***************************************************************/

  imiINFO( "Generating regressand training matrix ... " );
  VnlMatrixType regressandTrainingMatrix;

  if( noOfRegressandObservations == 0 )
  {
    imiERROR( " No regressand observations specified. Aborting computation." );
    return -1;
  }
  else
  {
    fileTypePattern = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( regressandObservationsFilenames[0] );
    if( strcmp( fileTypePattern.c_str(), "mha" ) == 0 )
        {
          imiINFO( "Generating regressand training matrix from individual motion fields ... " );
          imiSurrogateBasedMotionPredictionHelpers::GenerateMatrixFromVectorFields( regressandTrainingMatrix, regressandObservationsFilenames, maskFilename );

        }
    else if( (noOfRegressandObservations == 1) && (strcmp( fileTypePattern.c_str(), "mat" ) == 0) )
    {
      imiINFO( "Reading regressand training matrix not yet supported. Aborting computation." );
      return -1;
      //ReadMatrixFromMatlabFile( regressorTrainingMatrix, regressorObservationsFilenames[0] );
    }
    else if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      imiINFO( "Generating regressand training matrix from individual samples ... " );

      VnlMatrixType currentRegressandObservationMatrix;
      for( unsigned int obsCounter = 0; obsCounter < noOfRegressandObservations; obsCounter++ )
      {
        imiINFO( "  Reading observation " << regressandObservationsFilenames[obsCounter] << " ..." );
        vcl_ifstream fid;
        fid.open( regressandObservationsFilenames[obsCounter].c_str() );
        imiSurrogateBasedMotionPredictionHelpers::ReadMatrixFromMatlabFileDouble( currentRegressandObservationMatrix, regressandObservationsFilenames[obsCounter]);
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
      imiERROR( "MLR not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }
  }

  /* ***************************************************************
   *    MLR TRAINING PART:
   ****************************************************************/

  std::cout << std::endl;
  imiINFO( "Training OLS system matrix ... " );

  predictionFilter->SetRegressorTrainingMatrix( regressorTrainingMatrix );
  predictionFilter->SetRegressandTrainingMatrix( regressandTrainingMatrix );
  predictionFilter->SetMulticollinearityCheck( false );
  predictionFilter->TrainLSEstimator();

  /* ***************************************************************
   * If B is set, then write the system matrix to file.
   * Assumption: Matlab format.
   ****************************************************************/

  if( !outSystemMatrixFilename.empty() )
  {
    // Converting imiRealType-matrix to double (= matlab standard format):
    imiINFO( "  Writing system matrix to file ... " );

    VnlMatrixType lsEstimator;
    predictionFilter->GetTrainedEstimator( lsEstimator );

    vnl_matlab_filewrite outSystemMatrixFile( outSystemMatrixFilename.c_str() );
    outSystemMatrixFile.write( lsEstimator, "SystemMatrix" );
  }

  /* ***************************************************************
   * If Z_mean is set, then write Z_mean to file.
   * Assumption: Matlab format.
   ****************************************************************/

  if( !meanRegressorFilename.empty() )
  {
    imiINFO( "  Writing mean regressor to file ... " );

    VnlMatrixType meanRegressor;
    predictionFilter->GetMeanRegressor( meanRegressor );

    vnl_matlab_filewrite outMeanRegressorFile( meanRegressorFilename.c_str() );
    outMeanRegressorFile.write( meanRegressor, "MeanRegressor" );
  }

  /* ***************************************************************
   * If V_mean is set, then write V_mean to file.
   * Assumption: Matlab format.
   ****************************************************************/

  if( !meanRegressandFilename.empty() )
  {
    imiINFO( "  Writing mean regressand to file ... " );

    VnlMatrixType meanRegressand;
    predictionFilter->GetMeanRegressand( meanRegressand );

    vnl_matlab_filewrite outMeanRegressandFile( meanRegressandFilename.c_str() );
    outMeanRegressandFile.write( meanRegressand, "MeanRegressand" );
  }

  predictionFilter->Delete();

  imiINFO( "\n------------------------------------------" );
  imiINFO( "imiGenericMLRTraining finished." );
  imiINFO( "==========================================\n" );

  return 0;
} // end of main

