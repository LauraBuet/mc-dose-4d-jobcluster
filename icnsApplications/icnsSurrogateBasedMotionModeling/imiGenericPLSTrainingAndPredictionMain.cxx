/** \file imiGenericPLSTrainingAndPredictionMain.cxx
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
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiPLSMotionPrediction.h"
#include "imiMotionPredictionLOOCV.h"

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenericPLSTrainingAndPrediction ... \n\n";

  std::cout << "-V            List of regressand observations (.mat|.mha) | \n";
  std::cout << "-Z            List of regressor observations (.mat|.mha)\n";
  std::cout << "-I            Input measurement (.mat).\n";
  std::cout << "-O            Prediction output (.mat)|(.mha).\n";
  std::cout << "-M            Mask file to restrict the computations, but only usable with -O *.mha.\n";
  std::cout << "-c            Number of componentes used for prediction.\n";
  std::cout << "-l            Perform LOOCV to determine optimal number of components used for prediction.\n";
  std::cout << "-m            1(default)|2 1: Only the new regressor basis is used for prediction. 2: Use new regressors and regressand bases for prediction.\n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{
  imiINFO( "===================================================" );
  imiINFO( "====    imiGenericPLSTrainingAndPrediction     ====" );
  imiINFO( "===================================================" );
  imiINFO( "Reading parameters ...\n" );

  if( argc < 2 )
  {
    PrintHelp();
    return 1;
  }

  imiObject::SetGlobalDebugLevel( 5 );

  imiTime::imiInitTime();

  // TYPEDEFS:

  typedef imiRealType MatrixValueType;
  typedef vnl_vector<MatrixValueType> VnlVectorType;
  typedef vnl_matrix<MatrixValueType> VnlMatrixType;

  // VARIABLES AND CALL PARAMS:

  std::vector<std::string> regressorObservationsFilenames;
  bool regressorObservations_flag = false;
  unsigned int noOfRegressorObservations = 0;

  std::vector<std::string> regressandObservationsFilenames;
  bool regressandObservations_flag = false;
  unsigned int noOfRegressandObservations = 0;

  std::string maskFilename;

  std::vector<std::string> inputFilenames;
  bool inputObservations_flag = false;
  unsigned int noOfInputObservations = 0;

  std::vector<std::string> predictionFilenames;
  bool outputFilenames_flag = false;
  unsigned int noOfOutputFilenames = 0;

  bool majorFlagChange = false;

  unsigned int noOfComponents = 1;
  bool loocv=false;

  bool useBothBases=false;

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
    // Regressand related params:
    if( strcmp( argv[i], "-V" ) == 0 )
    {
      regressorObservations_flag = false;
      regressandObservations_flag = true;
      std::cout << "Reading filenames of regressand observations:" << std::endl;
      continue;
    }

    //Input measurement params:
    if( strcmp( argv[i], "-I" ) == 0 )
    {
      regressorObservations_flag = false;
      regressandObservations_flag = false;
      inputObservations_flag = true;
      outputFilenames_flag = false;
      std::cout << "Reading filenames of input measurements:" << std::endl;
      continue;
    }

    if( strcmp( argv[i], "-c" ) == 0 )
    {
      i++;
      noOfComponents = atoi(argv[i]);
      std::cout << "Number of components used for prediction:        " << noOfComponents << std::endl;
      majorFlagChange = true;
      continue;
    }

    if( strcmp( argv[i], "-l" ) == 0 )
    {
      loocv=true;
      std::cout << "Performing LOOCV to determine optimal number of components."<< std::endl;
      majorFlagChange = true;
      continue;
    }

    if( strcmp( argv[i], "-m" ) == 0 )
    {
      i++;
      if (atoi(argv[i])>1)
      {
        useBothBases=true;
        std::cout << "Using both bases for prediction."<< std::endl;
      } else
      {
        useBothBases=false;
        std::cout << "Only using the new regressor basis for prediction."<< std::endl;
      }
      majorFlagChange = true;
      continue;
    }

    //Prediction:
    if( strcmp( argv[i], "-O" ) == 0 )
    {
      regressorObservations_flag = false;
      regressandObservations_flag = false;
      inputObservations_flag = false;
      outputFilenames_flag = true;
      std::cout << "Reading filenames of prediction outputs:" << std::endl;
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
    if( inputObservations_flag )
    {
      inputFilenames.push_back( argv[i] );
      noOfInputObservations++;
      std::cout << "  Input measurement " << noOfInputObservations << ": " << inputFilenames[noOfInputObservations - 1] << std::endl;
      continue;
    }
    if( outputFilenames_flag )
    {
      predictionFilenames.push_back( argv[i] );
      noOfOutputFilenames++;
      std::cout << "  Output file " << noOfOutputFilenames << ": " << predictionFilenames[noOfOutputFilenames - 1] << std::endl;
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

  imiPLSMotionPrediction* predictionFilter = imiPLSMotionPrediction::New();
  predictionFilter->SetUseBothBases(useBothBases);


  /****************************************************************
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
        vnl_matlab_read_or_die( fid, currentRegressorObservationMatrix );
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
      imiERROR( "PLS not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }
  }

  /****************************************************************
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
      imiSurrogateBasedMotionPredictionHelpers::GenerateMatrixFromVectorFields(regressandTrainingMatrix,regressandObservationsFilenames,maskFilename);

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
        vnl_matlab_read_or_die( fid, currentRegressandObservationMatrix );
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
      imiERROR( "PLS not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }
  }

  /****************************************************************
   *    PLS TRAINING PART:
   ****************************************************************/

  if (loocv)
  {
    imiMotionPredictionLOOCV* loocvTest=imiMotionPredictionLOOCV::New();
    loocvTest->SetRegressandTrainingMatrix(regressandTrainingMatrix);
    loocvTest->SetRegressorTrainingMatrix(regressorTrainingMatrix);
    noOfComponents=loocvTest->PerformLOOCVPLS(1,std::min<unsigned int>(regressorTrainingMatrix.cols()-2,std::min<unsigned int>(regressorTrainingMatrix.rows(),regressandTrainingMatrix.rows())));
    imiINFO( "LOOCV: "<<noOfComponents<<" selected!");

  } else
  {

    if(noOfComponents > std::min<unsigned int>(regressorTrainingMatrix.rows(),regressandTrainingMatrix.rows()))
    {
      noOfComponents=std::min<unsigned int>(regressorTrainingMatrix.rows(),regressandTrainingMatrix.rows());
      imiWARNING("The number of eigenvectors specified exceeds the number of min(#regressor variables,#regressand variables). Setting number of eigenvectors to "<<noOfComponents<<" (min(#regressor variables,#regressand variables))!");
    }
  }

  predictionFilter->SetNumberOfComponents(noOfComponents);
  predictionFilter->SetRegressorTrainingMatrix( regressorTrainingMatrix );
  predictionFilter->SetRegressandTrainingMatrix( regressandTrainingMatrix );
  std::cout << std::endl;
  imiINFO( "Training PLS system matrix ... " );

  predictionFilter->TrainEstimator();


  /****************************************************************
   *    PLS PREDICTION PART:
   ****************************************************************/

  for( unsigned int i = 0; i < noOfInputObservations; i++ )
  {
  vcl_ifstream fid;

  VnlMatrixType measurement;
  VnlMatrixType prediction;

  imiINFO( "  Reading measurement (Z_hat) "<<i+1<<" ..." );
  fid.open( inputFilenames[i].c_str() );
  vnl_matlab_read_or_die( fid, measurement );
  fid.close();

  predictionFilter->PredictOutput(measurement,prediction);

  imiINFO( "Saving prediction output "<<i+1<<" ... " );

    fileTypePattern = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( predictionFilenames[i]  );

    if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      vnl_matlab_filewrite predictionFile( predictionFilenames[i].c_str() );
      predictionFile.write( prediction, "Prediction" );
    } else if ( strcmp( fileTypePattern.c_str(), "mha" ) == 0 )
    {
      imiSurrogateBasedMotionPredictionHelpers::SaveVnlVectorAsVectorField( prediction, predictionFilenames[i], regressandObservationsFilenames[0], maskFilename );
    } else
    {
      imiERROR( "Unsupported file type: " << fileTypePattern );
            return -1;
    }

  }




  predictionFilter->Delete();

  imiINFO( "\n------------------------------------------" );
  imiINFO( "imiGenericPLSTraining finished." );
  imiINFO( "==========================================\n" );

  return 0;
} // end of main
