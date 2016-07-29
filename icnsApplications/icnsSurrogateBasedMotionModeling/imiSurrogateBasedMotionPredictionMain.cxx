/** \file imiSurrogateBasedMotionPredictionMain.cxx
 *
 *  \b Initial \b Author: werner \n\n
 *  \b Copyright (C) 2011 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/


// System includes:
#include <iostream>
#include <string>
#include <vector>
#include <complex>
extern "C"
{
#include "getopt.h"
}

// ITK includes:
#include "itkCastImageFilter.h"
#include "itkImageMaskSpatialObject.h"
#include "itkImageFileWriter.h"

// VNL includes:
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiImageFieldReader.h"
#include "imiImageReader.h"
#include "imiImageWriter.h"
#include "imiMLRMotionPrediction.h"


using namespace imi;


void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiSurrogateBasedMotionPrediction ... \n\n";

  std::cout << "-V            Vector fields (in .mha format) or arbitrary regressand observations\n";
  std::cout << "              in standard Matlab format (.mat).\n";
  std::cout << "-V_mean       Mean vector field or regressand observation (.mha / .mat format).\n";
  std::cout << "              NB: If B is specified, V_mean will be loaded and used for\n";
  std::cout << "              prediction. Otherwise, V_mean will be saved.\n";
  std::cout << "-V_predicted  Filenames for the predicted regressand values.\n";
  std::cout << "-V_out        Store regressand observations under the given filenames.\n";
  std::cout << "-S            Regressor samples in in standard Matlab format (.mat).\n";
  std::cout << "-S_mean       Mean regressor observation (.mat format).\n";
  std::cout << "              NB: If B is specified, S_mean will be loaded and used for\n";
  std::cout << "              prediction. Otherwise, S_mean will be saved.\n";
  std::cout << "-S_diaphragm  If this flag is set, diaphragm motion will be extracted from the\n";
  std::cout << "              vector fields specified by -V.\n";
  std::cout << "-S_measured   Regressor measurement(s) used for prediction purposes.\n";
  std::cout << "-S_out        Store regressor observations under the given filenames.\n";
  std::cout << "              will be saved under the specified filename.\n";
  std::cout << "-B            Filename of system matrix for loading. If specified, no new system\n";
  std::cout << "              is computed.\n";
  std::cout << "-B_out        Filename of system matrix for saving.\n";
  std::cout << "-M            Filename for image to mask MLR training.\n";
  std::cout << "-N            Filename of diaphragm image (currently: only for saving).\n";
  std::cout << "-n_diaphragm  Number of points on diaphragm (used only if -S_diaphragm is set; default = 1).\n";
  std::cout << "-r_diaphragm  Diaphragm points, inner circle radius (used only if -S_diaphragm is set; default = 5).\n";
  std::cout << "-h            Print this help.\n";
  std::cout << "\n";
}

bool ExtractDiaphragmMotionFromRegressandFields( vnl_matrix<imiRealType>& diaphragmMotionMatrix,
                                                 std::vector<std::string> fieldFilenames,
                                                 std::vector<std::string> outRegressorFilenames,
                                                 std::string leftLungFilename,
                                                 std::string RightLungFilename,
                                                 std::string diaphragmFilename,
                                                 unsigned int noOfPoints,
                                                 float innerCircleRadius );
bool GenerateRegressandMatrixFromVectorFields( vnl_matrix<imiRealType>& regressandMatrix,
                                               std::vector<std::string> fieldFilenames,
                                               std::string maskFilename );
std::string GetFileEnding( std::string filename );
bool ReadMatrixFromMatlabFile( vnl_matrix<imiRealType>& vnlMatrix, std::string filename );
bool SavePredictedVnlVectorAsVectorField( vnl_matrix<imiRealType>& predictedOutputAsVnlVector,
                                          std::string outFieldFilename,
                                          std::string refFieldname,
                                          std::string maskFilename );


int main( int argc, char *argv[] )
{
  imiINFO( "==========================================" );
  imiINFO( "=== imiSurrogateBasedMotionPrediction ====" );
  imiINFO( "==========================================" );
  imiINFO( "Reading parameters ...\n" );

  if( argc < 1 )
  {
    PrintHelp();
    return 1;
  }

  // TYPEDEFS:

  typedef imiRealType MatrixValueType;
  typedef vnl_vector<MatrixValueType> VnlVectorType;
  typedef vnl_matrix<MatrixValueType> VnlMatrixType;

  // VARIABLES AND CALL PARAMS:

  std::vector<std::string> regressorObservationsFilenames;
  bool regressorObservations_flag = false;
  unsigned int noOfRegressorObservations = 0;
  std::string meanRegressorFilename;
  std::vector<std::string> regressorMeasurementsFilenames;
  bool regressorMeasurements_flag = false;
  unsigned int noOfRegressorMeasurements = 0;
  std::vector<std::string> outGeneratedRegressorObservationsFilenames;
  bool outGeneratedRegressorObservationsFilenames_flag = false;
  unsigned int noOfOutGeneratedRegressorObservationsFilenames = 0;

  std::string maskFilename;

  std::string diaphragmFilename;
  std::string leftLungFilename;
  std::string rightLungFilename;
  unsigned int noOfDiaphragmPoints = 1;
  float radiusDiaphragmCircle = 5;

  std::vector<std::string> regressandObservationsFilenames;
  bool regressandObservations_flag = false;
  unsigned int noOfRegressandObservations = 0;
  std::string meanRegressandFilename;
  std::vector<std::string> outGeneratedRegressandObservationsFilenames;
  bool outGeneratedRegressandObservationsFilenames_flag = false;
  unsigned int noOfOutGeneratedRegressandObservationsFilenames = 0;

  std::vector<std::string> predictedOutputsFilenames;
  bool predictedOutputs_flag = false;
  unsigned int noOfPredictedOutputsFilenames = 0;

  std::string inSystemMatrixFilename;
  std::string outSystemMatrixFilename;

  bool extractDiaphragmMotionFromVectorfields = false;
  bool majorFlagChange = false;

  // Running through cmd line params:

  std::cout << std::endl;
  for (int i = 1; i < argc; i++)
  {
    if(strcmp(argv[i], "-h") == 0)
    {
      PrintHelp();
      return -1;
    }

    if( majorFlagChange )
    {
      regressorObservations_flag = false;
      regressorMeasurements_flag = false;
      regressandObservations_flag = false;
      predictedOutputs_flag = false;
      majorFlagChange = false;
      outGeneratedRegressorObservationsFilenames_flag = false;
      outGeneratedRegressandObservationsFilenames_flag = false;
    }

    // Regressor related params:
    if(strcmp(argv[i], "-S") == 0)
    {
      regressorObservations_flag = true;
      regressorMeasurements_flag = false;
      regressandObservations_flag = false;
      predictedOutputs_flag = false;
      outGeneratedRegressorObservationsFilenames_flag = false;
      outGeneratedRegressandObservationsFilenames_flag = false;
      std::cout << "Reading filenames of regressor observations:" << std::endl;
      continue;
    }
    if(strcmp(argv[i], "-S_mean") == 0)
    {
      i++;
      meanRegressorFilename = argv[i];
      std::cout << "Mean regressor (S_mean) filename:        " << meanRegressorFilename << std::endl;
      majorFlagChange = true;
      continue;
    }
    if(strcmp(argv[i], "-S_measured") == 0)
    {
      regressorObservations_flag = false;
      regressorMeasurements_flag = true;
      regressandObservations_flag = false;
      predictedOutputs_flag = false;
      outGeneratedRegressorObservationsFilenames_flag = false;
      outGeneratedRegressandObservationsFilenames_flag = false;
      std::cout << "Reading filenames of regressor measurements:" << std::endl;
      continue;
    }
    if(strcmp(argv[i], "-S_diaphragm") == 0)
    {
      extractDiaphragmMotionFromVectorfields = true;
      leftLungFilename = argv[++i];
      rightLungFilename = argv[++i];
      std::cout << "Regressor will be diaphragm motion extracted from regressand vector fields!" << std::endl;
      std::cout << "  Left lung file name:                   " << leftLungFilename << std::endl;
      std::cout << "  Right lung file name:                  " << rightLungFilename << std::endl;
      majorFlagChange = true;
      continue;
    }
    if(strcmp(argv[i], "-n_diaphragm") == 0)
    {
      i++;
      noOfDiaphragmPoints = atoi(argv[i]);
      std::cout << "Number of points on diaphragm:           " << noOfDiaphragmPoints << std::endl;
      majorFlagChange = true;
      continue;
    }
    if(strcmp(argv[i], "-r_diaphragm") == 0)
    {
      i++;
      radiusDiaphragmCircle = atof(argv[i]);
      std::cout << "Radius of inner circle on diaphragm:     " << radiusDiaphragmCircle << std::endl;
      majorFlagChange = true;
      continue;
    }
    if(strcmp(argv[i], "-S_out") == 0)
    {
      regressorObservations_flag = false;
      regressorMeasurements_flag = false;
      regressandObservations_flag = false;
      predictedOutputs_flag = false;
      outGeneratedRegressorObservationsFilenames_flag = true;
      outGeneratedRegressandObservationsFilenames_flag = false;
      std::cout << "Reading filenames to save generated regressors at:" << std::endl;
      continue;
    }

    // Regressand related params:
    if(strcmp(argv[i], "-V") == 0)
    {
      regressorObservations_flag = false;
      regressorMeasurements_flag = false;
      regressandObservations_flag = true;
      predictedOutputs_flag = false;
      outGeneratedRegressorObservationsFilenames_flag = false;
      outGeneratedRegressandObservationsFilenames_flag = false;
      std::cout << "Reading filenames of regressand observations:" << std::endl;
      continue;
    }
    if(strcmp(argv[i], "-V_mean") == 0)
    {
      i++;
      meanRegressandFilename = argv[i];
      std::cout << "Mean regressand (V_mean) filename:       " << meanRegressandFilename << std::endl;
      majorFlagChange = true;
      continue;
    }
    if(strcmp(argv[i], "-V_predicted") == 0)
    {
      regressorObservations_flag = false;
      regressorMeasurements_flag = false;
      regressandObservations_flag = false;
      predictedOutputs_flag = true;
      outGeneratedRegressorObservationsFilenames_flag = false;
      outGeneratedRegressandObservationsFilenames_flag = false;
      std::cout << "Reading filenames of outputs to be predicted:" << std::endl;
      continue;
    }
    if(strcmp(argv[i], "-V_out") == 0)
    {
      regressorObservations_flag = false;
      regressorMeasurements_flag = false;
      regressandObservations_flag = false;
      predictedOutputs_flag = false;
      outGeneratedRegressorObservationsFilenames_flag = false;
      outGeneratedRegressandObservationsFilenames_flag = true;
      std::cout << "Reading filenames to save generated regressands at:" << std::endl;
      continue;
    }

    // Sytem matrix related params:
    if(strcmp(argv[i], "-B") == 0)
    {
      i++;
      inSystemMatrixFilename = argv[i];
      std::cout << "System matrix to be loaded (B_in):       " << inSystemMatrixFilename << std::endl;
      majorFlagChange = true;
      continue;
    }
    if(strcmp(argv[i], "-B_out") == 0)
    {
      i++;
      outSystemMatrixFilename = argv[i];
      std::cout << "System matrix will be stored at (B_out): " << outSystemMatrixFilename << std::endl;
      majorFlagChange = true;
      continue;
    }
    if(strcmp(argv[i], "-M") == 0)
    {
      i++;
      maskFilename = argv[i];
      std::cout << "  Mask file name:                        " << maskFilename << std::endl;
      majorFlagChange = true;
      continue;
    }
    if(strcmp(argv[i], "-N") == 0)
    {
      i++;
      diaphragmFilename = argv[i];
      std::cout << "  Diaphragm image file name:             " << diaphragmFilename << std::endl;
      majorFlagChange = true;
      continue;
    }

    // Now reading filenames corresponding to variable length params:
    if( regressorObservations_flag )
    {
      regressorObservationsFilenames.push_back( argv[i] );
      noOfRegressorObservations++;
      std::cout << "  Regressor observation " << noOfRegressorObservations << ":  " << regressorObservationsFilenames[noOfRegressorObservations-1] << std::endl;
      continue;
    }
    if( regressorMeasurements_flag )
    {
      regressorMeasurementsFilenames.push_back( argv[i] );
      noOfRegressorMeasurements++;
      std::cout << "  Regressor measurement " << noOfRegressorMeasurements << ":  " << regressorMeasurementsFilenames[noOfRegressorMeasurements-1] << std::endl;
      continue;
    }
    if( regressandObservations_flag )
    {
      regressandObservationsFilenames.push_back( argv[i] );
      noOfRegressandObservations++;
      std::cout << "  Regressand observation " << noOfRegressandObservations << ": " << regressandObservationsFilenames[noOfRegressandObservations-1] << std::endl;
      continue;
    }
    if( predictedOutputs_flag )
    {
      predictedOutputsFilenames.push_back( argv[i] );
      noOfPredictedOutputsFilenames++;
      std::cout << "  Predicted output FN " << noOfPredictedOutputsFilenames << ":    " << predictedOutputsFilenames[noOfPredictedOutputsFilenames-1] << std::endl;
      continue;
    }
    if( outGeneratedRegressorObservationsFilenames_flag )
    {
      outGeneratedRegressorObservationsFilenames.push_back( argv[i] );
      noOfOutGeneratedRegressorObservationsFilenames++;
      std::cout << "  Output regressor FN " << noOfOutGeneratedRegressorObservationsFilenames << ":    " << outGeneratedRegressorObservationsFilenames[noOfOutGeneratedRegressorObservationsFilenames-1] << std::endl;
      continue;
    }
    if( outGeneratedRegressandObservationsFilenames_flag )
    {
      outGeneratedRegressandObservationsFilenames.push_back( argv[i] );
      noOfOutGeneratedRegressandObservationsFilenames++;
      std::cout << "  Output regressand FN " << noOfOutGeneratedRegressandObservationsFilenames << ":    " << outGeneratedRegressandObservationsFilenames[noOfOutGeneratedRegressandObservationsFilenames-1] << std::endl;
      continue;
    }

    // Anything not accounted for?
    std::cout << "Unknown param at " << i << ": " <<  argv[i] << std::endl;
  }


  /* ******************************************************************
   * The following section has to be considered as "proof of principle"
   * of the implementation provided. Numbers are based on an example
   * given in the textbook Fahrmeir et al.:"Multivariate statistische
   * Verfahren."
   * ******************************************************************/

  /*
  imiINFO( "------------------------------------------" );
  imiINFO( " TEST OF FUNCTION: Fahrmeir-Example   " );
  imiINFO( "------------------------------------------" );

  VnlMatrixType regressorMeasurement( 6, 1 );
  MatrixValueType regressorCol[] = { 1.0,  5.0,  14.0,  14.0, 30.0, 27.0};
  regressorMeasurement.set_column( 0, regressorCol );

  VnlMatrixType predictedOutput( 3, 1 );

  imiMLRMotionPrediction* predictionFilter = imiMLRMotionPrediction::New();
    predictionFilter->LoadExampleFahrmeir( regressorMeasurement );
    predictionFilter->TrainLSEstimator();
    predictionFilter->PredictOutput( regressorMeasurement, predictedOutput );

  imiINFO( "--- RESULTS w/o collinearity check -------\n" );

  imiINFO( " Regressor Measurement:\n" );
  vcl_cerr << regressorMeasurement;
  imiINFO( " Predicted Output:\n" );
  vcl_cerr << predictedOutput;

  imiINFO( "--- RESULTS with collinearity check ------\n" );

    predictionFilter->MulticollinearityCheckOn();
    predictionFilter->TrainLSEstimator();
    predictionFilter->PredictOutput( regressorMeasurement, predictedOutput );

  imiINFO( " Regressor Measurement:\n" );
  vcl_cerr << regressorMeasurement;
  imiINFO( " Predicted Output:\n" );
  vcl_cerr << predictedOutput;

  predictionFilter->Delete();
  */

  /* ******************************************************************
   * Now getting to the intended use: Assuming dense vector fields
   * as regressand observations and trying to predict the fields
   * based on sparse motion measurements / indicators.
   * ******************************************************************/

  imiMLRMotionPrediction* predictionFilter = imiMLRMotionPrediction::New();

  // If B is not specified, the prediction filter has to be trained.
  if( inSystemMatrixFilename.empty() )
  {
    imiINFO( "Parameter B_in not set --> Prediction filter has to be trained ...");

    /* ***************************************************************
     *    REGRESSOR PART:
     *****************************************************************
     * Current implementation: Centering inside prediction filter. Thus:
     * First step: Generating regressor matrix for training purposes:
     *
     * If -S is not set or no regressor observation files are specified,
     * it will be checked whether -S_diaphragm is set. Then extract
     * diaphragm motion from fields. Otherwise abort.
     * --> USE CASE 1
     *
     * If -S is set and only a single filename with ending .mat is given,
     * then it is assumed that the file contains the regressand matrix,
     * which will be loaded.
     * --> USE CASE 2
     *
     * If -S is set and a series of filenames is given, they are stacked
     * together into a single regressor matrix.
     * --> USE CASE 3
     *****************************************************************/

    VnlMatrixType regressorTrainingMatrix;
    std::string fileTypePattern;

    // Implementation of USE CASE 01:
    if( noOfRegressorObservations == 0 )
    {
      if( extractDiaphragmMotionFromVectorfields )
      {
        imiINFO( "  Extracting diaphragm motion from regressands." );
        imiINFO( "  NB: Regressands are assumed to be vector fields, given as .mha files." );

        if( !ExtractDiaphragmMotionFromRegressandFields( regressorTrainingMatrix,
                                                         regressandObservationsFilenames,
                                                         outGeneratedRegressorObservationsFilenames,
                                                         leftLungFilename,
                                                         rightLungFilename,
                                                         diaphragmFilename,
                                                         noOfDiaphragmPoints,
                                                         radiusDiaphragmCircle ) )
        {
          imiERROR( "Failed to extract diaphragm motion!" << std::endl) ;
          return -1;
        }
      }
      else
      {
        imiERROR( " Neither regressor observations are specified nor -S_diaphragm set. Aborting computation." );
        return -1;
      }
    }
    else
    {
      fileTypePattern = GetFileEnding( regressorObservationsFilenames[0] );

      // Implementation of USE CASE 02:
      if( (noOfRegressorObservations == 1) && ( strcmp( fileTypePattern.c_str(), "mat" ) == 0 ) )
      {
        imiINFO( "Reading regressor training matrix ..." );
        ReadMatrixFromMatlabFile( regressorTrainingMatrix, regressorObservationsFilenames[0] );
      }
      // Implementation of USE CASE 03:
      else if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
      {
        imiINFO( "Generating regressor training matrix from individual samples ... " );

        VnlMatrixType currentRegressorObservationMatrix;
        for( unsigned int obsCounter = 0; obsCounter < noOfRegressorObservations; obsCounter++ )
        {
          imiINFO( "  Reading observation " << regressorObservationsFilenames[obsCounter] << " ..." );
          ReadMatrixFromMatlabFile( currentRegressorObservationMatrix, regressorObservationsFilenames[obsCounter] );
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
        imiERROR( "Prediction not implemented for regressor file ending: " << fileTypePattern );
        return -1;
      }
    }

    /* ***************************************************************
     * If -S_out is set, write regressor observations into a file.
     *
     * OPTION 1: Only a single FN of ending .mat is given
     * --> Write regressor matrix to file.
     *
     * OPTION 2: If a series of FNs is given:
     * --> Case 1: number of out names = observations
     * --> Write observations to single files.
     * --> Case 2: number of out names = observations + 1
     * --> Write observations to single files + write regressor
     * matrix to last FN
     *****************************************************************/

    if( !outGeneratedRegressorObservationsFilenames.empty() )
    {
      // Converting imiRealType-matrix to double (= matlab standard format):
      vnl_matrix<double> tempRegressorMatrix( regressorTrainingMatrix.rows(),
                                              regressorTrainingMatrix.cols() );
      for( unsigned int k = 0; k < tempRegressorMatrix.rows(); k++ )
      {
        for( unsigned int l = 0; l < tempRegressorMatrix.cols(); l++ )
        {
          tempRegressorMatrix[k][l] = static_cast<double>( regressorTrainingMatrix[k][l] );
        }
      }

      // Now writing files:
      // OPTION 1:
      if( outGeneratedRegressorObservationsFilenames.size() == 1 )
      {
        imiINFO( "Writing regressor matrix" << outGeneratedRegressorObservationsFilenames[0] << " ..." );
        vnl_matlab_filewrite outRegressorFile( outGeneratedRegressorObservationsFilenames[0].c_str() );
        outRegressorFile.write( tempRegressorMatrix, "RegressorMatrix" );
      }
      // OPTION 2:
      else if( (outGeneratedRegressorObservationsFilenames.size() == regressorTrainingMatrix.cols()) ||
               (outGeneratedRegressorObservationsFilenames.size() == regressorTrainingMatrix.cols()+1) )
      {
        // Saving individual matrices:
        vnl_matrix<double> currentRegressorObservationMatrix( regressorTrainingMatrix.rows(), 1 );
        for( unsigned int obsCounter = 0; obsCounter < regressorTrainingMatrix.cols(); obsCounter++ )
        {
          imiINFO( "Writing regressor " << outGeneratedRegressorObservationsFilenames[obsCounter] << " ..." );
          currentRegressorObservationMatrix.set_column( 0, tempRegressorMatrix.get_column( obsCounter ) );
          vnl_matlab_filewrite outRegressorFile( outGeneratedRegressorObservationsFilenames[obsCounter].c_str() );
          outRegressorFile.write( currentRegressorObservationMatrix, "RegressorMatrix" );
        }
        // Saving regressor matrix (if sought ...):
        if( outGeneratedRegressorObservationsFilenames.size() == regressorTrainingMatrix.cols()+1 )
        {
          imiINFO( "Writing regressor matrix" << outGeneratedRegressorObservationsFilenames[regressorTrainingMatrix.cols()] << " ..." );
          vnl_matlab_filewrite outRegressorFile( outGeneratedRegressorObservationsFilenames[regressorTrainingMatrix.cols()].c_str() );
          outRegressorFile.write( tempRegressorMatrix, "RegressorMatrix" );
        }
      }
      else
      {
        imiWARNING( "Number of regressor output FNs wrong! Aborting saving ... " );
      }
    }

    /* ***************************************************************
     *    REGRESSAND PART:
     *****************************************************************
     * STANDARD USE CASE: Regressand observations are given as .mha
     * vector fields.
     * ALTERNATIVE 1: Regressand matrix to be loaded (no of files == 1
     * and ending == .mat).
     * ***************************************************************/

    imiINFO( "Generating regressand training matrix ... " );
    VnlMatrixType regressandTrainingMatrix;

    // Now generating regressand data:
    // STANDARD USE CASE:
    fileTypePattern = GetFileEnding( regressandObservationsFilenames[0] );
    if( fileTypePattern.find("mha",0) != std::string::npos )
    {
      if( regressandObservationsFilenames.size() != regressorTrainingMatrix.cols() )
      {
        imiERROR( "Numbers of regressand and regressor observations are not the same!" );
        return -1;
      }
      if( !GenerateRegressandMatrixFromVectorFields( regressandTrainingMatrix, regressandObservationsFilenames, maskFilename ) )
      {
        imiERROR( "Regressand Matrix generation failed!" );
        return -1;
      }
    }
    // ALTERNATIVE 1:
    else if( fileTypePattern.find("mat",0) != std::string::npos )
    {
      if( noOfRegressandObservations == 1 )
      {
        imiINFO( "Reading regressand training matrix ..." );
        ReadMatrixFromMatlabFile( regressandTrainingMatrix, regressandObservationsFilenames[0] );
      }
      else if( noOfRegressandObservations == regressorTrainingMatrix.cols() )
      {
        imiINFO( "Generating regressand training matrix from individual samples ... " );

        VnlMatrixType currentRegressandObservationMatrix;
        for( unsigned int obsCounter = 0; obsCounter < noOfRegressandObservations; obsCounter++ )
        {
          imiINFO( "  Reading observation " << regressandObservationsFilenames[obsCounter] << " ..." );
          ReadMatrixFromMatlabFile( currentRegressandObservationMatrix, regressandObservationsFilenames[obsCounter] );
          if( obsCounter == 0 )
          {
            regressandTrainingMatrix.set_size( currentRegressandObservationMatrix.cols(), regressandObservationsFilenames.size() );
          }
          regressandTrainingMatrix.set_column( obsCounter, currentRegressandObservationMatrix.get_column( 0 ) );
        }
      }
      else
      {
        imiERROR( "Numbers of regressand and regressor observations are not the same!" );
        return -1;
      }

    }
    else
    {
      imiERROR( "Prediction not implemented for regressand file ending: " << fileTypePattern );
      return -1;
    }

    /* ***************************************************************
     * If -V_out is set, write regressand observations into a file.
     *
     * OPTION 1: Only a single FN of ending .mat is given
     * --> Write regressand matrix to file.
     *
     * OPTION 2: If a series of FNs is given:
     * --> Case 1: number of out names = observations
     * --> Write observations to single files.
     * --> Case 2: number of out names = observations + 1
     * --> Write observations to single files + write regressor
     * matrix to last FN
     *****************************************************************/

    if( !outGeneratedRegressandObservationsFilenames.empty() )
    {
      // Converting imiRealType-matrix to double (= matlab standard format):
      vnl_matrix<double> tempRegressandMatrix( regressandTrainingMatrix.rows(),
                                               regressandTrainingMatrix.cols() );
      for( unsigned int k = 0; k < tempRegressandMatrix.rows(); k++ )
      {
        for( unsigned int l = 0; l < tempRegressandMatrix.cols(); l++ )
        {
          tempRegressandMatrix[k][l] = static_cast<double>( regressandTrainingMatrix[k][l] );
        }
      }

      // Now writing files:
      // OPTION 1:
      if( outGeneratedRegressandObservationsFilenames.size() == 1 )
      {
        imiINFO( "Writing regressand matrix" << regressandObservationsFilenames[0] << " ..." );
        vnl_matlab_filewrite outRegressandFile( outGeneratedRegressandObservationsFilenames[0].c_str() );
        outRegressandFile.write( tempRegressandMatrix, "RegressandMatrix" );
      }
      // OPTION 2:
      else if( (outGeneratedRegressandObservationsFilenames.size() == regressandTrainingMatrix.cols()) ||
               (outGeneratedRegressandObservationsFilenames.size() == regressandTrainingMatrix.cols()+1) )
      {
        // Saving individual matrices:
        vnl_matrix<double> currentRegressandObservationMatrix( regressandTrainingMatrix.rows(), 1 );
        for( unsigned int obsCounter = 0; obsCounter < noOfRegressandObservations; obsCounter++ )
        {
          imiINFO( "Writing regressand " << regressandObservationsFilenames[obsCounter] << " ..." );
          currentRegressandObservationMatrix.set_column( 0, tempRegressandMatrix.get_column( obsCounter ) );
          vnl_matlab_filewrite outRegressandFile( outGeneratedRegressandObservationsFilenames[obsCounter].c_str() );
          outRegressandFile.write( currentRegressandObservationMatrix, "RegressandMatrix" );
        }
        // Saving regressand matrix (if sought ...):
        if( outGeneratedRegressandObservationsFilenames.size() == regressandTrainingMatrix.cols()+1 )
        {
          imiINFO( "Writing regressand matrix" << regressandObservationsFilenames[regressandTrainingMatrix.cols()] << " ..." );
          vnl_matlab_filewrite outRegressandFile( outGeneratedRegressandObservationsFilenames[regressandTrainingMatrix.cols()].c_str() );
          outRegressandFile.write( tempRegressandMatrix, "RegressandMatrix" );
        }
      }
      else
      {
        imiWARNING( "Number of regressand output FNs wrong! Aborting saving ... " );
      }
    }

    /* ***************************************************************
     *    MLR TRAINING PART:
     ****************************************************************/

    std::cout << std::endl;
    imiINFO( "Training OLS system matrix ... " );

    predictionFilter->SetRegressorTrainingMatrix( regressorTrainingMatrix );
    predictionFilter->SetRegressandTrainingMatrix( regressandTrainingMatrix );
    predictionFilter->SetMulticollinearityCheck( true );
    predictionFilter->TrainLSEstimator();

    /* ***************************************************************
     * If B_out is set, then write the system matrix to file.
     * Assumption: Matlab format.
     ****************************************************************/

    if( !outSystemMatrixFilename.empty() )
    {
      // Converting imiRealType-matrix to double (= matlab standard format):
      imiINFO( "  Writing system matrix to file ... " );

      VnlMatrixType lsEstimator;
      predictionFilter->GetLSEstimator( lsEstimator );

      vnl_matrix<double> tempSystemMatrix( lsEstimator.rows(), lsEstimator.cols() );
      for( unsigned int k = 0; k < tempSystemMatrix.rows(); k++ )
      {
        for( unsigned int l = 0; l < tempSystemMatrix.cols(); l++ )
        {
          tempSystemMatrix[k][l] = static_cast<double>( lsEstimator[k][l] );
        }
      }

      vnl_matlab_filewrite outSystemMatrixFile( outSystemMatrixFilename.c_str() );
      outSystemMatrixFile.write( tempSystemMatrix, "SystemMatrix" );
    }

    /* ***************************************************************
     * If S_mean is set, then write S_mean to file.
     * Assumption: Matlab format.
     ****************************************************************/

    if( !meanRegressorFilename.empty() )
    {
      // Converting imiRealType-matrix to double (= matlab standard format):
      imiINFO( "  Writing mean regressor to file ... " );

      VnlMatrixType meanRegressor;
      predictionFilter->GetMeanRegressor( meanRegressor );

      vnl_matrix<double> tempMatrix( meanRegressor.rows(), meanRegressor.cols() );
      for( unsigned int k = 0; k < meanRegressor.rows(); k++ )
      {
        for( unsigned int l = 0; l < meanRegressor.cols(); l++ )
        {
          tempMatrix[k][l] = static_cast<double>( meanRegressor[k][l] );
        }
      }

      vnl_matlab_filewrite outMeanRegressorFile( meanRegressorFilename.c_str() );
      outMeanRegressorFile.write( tempMatrix, "MeanRegressor" );
    }

    /* ***************************************************************
     * If S_mean is set, then write S_mean to file.
     * Assumption: Matlab format.
     * TODO Incorporate possibility of direct .mha file extraction.
     ****************************************************************/

    if( !meanRegressandFilename.empty() )
    {
      // Converting imiRealType-matrix to double (= matlab standard format):
      imiINFO( "  Writing mean regressand to file ... " );

      VnlMatrixType meanRegressand;
      predictionFilter->GetMeanRegressand( meanRegressand );

      vnl_matrix<double> tempMatrix( meanRegressand.rows(), meanRegressand.cols() );
      for( unsigned int k = 0; k < meanRegressand.rows(); k++ )
      {
        for( unsigned int l = 0; l < meanRegressand.cols(); l++ )
        {
          tempMatrix[k][l] = static_cast<double>( meanRegressand[k][l] );
        }
      }

      vnl_matlab_filewrite outMeanRegressandFile( meanRegressandFilename.c_str() );
      outMeanRegressandFile.write( tempMatrix, "MeanRegressand" );
    }
  }
  else
  {
    imiINFO( "Parameter B set. Using the specified matrix for prediction.");
  }

  /* ******************************************************************
   * APPLYING PREDICTION FILTER:
   * ******************************************************************
   * FIRST STEP: Providing prediction information:
   * - B
   * - S_mean
   * - V_mean
   * If -B is specified, the matrix and the mean vals will be loaded.
   * TODO V_mean is currently assumed to be matlab format.
   * ******************************************************************/

  VnlMatrixType systemMatrix;
  VnlMatrixType meanRegressor;
  VnlMatrixType meanRegressand;

  if( !inSystemMatrixFilename.empty() )
  {
    imiINFO( "Loading system matrix and other input from file ... " );
    if( (meanRegressandFilename.empty()) || (meanRegressorFilename.empty()) )
    {
      imiERROR( "Either mean regressand or regressor file not given! Aborting!" );
      return -1;
    }

    // Loading system matrix:
    std::cout << "  System matrix:   " << inSystemMatrixFilename << std::endl;
    ReadMatrixFromMatlabFile( systemMatrix, inSystemMatrixFilename );

    // Loading mean regressor:
    std::cout << "  Mean regressor:  " << meanRegressorFilename << std::endl;
    ReadMatrixFromMatlabFile( meanRegressor, meanRegressorFilename );

    // Loading mean regressand:
    std::cout << "  Mean regressand:  " << meanRegressandFilename << std::endl;
    ReadMatrixFromMatlabFile( meanRegressand, meanRegressandFilename );

    // Updating prediction filter:
    predictionFilter->SetLSEstimator( systemMatrix );
    predictionFilter->SetMeanRegressor( meanRegressor );
    predictionFilter->SetMeanRegressand( meanRegressand );
  }

  /* ******************************************************************
   * SECOND STEP: Prediction
   * Assuming regressor observation given as .mat.
   * Prediction will be saved unter specified name.
   * If ending == .mat --> Matlab
   * If ending == .mha --> Vector-Field
   * ******************************************************************/

  if( (noOfRegressorMeasurements > 0) && (noOfRegressorMeasurements == noOfPredictedOutputsFilenames) )
  {
    // File ending?
    std::string fileTypePattern = GetFileEnding( regressorMeasurementsFilenames[0] );
    std::string outputFileTypePattern = GetFileEnding( predictedOutputsFilenames[0] );

    // Iterate over number of measurements ...
    for( unsigned int i = 0; i < noOfRegressorMeasurements; i++ )
    {
      VnlMatrixType regressorMeasurementMatrix;
      VnlMatrixType predictedOutput;

      // Reading measurement from file:
      if( fileTypePattern.find("mat",0) != std::string::npos )
      {
        ReadMatrixFromMatlabFile( regressorMeasurementMatrix, regressorMeasurementsFilenames[i] );
        vcl_cerr << "MATRIX: " << regressorMeasurementMatrix << std::endl;
      }
      else
      {
        imiERROR( "No implementation for regressor measurement FN ending ... " << fileTypePattern );
        return -1;
      }

      // Runnning prediction:
      imiINFO( "Applying OLS system for prediction ... " );
      predictionFilter->PredictOutput( regressorMeasurementMatrix, predictedOutput );
      //predictionFilter->GetMeanRegressand( predictedOutput );

      // Saving output:
      if( outputFileTypePattern.find("mat",0) != std::string::npos )
      {
        imiINFO( "Saving prediction output ... as .mat" );
        vnl_matrix<double> tempMatrix( predictedOutput.rows(), predictedOutput.cols() );
        for( unsigned int k = 0; k < predictedOutput.rows(); k++ )
        {
          for( unsigned int l = 0; l < predictedOutput.cols(); l++ )
          {
            tempMatrix[k][l] = static_cast<double>( predictedOutput[k][l] );
          }
        }
        vnl_matlab_filewrite outputFile( predictedOutputsFilenames[i].c_str() );
        outputFile.write( tempMatrix, "PredictedOutput" );
      }
      else if( outputFileTypePattern.find("mha",0) != std::string::npos )
      {
        imiINFO( "Saving prediction output ... as .mha" );
        if( !regressandObservationsFilenames.empty() )
        {
          SavePredictedVnlVectorAsVectorField( predictedOutput, predictedOutputsFilenames[i], regressandObservationsFilenames[0], maskFilename );
        }
        else if( !meanRegressandFilename.empty() )
        {
          SavePredictedVnlVectorAsVectorField( predictedOutput, predictedOutputsFilenames[i], meanRegressandFilename, maskFilename );
        }
      }
      else
      {
        imiERROR( "No implementation for prediction output FN ending ... " << outputFileTypePattern );
        return -1;
      }
    }

    /* *****************************************************************/

    predictionFilter->Delete();
  }

  imiINFO( "\n------------------------------------------" );
  imiINFO( "imiSurrogateBasedMotionPrediction finished." );
  imiINFO( "==========================================\n" );

  return 0;
} // end of main


bool ExtractDiaphragmMotionFromRegressandFields( vnl_matrix<imiRealType>& diaphragmMotionMatrix,
                                                 std::vector<std::string> fieldFilenames,
                                                 std::vector<std::string> outRegressorFilenames,
                                                 std::string leftLungFilename,
                                                 std::string rightLungFilename,
                                                 std::string diaphragmFilename,
                                                 unsigned int noOfPoints,
                                                 float innerCircleRadius )
{
  //---------------------------------------------------
  // First step: Loading mask image and detecting inner
  // lung voxels with z-gradient component != 0 and
  // no lung voxels neighbored in inferior direction.
  // Additionally, check if in x and y direction
  // some inner lung neighbors exist.
  //---------------------------------------------------

  //---------------------------------------------------
  // Loading image:
  //---------------------------------------------------

  imiINFO( "  Extracting diaphragm from fields ... " );

  ImageType::Pointer diaphragmImage;
  std::vector<ImageType::IndexType> diaphragmIndices;

  for(unsigned int lung = 0; lung < 2; lung++ )
  {
    std::string currentLungFilename;
    if( lung == 0 )
    {
      imiINFO( "  Loading left lung image ... " );
      std::cout << "\n  " << leftLungFilename;
      currentLungFilename = leftLungFilename;
    }
    else
    {
      imiINFO( "  Loading right lung image ... " );
      std::cout << "\n  " << rightLungFilename;
      currentLungFilename = rightLungFilename;
    }

    ImagePointerType maskImage;
    if( !currentLungFilename.empty() )
    {
      typedef itk::ImageFileReader<ImageType> MaskReaderType;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader->SetFileName( currentLungFilename );
      try
      {
        maskReader->Update();
      }
      catch (itk::ExceptionObject & excep)
      {
        imiERROR(" Load failed with exception:" << excep << "\n");
      }
      maskImage = maskReader->GetOutput();
      maskImage->Update();
      maskImage->DisconnectPipeline();
      if( lung == 0 ) diaphragmImage = maskImage;
    }
    else
    {
      imiERROR( "No mask filename given. Cannot extract diaphragm motion." );
      return false;
    }

    //---------------------------------------------------
    // Check, if mask image contains any non-zero entries
    // if yes: proceed. if no: continue loop.
    //---------------------------------------------------

    typedef itk::ImageRegionIterator<ImageType> SegmentationIteratorType;
    SegmentationIteratorType maskIt( maskImage, maskImage->GetLargestPossibleRegion() );
    unsigned int noOfNonZeroMaskPoints = 0;

    for( maskIt.GoToBegin(); !maskIt.IsAtEnd(); ++maskIt )
    {
      if( maskIt.Get() != 0.0 ) noOfNonZeroMaskPoints++;
    }

    if( noOfNonZeroMaskPoints == 0 ) continue;

    //---------------------------------------------------
    // Compute bounding box of lung segmentation:
    // Nasty itkImageMaskSpatialObject can only handle
    // uchar images ...
    //---------------------------------------------------

    typedef itk::CastImageFilter<ImageType,ImageUChar3DType> CastImageFilterType;
    CastImageFilterType::Pointer castFilter = CastImageFilterType::New();
    castFilter->SetInput( maskImage );

    try
    {
      castFilter->Update();
    }
    catch ( itk::ExceptionObject &err)
    {
      imiERROR( "Can not execute casting! Exception error: " << err );
      return false;
    }

    typedef itk::ImageMaskSpatialObject<Dimension> ImageMaskSpatialObject;
    ImageMaskSpatialObject::Pointer maskSO = ImageMaskSpatialObject::New();
    maskSO->SetImage( castFilter->GetOutput() );
    ImageType::RegionType boundingBox = maskSO->GetAxisAlignedBoundingBoxRegion();

    std::cout << "\n  Bounding box region:\n  " << boundingBox;

    //---------------------------------------------------
    // First guess approximation: Hemidiaphragm location
    // in plane is approx. the center of the bounding box
    // xy-plane.
    //---------------------------------------------------

    imiINFO( "  Determining points on diaphragm ... ");

    ImageType::IndexType hemidiaphragm;
    hemidiaphragm[0] = boundingBox.GetIndex()[0]+boundingBox.GetSize()[0]/2.0;
    hemidiaphragm[1] = boundingBox.GetIndex()[1]+boundingBox.GetSize()[1]*1.8/3.0;
    hemidiaphragm[2] = boundingBox.GetIndex()[2];

    // OK, center point is extracted.
    // Now generate circles around this point
    std::vector<ImageType::IndexType> pointIndicesOnBoundingBox;
    pointIndicesOnBoundingBox.push_back( hemidiaphragm );

    unsigned int noOfCircles = static_cast<unsigned int>( log2( (noOfPoints+3)/2 ) );

    if( noOfCircles > 0 )
    {
      std::cout << "\n  Orig pos: " << hemidiaphragm << std::endl;

      noOfCircles--;
      float circleRadius = 0;
      float circleRadiusStep = innerCircleRadius;
      float initialArc = M_PI/4.0;

      for( unsigned int circleCounter = 0; circleCounter < noOfCircles; circleCounter++ )
      {
        std::cout << "  Circle: " << circleCounter+1 << std::endl;

        circleRadius += circleRadiusStep;
        float arcStep = 2*M_PI/pow(2,(circleCounter+2));
        std::cout << "  Arc step: " << arcStep << std::endl;

        for( unsigned int circlePointCounter = 0; circlePointCounter < pow(2,(circleCounter+2)); circlePointCounter++ )
        {
          ImageType::IndexType additionalPoint = hemidiaphragm;
          additionalPoint[0] += std::polar( circleRadius, initialArc+circlePointCounter*arcStep ).real();
          additionalPoint[1] += std::polar( circleRadius, initialArc+circlePointCounter*arcStep ).imag();
          std::cout << "  Adding pos: " << additionalPoint << std::endl;
          pointIndicesOnBoundingBox.push_back( additionalPoint );
        }
      }
    }

    //---------------------------------------------------
    // In z-direction search for the first lung voxel
    // along a ray.
    // Result: a vector of diaphragm voxel indices.
    //---------------------------------------------------

    ImageType::IndexType currentDiaphragmPos;
    for( unsigned int i = 0; i < pointIndicesOnBoundingBox.size(); i++ )
    {
      currentDiaphragmPos = pointIndicesOnBoundingBox[i];

      while( maskImage->GetLargestPossibleRegion().IsInside( currentDiaphragmPos ) )
      {
        currentDiaphragmPos[2]++;
        if( ( maskImage->GetPixel( currentDiaphragmPos ) ) != 0.0 )
        {
          maskImage->SetPixel( currentDiaphragmPos, 500 );
          break;
        }
        maskImage->SetPixel( currentDiaphragmPos, 1000 );
      }
    }

    //---------------------------------------------------
    // Generating diaphragm image.
    // Meanwhile: Generate and fill a vector with the
    // voxel positions ...
    //---------------------------------------------------

    SegmentationIteratorType diaphragmIt( diaphragmImage, diaphragmImage->GetLargestPossibleRegion() );

    unsigned int noOfDiaphragmVoxels = 0;

    maskIt.GoToBegin();
    diaphragmIt.GoToBegin();
    while ( !maskIt.IsAtEnd() )
    {
      if( (maskIt.Get() != 0) && (maskIt.Get() != 500) && (maskIt.Get() != 1000) )
      {
        maskIt.Set( 0 );
        diaphragmIt.Set( 0 );
      }
      else if( maskIt.Get() == 500 )
      {
        maskIt.Set( 1 );
        diaphragmIt.Set( 1 );
        diaphragmIndices.push_back( maskIt.GetIndex() );
        ++noOfDiaphragmVoxels;
      }
      else if( maskIt.Get() == 1000 )
      {
        diaphragmIt.Set( 2 );
      }
      ++maskIt;
      ++diaphragmIt;
    }
  }

  //---------------------------------------------------
  // Saving diaphragm image.
  //---------------------------------------------------

  /*
  typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
  FieldReaderType::Pointer fieldReader = FieldReaderType::New();
  fieldReader->SetFileName( fieldFilenames[0] );
  try
  {
    fieldReader->Update();
  }
  catch (itk::ExceptionObject & excep)
  {
    imiERROR(" Load failed with exception:" << excep << "\n");
    return false;
  }
  DisplacementFieldPointerType inputField = fieldReader->GetOutput();
  inputField->Update();

  for( unsigned int i=0; i < diaphragmIndices.size(); i++ )
  {
    ImagePixelType diaphragmMotionMagnitude = sqrt( inputField->GetPixel(diaphragmIndices[i])[0]*inputField->GetPixel(diaphragmIndices[i])[0] +
                           inputField->GetPixel(diaphragmIndices[i])[1]*inputField->GetPixel(diaphragmIndices[i])[1] +
                           inputField->GetPixel(diaphragmIndices[i])[2]*inputField->GetPixel(diaphragmIndices[i])[2] );
    diaphragmImage->SetPixel( diaphragmIndices[i], diaphragmMotionMagnitude );
  }
  */

  if( !diaphragmFilename.empty() )
  {
    imiINFO( "  Saving diaphragm image ..." );
    typedef itk::ImageFileWriter<ImageType> ImageWriterType;
    ImageWriterType::Pointer imageWriter = ImageWriterType::New();

    imageWriter->SetInput( diaphragmImage );
    imageWriter->SetFileName( diaphragmFilename.c_str() );

    try
    {
      imageWriter->Update();
    }
    catch (itk::ExceptionObject & excep)
    {
      imiERROR(" Save failed with exception:" << excep << "\n");
      return false;
    }
  }

  //---------------------------------------------------
  // Generating regressand matrix
  //---------------------------------------------------

  imiINFO( "  Filling regressor matrix ... ")

  unsigned int noOfColumns = diaphragmIndices.size()*Dimension; // 2 lungs, given no of points, image dimension
  diaphragmMotionMatrix.set_size( noOfColumns, fieldFilenames.size() );

  typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
  for( unsigned int i = 0; i < fieldFilenames.size(); i++ )
  {
    imiINFO( "  Field: " << fieldFilenames[i] );

    FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( fieldFilenames[i] );
    try
    {
      fieldReader->Update();
    }
    catch (itk::ExceptionObject & excep)
    {
      imiERROR(" Load failed with exception:" << excep << "\n");
      return false;
    }
    DisplacementFieldPointerType inputField = fieldReader->GetOutput();
    inputField->Update();

    unsigned int columnCounter = 0;
    for( unsigned int k = 0; k < diaphragmIndices.size(); k++ )
    {
      for( unsigned int d = 0; d < Dimension; d++ )
      {
        diaphragmMotionMatrix[columnCounter][i] = (inputField->GetPixel( diaphragmIndices[k] ))[d];
        columnCounter++;
      }
    }
  }
  vcl_cerr << "\n  Diaphragm matrix:\n" << diaphragmMotionMatrix;
  return true;
}


bool GenerateRegressandMatrixFromVectorFields( vnl_matrix<imiRealType>& regressandMatrix, std::vector<std::string> fieldFilenames, std::string maskFilename )
{
  /* ******************************************************************
   * If a mask file is given (and this is highly recommended for
   * memory efficiency reasons ...), then use this file for masking
   * the relevant regions.
   * ******************************************************************/

  bool useMask = false;
  unsigned int noOfMaskVoxels = 0;
  ImagePointerType maskImage;

  typedef itk::ImageRegionConstIterator<ImageType> ConstSegmentationIteratorType;
  ConstSegmentationIteratorType maskIt;

  if( !maskFilename.empty() )
  {
    typedef itk::ImageFileReader<ImageType> MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( maskFilename );
    try
    {
      maskReader->Update();
    }
    catch (itk::ExceptionObject & excep)
    {
      imiERROR(" Load failed with exception:" << excep << "\n");
    }
    useMask = true;
    maskImage = maskReader->GetOutput();
    maskImage->Update();

    maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );
    maskIt.GoToBegin();
    while( !maskIt.IsAtEnd() )
    {
      if( maskIt.Get() != 0 ) noOfMaskVoxels++;
      ++maskIt;
    }

    imiINFO(  "Using specified mask for training ..." );
    std::cout << "\n  Number of lung voxels: " << noOfMaskVoxels << std::endl;
  }

  /* ******************************************************************
   * Now coming up for the core functionality:
   * ******************************************************************/

  typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;

  unsigned int noOfRegressandObservations = fieldFilenames.size();
  unsigned int noOfRegressandVectorEntries = 0;

  DisplacementFieldType::SizeType fieldSize;
  ImageType::SizeType maskSize;

  for( unsigned int i = 0; i < noOfRegressandObservations; i++ )
  {
    FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( fieldFilenames[i] );
    try
    {
      fieldReader->Update();
    }
    catch (itk::ExceptionObject & excep)
    {
      imiERROR(" Load failed with exception:" << excep << "\n");
      return -1;
    }
    DisplacementFieldPointerType inputField = fieldReader->GetOutput();

    // If first field is loaded, allocate vnl regressand matrix:
    if( i == 0 )
    {
      if( useMask )
      {
        noOfRegressandVectorEntries = noOfMaskVoxels*Dimension;
        maskSize = maskImage->GetLargestPossibleRegion().GetSize();
        fieldSize = inputField->GetLargestPossibleRegion().GetSize();
      }
      else
      {
        fieldSize = inputField->GetLargestPossibleRegion().GetSize();
        unsigned int noOfVoxels = 1;
        for( unsigned int d = 0; d < Dimension; d++ )
        {
          noOfVoxels *= fieldSize[d];
        }
        noOfRegressandVectorEntries = noOfVoxels*Dimension;
      }

      // Allocate matrix:
      imiINFO("  Allocate regressand VNL matrix (size: " << noOfRegressandVectorEntries << " x " << noOfRegressandObservations << ") ..." );
      if( !regressandMatrix.set_size( noOfRegressandVectorEntries, noOfRegressandObservations) )
      {
        imiERROR( "  Cannot allocate matrix of requested size !");
        return -1;
      }
      regressandMatrix.fill( 0.0 );
    }
    // Otherwise check size. Has to be the same for all fields.
    else
    {
      if( inputField->GetLargestPossibleRegion().GetSize() != fieldSize )
      {
        imiERROR("  Loaded fields have different size!");
        return -1;
      }
      if( useMask && (inputField->GetLargestPossibleRegion().GetSize() != maskSize) )
      {
        imiERROR("  Loaded field and mask have different size!");
        return -1;
      }
    }

    // Now, write field information into vnl matrix:
    imiINFO("  Converting itk field to vnl matrix ...");

    typedef itk::ImageRegionConstIterator<DisplacementFieldType> ConstFieldIteratorType;
    ConstFieldIteratorType inputIt( inputField, inputField->GetLargestPossibleRegion() );
    inputIt.GoToBegin();
    unsigned int vector_pos = 0;
    DisplacementFieldType::PixelType vectorValue;

    // First version (= default version): No mask used ...
    if( !useMask )
    {
      while ( !inputIt.IsAtEnd() )
      {
        vectorValue = inputIt.Get();
        for (unsigned int d = 0; d < Dimension; d++)
        {
          regressandMatrix[vector_pos][i] = vectorValue[d];
          vector_pos++;
        }
        ++inputIt;
      }
    }
    // Second version: Use mask ...
    else
    {
      maskIt.GoToBegin();
      while ( !inputIt.IsAtEnd() )
      {
        if( maskIt.Get() != 0.0 )
        {
          vectorValue = inputIt.Get();
          for (unsigned int d = 0; d < Dimension; d++)
          {
            regressandMatrix[vector_pos][i] = vectorValue[d];
            vector_pos++;
          }
        }
        ++inputIt;
        ++maskIt;
      }
    }

    // Final check:
    if (vector_pos != regressandMatrix.rows())
    {
      imiERROR("  Matrix not correctly filled with field values!");
      return false;
    }
  } // End for-loop over fields.
  return true;
}


std::string GetFileEnding( std::string filename )
{
  std::string pattern = ".";
  std::string::size_type pos = filename.rfind( pattern, filename.length() );
  std::string fileTypePattern = filename.substr(pos+1,3);

  return fileTypePattern;
}

bool ReadMatrixFromMatlabFile( vnl_matrix<imiRealType>& vnlMatrix, std::string filename )
{
  imiINFO( "  Reading file " << filename << " ...\n" );

  vcl_ifstream fid;
  fid.open( filename.c_str() );
  vnl_matrix<double> tempDoubleMatrix;
  vnl_matlab_read_or_die( fid, tempDoubleMatrix );

  vnlMatrix.set_size( tempDoubleMatrix.rows(), tempDoubleMatrix.cols() );
  for( unsigned int k = 0; k < tempDoubleMatrix.rows(); k++ )
  {
    for( unsigned int l = 0; l < tempDoubleMatrix.cols(); l++ )
    {
      vnlMatrix[k][l] = static_cast<imiRealType>( tempDoubleMatrix[k][l] );
    }
  }
  //vcl_cerr << vnlMatrix;
  return true;
}

bool SavePredictedVnlVectorAsVectorField( vnl_matrix<imiRealType>& predictedOutputAsVnlVector,
                                          std::string outFieldFilename,
                                          std::string refFieldname,
                                          std::string maskFilename )
{
  //---------------------------------------------------
  // Load reference field which serves as a template.
  //---------------------------------------------------

  imiINFO( "  Loading template field ..." );

  typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
  FieldReaderType::Pointer fieldReader = FieldReaderType::New();
  fieldReader->SetFileName( refFieldname.c_str() );
  try
  {
    fieldReader->Update();
  }
  catch (itk::ExceptionObject & excep)
  {
    imiERROR(" Load failed with exception:" << excep << "\n");
    return -1;
  }
  DisplacementFieldPointerType outputField = fieldReader->GetOutput();
  DisplacementFieldType::SizeType fieldSize = outputField->GetLargestPossibleRegion().GetSize();

  //---------------------------------------------------
  // If specified: load mask for conversion:
  //---------------------------------------------------

  bool useMask = false;
  unsigned int noOfMaskVoxels = 0;
  ImagePointerType maskImage;

  typedef itk::ImageRegionConstIterator<ImageType> ConstSegmentationIteratorType;
  ConstSegmentationIteratorType maskIt;

  if( !maskFilename.empty() )
  {
    typedef itk::ImageFileReader<ImageType> MaskReaderType;
    MaskReaderType::Pointer maskReader = MaskReaderType::New();
    maskReader->SetFileName( maskFilename );
    try
    {
     maskReader->Update();
    }
    catch (itk::ExceptionObject & excep)
    {
     imiERROR("  Load failed with exception:" << excep << "\n");
    }
    useMask = true;
    maskImage = maskReader->GetOutput();
    maskImage->Update();

    maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );
    maskIt.GoToBegin();
    while( !maskIt.IsAtEnd() )
    {
     if( maskIt.Get() != 0 ) noOfMaskVoxels++;
     ++maskIt;
    }

    imiINFO( "  Using specified mask for conversion ..." );
    std::cout << "\n  Number of lung voxels: " << noOfMaskVoxels << std::endl;
  }

  //---------------------------------------------------
  // Convert vnl matrix to mha vector.
  //---------------------------------------------------

  imiINFO( "  Converting vnl to mha ..." );

  typedef itk::ImageRegionIterator<DisplacementFieldType> FieldIteratorType;
  FieldIteratorType outputIt( outputField, outputField->GetLargestPossibleRegion());
  outputIt.GoToBegin();

  unsigned int urow = 0;
  DisplacementFieldType::PixelType vectorValue;

  // Default version: no lung mask used to fill the field:
  if( !useMask )
  {
    while (!outputIt.IsAtEnd())
    {
      vectorValue[0] = predictedOutputAsVnlVector[urow][0];
      vectorValue[1] = predictedOutputAsVnlVector[urow + 1][0];
      vectorValue[2] = predictedOutputAsVnlVector[urow + 2][0];
      urow += 3;
      outputIt.Set(vectorValue);
      ++outputIt;
    }
  }
  // Alternative version: lung mask used to fill the field:
  // If mask == 0, fill (0,0,0) into field. Otherwise
  // fill in vnl vector entries.
  else
  {
    maskIt.GoToBegin();
    while (!outputIt.IsAtEnd())
    {
      if( maskIt.Get() != 0 )
      {
        vectorValue[0] = predictedOutputAsVnlVector[urow][0];
        vectorValue[1] = predictedOutputAsVnlVector[urow + 1][0];
        vectorValue[2] = predictedOutputAsVnlVector[urow + 2][0];
        urow += 3;
        outputIt.Set(vectorValue);
      }
      else
      {
        vectorValue[0] = 0;
        vectorValue[1] = 0;
        vectorValue[2] = 0;
        outputIt.Set(vectorValue);
      }
      ++outputIt;
      ++maskIt;
    }
  }
  outputField->Modified();

  //---------------------------------------------------
  // Save output field as specified.
  //---------------------------------------------------

  imiINFO( "  Saving mha ..." );

  typedef itk::Image<DisplacementFieldPixelType, Dimension> DisplacementFieldType;
  typedef itk::ImageFileWriter< DisplacementFieldType > FieldWriterType;

  FieldWriterType::Pointer displacementFieldWriter = FieldWriterType::New();

  displacementFieldWriter->SetInput( outputField );
  displacementFieldWriter->SetFileName( outFieldFilename.c_str() );
  try
  {
    displacementFieldWriter->Update();
  }
  catch (itk::ExceptionObject & excep)
  {
    imiERROR(" Save failed with exception:" << excep << "\n");
    return false;
  }

  return true;
}
