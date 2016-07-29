/** \file imiGenericCCATrainingAndPredictionMain.cxx
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
#include "vnl/algo/vnl_svd.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// Project includes:
#include "imiCCAMotionPrediction.h"
#include "imiMotionPredictionLOOCV.h"
#include "imiMacro.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace imi;


void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenericCCATrainingAndPrediction ... \n\n";

  std::cout << "-V            List of regressand observations (.mat|.mha) | \n";
  std::cout << "-Z            List of regressor observations (.mat|.mha)\n";
  std::cout << "-I            Input measurement (.mat).\n";
  std::cout << "-O            Prediction output (.mat)|(.mha).\n";
  std::cout << "-M            Mask file to restrict the computations, but only usable with -O *.mha.";
  std::cout << "-m            1(default)|2 1: Only the new regressor basis is used for prediction. 2: Use new regressors and regressand bases for prediction.";
  std::cout << "-c            Number of componentes used for prediction.\n";
  std::cout << "-a            Variability threshold for regressand PCA.\n";
  std::cout << "-o            Variability threshold for regressor PCA.\n";
  std::cout << "-l            Perform LOOCV to determine optimal number of components used for prediction.\n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{
  imiINFO( "===================================================" );
  imiINFO( "====    imiGenericCCATrainingAndPrediction     ====" );
  imiINFO( "===================================================" );
  imiINFO( "Reading parameters ...\n" );

  if( argc < 2 )
  {
    PrintHelp();
    return 1;
  }

  imiObject::SetGlobalDebugLevel( 8 );

  // Initialize time for log-output
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

  bool useBothBases=false;

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

  imiRealType regressorPCAvariabilityThreshold=1.00;
  imiRealType regressandPCAvariabilityThreshold=1.00;

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
      inputObservations_flag = false;
      outputFilenames_flag = false;
      majorFlagChange = false;
    }

    // Regressor related params:
    if( strcmp( argv[i], "-Z" ) == 0 )
    {
      regressorObservations_flag = true;
      regressandObservations_flag = false;
      inputObservations_flag = false;
      outputFilenames_flag = false;
      std::cout << "Reading filenames of regressor observations:" << std::endl;
      continue;
    }
    // Regressand related params:
    if( strcmp( argv[i], "-V" ) == 0 )
    {
      regressorObservations_flag = false;
      regressandObservations_flag = true;
      inputObservations_flag = false;
      outputFilenames_flag = false;
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

    if( strcmp( argv[i], "-c" ) == 0 )
    {
      i++;
      noOfComponents = atoi(argv[i]);
      std::cout << "Number of components used for prediction:        " << noOfComponents << std::endl;
      majorFlagChange = true;
      continue;
    }

    if( strcmp( argv[i], "-a" ) == 0 )
    {
      i++;
      regressandPCAvariabilityThreshold = atof(argv[i]);
      std::cout << "Variability threshold for regressand PCA:        " << regressandPCAvariabilityThreshold << std::endl;
      majorFlagChange = true;
      continue;
    }

    if( strcmp( argv[i], "-o" ) == 0 )
    {
      i++;
      regressorPCAvariabilityThreshold = atof(argv[i]);
      std::cout << "Variability threshold for regressor PCA:        " << regressorPCAvariabilityThreshold << std::endl;
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
  if( noOfInputObservations != noOfOutputFilenames )
  {
    imiERROR( "Numbers of input observations and output filenames are not the same!" );
    return -1;
  }

  imiCCAMotionPrediction* predictionFilter = imiCCAMotionPrediction::New();
  //predictionFilter->SetNumberOfComponents(noOfComponents);
  predictionFilter->SetUseBothBases(useBothBases);

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
      imiERROR( "CCA regression not implemented for regressor file ending: " << fileTypePattern );
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
      imiERROR( "CCA regression not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }
  }

  /* ***************************************************************
   *    CCA TRAINING PART:
   ****************************************************************/

  VnlMatrixType regressorPCAComponents;
  VnlMatrixType regressorPCAmeanVector;


  VnlMatrixType regressandPCAComponents;
  VnlMatrixType regressandPCAmeanVector;

  //START OF REGRESSOR PCA

  regressorPCAmeanVector.set_size( regressorTrainingMatrix.rows(), 1 );

  for( unsigned int i = 0; i < regressorTrainingMatrix.rows(); i++ )
  {
    regressorPCAmeanVector( i, 0 ) = (regressorTrainingMatrix.get_row( i )).mean();
  }

  // Compute mean regressand matrix entries:

  imiINFO( "  Computing regressand matrix with centered entries ... " );
  VnlMatrixType regressorPCAcenteredMatrix = regressorTrainingMatrix;
  for( unsigned int i = 0; i < regressorTrainingMatrix.rows(); i++ )
  {
    regressorPCAcenteredMatrix.set_row( i, (regressorTrainingMatrix.get_row( i ) - regressorPCAmeanVector( i, 0 )) );
  }

  imiINFO( "  Computing eigenvectors of the regressor covariance matrix ZZ^T (SVD of Z)..." << std::endl );

  vnl_svd<MatrixValueType> svdSolver( regressorPCAcenteredMatrix );

  VnlMatrixType VMatrix = svdSolver.V();
  VnlMatrixType UMatrix = svdSolver.U();
  VnlDiagMatrixType SigmaMatrix = svdSolver.W();

  //std::cout<<"W"<<svdSolver.W()<<std::endl;

  SigmaMatrix = SigmaMatrix * SigmaMatrix;

  // Calculate sum of Eigenvalues and relative Eigenvalues
  double sumOfEigenvalues = 0;
  VnlVectorType RelativeEigenValues( SigmaMatrix.size() );
  for( unsigned int i = 0; i < SigmaMatrix.size(); i++ )
  {
    sumOfEigenvalues += SigmaMatrix[i];
//    std::cout<<"eigenvals: "<<SigmaMatrix[i]<<std::endl;
  }
  for( unsigned int i = 0; i < SigmaMatrix.size(); i++ )
  {
    RelativeEigenValues[i] = SigmaMatrix[i] / sumOfEigenvalues;
  }

  //determine the number of modes needed to describe x% of the variability
  unsigned int regressorPCAnumOfModesUsed = 0;
  double temp = 0;

  for( unsigned int i = 0; (i < SigmaMatrix.size() - 1) && temp < regressorPCAvariabilityThreshold; i++, regressorPCAnumOfModesUsed++ )
  {
    temp += RelativeEigenValues[i];
  }


  imiINFO( "Get first "<<regressorPCAnumOfModesUsed<<" Eigenvektors (cols of U) ... " );
  regressorPCAComponents = VnlMatrixType( regressorPCAcenteredMatrix.rows(), regressorPCAnumOfModesUsed );

  for( unsigned int j = 0; j < regressorPCAComponents.cols(); j++ )
  {
    for( unsigned int i = 0; i < regressorPCAComponents.rows(); i++ )
    {
      regressorPCAComponents[i][j] = svdSolver.U( i, j );
    }
  }

  //Projection of regressor data matrix on pca "subspace"
  regressorTrainingMatrix = regressorPCAComponents.transpose() * regressorPCAcenteredMatrix;

  /*for (unsigned int row = 0; row<regressorTrainingMatrix.rows();row++)
  {
    for (unsigned int col = 0; col<regressorTrainingMatrix.cols();col++)
    {
      regressorTrainingMatrix[row][col]=regressorTrainingMatrix[row][col]/sqrt(SigmaMatrix[row]);
    }
  }*/

  //END OF REGRESSOR PCA

  //START OF REGRESSAND PCA

  regressandPCAmeanVector.set_size( regressandTrainingMatrix.rows(), 1 );

  for( unsigned int i = 0; i < regressandTrainingMatrix.rows(); i++ )
  {
    regressandPCAmeanVector( i, 0 ) = (regressandTrainingMatrix.get_row( i )).mean();
  }

  // Compute mean regressand matrix entries:

  imiINFO( "  Computing regressand matrix with centered entries ... " );
  VnlMatrixType regressandPCAcenteredMatrix = regressandTrainingMatrix;
  for( unsigned int i = 0; i < regressandTrainingMatrix.rows(); i++ )
  {
    regressandPCAcenteredMatrix.set_row( i, (regressandTrainingMatrix.get_row( i ) - regressandPCAmeanVector( i, 0 )) );
  }

  imiINFO( "  Computing eigenvectors of the implicit regressand covariance matrix V^T*V..." << std::endl );

  vnl_svd<MatrixValueType> svdSolver2( regressandPCAcenteredMatrix.transpose()*regressandPCAcenteredMatrix );

  VnlMatrixType VMatrix2 = svdSolver2.V();
  VnlMatrixType UMatrix2 = svdSolver2.U();
  VnlDiagMatrixType SigmaMatrix2 = svdSolver2.W();



  // Calculate sum of Eigenvalues and relative Eigenvalues
  double sumOfEigenvalues2 = 0;
  VnlVectorType RelativeEigenValues2( SigmaMatrix2.size() );
  for( unsigned int i = 0; i < SigmaMatrix2.size(); i++ )
  {
    sumOfEigenvalues2 += SigmaMatrix2[i];
//    std::cout<<"eigenvals: "<<SigmaMatrix2[i]<<std::endl;
  }
  for( unsigned int i = 0; i < SigmaMatrix2.size(); i++ )
  {
    RelativeEigenValues2[i] = SigmaMatrix2[i] / sumOfEigenvalues2;
  }

  //determine the number of modes needed to describe x% of the variability
  unsigned int regressandPCAnumOfModesUsed = 0;
  temp = 0;

  for( unsigned int i = 0; (i < SigmaMatrix2.size() - 1) && temp < regressandPCAvariabilityThreshold; i++, regressandPCAnumOfModesUsed++ )
  {
    temp += RelativeEigenValues2[i];
  }


  imiINFO("Get first "<<regressandPCAnumOfModesUsed<<" Eigenvektors (cols of U) ... " );
  regressandPCAComponents = VnlMatrixType( regressandPCAcenteredMatrix.rows(), regressandPCAnumOfModesUsed );


  VnlMatrixType newEigenvectors=regressandPCAcenteredMatrix*UMatrix2;


  for( unsigned int j = 0; j < regressandPCAComponents.cols(); j++ )
  {
    for( unsigned int i = 0; i < regressandPCAComponents.rows(); i++ )
    {
      regressandPCAComponents[i][j] = newEigenvectors[ i][j]/sqrt(SigmaMatrix2[j]);
    }
  }

  //Projection of regressand data matrix on pca "subspace"
  regressandTrainingMatrix = regressandPCAComponents.transpose() * regressandPCAcenteredMatrix;

  /*for (unsigned int row = 0; row<regressandTrainingMatrix.rows();row++)
  {
    for (unsigned int col = 0; col<regressandTrainingMatrix.cols();col++)
    {
      regressandTrainingMatrix[row][col]=regressandTrainingMatrix[row][col]/sqrt(SigmaMatrix2[row]);
    }
  }*/


  //END OF REGRESSAND PCA*/*/


  if (loocv)
  {
    imiMotionPredictionLOOCV* loocvTest=imiMotionPredictionLOOCV::New();
    loocvTest->SetRegressandTrainingMatrix(regressandTrainingMatrix);
    loocvTest->SetRegressorTrainingMatrix(regressorTrainingMatrix);
    noOfComponents=loocvTest->PerformLOOCVCCA(1,std::min<unsigned int>(regressorTrainingMatrix.cols(),std::min<unsigned int>(regressorTrainingMatrix.rows(),regressandTrainingMatrix.rows())));
    imiINFO( "LOOCV: "<<noOfComponents<<" selected!");

  } else
  {

    if(noOfComponents > std::min<unsigned int>(regressorTrainingMatrix.cols()-1,std::min<unsigned int>(regressorTrainingMatrix.rows(),regressandTrainingMatrix.rows())))
    {
      noOfComponents=std::min<unsigned int>(regressorTrainingMatrix.cols()-1,std::min<unsigned int>(regressorTrainingMatrix.rows(),regressandTrainingMatrix.rows()));
      imiWARNING("The number of eigenvectors specified exceeds the number of min(#regressor variables,#regressand variables). Setting number of eigenvectors to "<<noOfComponents<<" (min(#regressor variables,#regressand variables))!");
    }
  }

  predictionFilter->SetNumberOfComponents(noOfComponents);

  std::cout << std::endl;
  imiINFO( "Training CCA system matrix ... " );

  predictionFilter->SetRegressorTrainingMatrix( regressorTrainingMatrix );
  predictionFilter->SetRegressandTrainingMatrix( regressandTrainingMatrix );
  predictionFilter->TrainEstimator();


  /* ***************************************************************
   *    CCA PREDICTION PART:
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

  measurement=regressorPCAComponents.transpose()*(measurement-regressorPCAmeanVector);

  /*for (unsigned int row = 0; row<measurement.rows();row++)
  {
      measurement[row][0]=measurement[row][0]/sqrt(SigmaMatrix[row]);
   
  }*/



  predictionFilter->PredictOutput(measurement,prediction);

  /*for (unsigned int row = 0; row<prediction.rows();row++)
  {
      prediction[row][0]=prediction[row][0]*sqrt(SigmaMatrix2[row]);
   
  }*/

  prediction=regressandPCAmeanVector+regressandPCAComponents*prediction;

  imiINFO( "Saving prediction output "<<i+1<<" ... " );

    fileTypePattern = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( predictionFilenames[i] );

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
  imiINFO( "imiGenericCCATraining finished." );
  imiINFO( "==========================================\n" );

  return 0;
} // end of main


