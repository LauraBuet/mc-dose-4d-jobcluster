/** \file imiGenericKernelMLRTrainingAndPredictionMain.cxx
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
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"
#include "imiMLRMotionPrediction.h"
#include "imiMotionPredictionSFS.h"


using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenericMLRTrainingAndPredictionMethods ... \n\n";

  std::cout << "-V            List of regressand observations (.mat|.mha) | \n";
  std::cout << "-Z            List of regressor observations (.mat|.mha)\n";
  std::cout << "-I            Input measurement (.mat).\n";
  std::cout << "-O            Prediction output (.mat)|(.mha).\n";
  std::cout << "-M            Mask file to restrict the computations, but only usable with -O *.mha.";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{
  imiINFO( "===========================================================" );
  imiINFO( "====    imiGenericMLRTrainingAndPredictionMethods     ====" );
  imiINFO( "===========================================================" );
  imiINFO( "Reading parameters ...\n" );

  if( argc < 2 )
  {
    PrintHelp();
    return 1;
  }

  imiObject::SetGlobalDebugLevel( 8 );

  imiTime::imiInitTime();

  // TYPEDEFS:

  typedef imiRealType MatrixValueType;
  typedef vnl_vector<MatrixValueType> VnlVectorType;
  typedef vnl_matrix<MatrixValueType> VnlMatrixType;

  // VARIABLES AND CALL PARAMS:

  std::vector<std::vector<std::string> > regressorObservationsFilenames;

  std::vector<std::string> regressandObservationsFilenames;

  float ridgeParameter = 0.0;
  float minRidgeParameter = 0.0;
  float maxRidgeParameter = 0.0;
  float multiplicatorRidgeParameter = 0.0;

  std::string maskFilename;

  std::vector<std::vector<std::string> > inputFilenames;

  std::vector<std::string> predictionFilenames;


  bool useSFS = false;
  unsigned int maxNumOfSFSIterations = 1;
  bool combineInputs=false;

  char c;
  int cnt = 0;

  // Reading parameters
  while( (c = getopt( argc, argv, "Z:V:I:O:M:hf:::c?" )) != -1 )
  {
    switch( c )
    {

      case 'Z':
        optind--;
        cnt = 0;
        std::cout << std::endl;
        regressorObservationsFilenames.push_back( std::vector<std::string>() );
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          regressorObservationsFilenames[regressorObservationsFilenames.size() - 1].push_back( argv[optind] );
          std::cout << "  Regressor input [" << regressorObservationsFilenames.size() - 1 << "][" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'V':
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          regressandObservationsFilenames.push_back( argv[optind] );
          std::cout << "  Regressand input [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'I':
        optind--;
        cnt = 0;
        std::cout << std::endl;
        inputFilenames.push_back( std::vector<std::string>() );
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputFilenames[inputFilenames.size() - 1].push_back( argv[optind] );
          std::cout << "  Input [" << inputFilenames.size() - 1 << "][" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        //reads in the output filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          predictionFilenames.push_back( argv[optind] );
          std::cout << "  Output [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'M':
        optind--;
        maskFilename = argv[optind];
        std::cout << "Mask filename " << maskFilename << std::endl;
        break;
      case 'c':
        combineInputs=true;
        std::cout << "Combine inputs:        true" << std::endl;
        break;
      case 'f':
        useSFS = true;
        std::cout << "Use subset selection: true " << std::endl;
        minRidgeParameter = atof( argv[optind] );
        std::cout << "Min ridge parameter:        " << minRidgeParameter << std::endl;
        maxRidgeParameter = atof( argv[++optind] );
        std::cout << "Max ridge parameter:        " << maxRidgeParameter << std::endl;
        multiplicatorRidgeParameter = atof( argv[++optind] );
        std::cout << "Multiplicator ridge parameter:        " << multiplicatorRidgeParameter << std::endl;
        maxNumOfSFSIterations = atoi( argv[++optind] );
        std::cout << "Maximum number of iterations:        " << maxNumOfSFSIterations << std::endl;
        optind++;
        break;
      case 'h':
      case '?':
        PrintHelp();
        exit( 1 );
        break;
      default:
        imiERROR( "Argument "<<(char)c<<" not processed !\n" );
        break;
    }
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

  std::vector<VnlMatrixType*> regressorMatrices;

  for( unsigned int surrogate = 0; surrogate < regressorObservationsFilenames.size(); surrogate++ )
  {

    VnlMatrixType* regressorTrainingMatrix = new VnlMatrixType();
    std::string fileTypePattern;

    fileTypePattern = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( regressorObservationsFilenames[surrogate][0] );

    if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      imiINFO( "Generating regressor training matrix from individual samples ... " );

      VnlMatrixType currentRegressorObservationMatrix;
      for( unsigned int obsCounter = 0; obsCounter < regressorObservationsFilenames[surrogate].size(); obsCounter++ )
      {
        imiINFO( "  Reading observation " << regressorObservationsFilenames[surrogate][obsCounter] << " ..." );
        vcl_ifstream fid;
        fid.open( regressorObservationsFilenames[surrogate][obsCounter].c_str() );
        vnl_matlab_read_or_die( fid, currentRegressorObservationMatrix );
        if( obsCounter == 0 )
        {
          regressorTrainingMatrix->set_size( currentRegressorObservationMatrix.rows(), regressorObservationsFilenames[surrogate].size() );
          regressorTrainingMatrix->fill( 0.0 );
        }
        regressorTrainingMatrix->set_column( obsCounter, currentRegressorObservationMatrix.get_column( 0 ) );
      }

//        std::cout<<"regressorTrainingMatrix"<<*regressorTrainingMatrix<<std::endl;
      regressorMatrices.push_back( regressorTrainingMatrix );
    }
    else
    {
      imiERROR( "MLR not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }

  }

  if (combineInputs)
  {
    VnlMatrixType* combinedMatrix = new VnlMatrixType();
    for (unsigned int i=0; i<regressorMatrices.size();i++)
    {
      combinedMatrix = imiMotionPredictionSFS::CombineMatrices(combinedMatrix,regressorMatrices[i]);
      delete regressorMatrices[i];
    }
    regressorMatrices.clear();
    regressorMatrices.push_back(combinedMatrix);
  }

  /* ***************************************************************
   *    REGRESSAND PART:
   *****************************************************************
   * STANDARD USE CASE: Regressand observations are given as .mat
   * files.
   * ***************************************************************/

  imiINFO( "Generating regressand training matrix ... " );
  VnlMatrixType regressandTrainingMatrix;

  if( regressandObservationsFilenames.size() == 0 )
  {
    imiERROR( " No regressand observations specified. Aborting computation." );
    return -1;
  }
  else
  {
    std::string fileTypePattern;
    fileTypePattern = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( regressandObservationsFilenames[0] );

    if( strcmp( fileTypePattern.c_str(), "mha" ) == 0 )
    {
      imiINFO( "Generating regressand training matrix from individual motion fields ... " );
      imiSurrogateBasedMotionPredictionHelpers::GenerateMatrixFromVectorFields( regressandTrainingMatrix, regressandObservationsFilenames, maskFilename );

    }
    else if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      imiINFO( "Generating regressand training matrix from individual samples ... " );

      VnlMatrixType currentRegressandObservationMatrix;
      for( unsigned int obsCounter = 0; obsCounter < regressandObservationsFilenames.size(); obsCounter++ )
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
      imiERROR( "MLR not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }
  }

  /* ***************************************************************
   *    MLR TRAINING PART:
   ****************************************************************/

  std::cout << std::endl;
  imiINFO( "Training OLS system matrix ... " );

  VnlMatrixType regressorTrainingMatrix( *regressorMatrices[0] );


  /*VnlMatrixType* combined = imiMotionPredictionSFS::CombineMatrices(&regressorTrainingMatrix,regressorMatrices[1]);
  combined = imiMotionPredictionSFS::CombineMatrices(combined,regressorMatrices[2]);
  combined = imiMotionPredictionSFS::CombineMatrices(combined,regressorMatrices[3]);
  combined = imiMotionPredictionSFS::CombineMatrices(combined,regressorMatrices[4]);
  combined = imiMotionPredictionSFS::CombineMatrices(combined,regressorMatrices[5]);
  combined = imiMotionPredictionSFS::CombineMatrices(combined,regressorMatrices[6]);
  combined = imiMotionPredictionSFS::CombineMatrices(combined,regressorMatrices[7]);
  combined = imiMotionPredictionSFS::CombineMatrices(combined,regressorMatrices[8]);

  regressorTrainingMatrix = VnlMatrixType(*combined);*/
  predictionFilter->SetMulticollinearityCheck( true );

  imiMotionPredictionSFS* sfs=NULL;
  if( useSFS )
  {
    sfs = imiMotionPredictionSFS::New();

    sfs->SetRidgeParameters( minRidgeParameter, maxRidgeParameter, multiplicatorRidgeParameter );
    sfs->SetRegressandMatrix( &regressandTrainingMatrix );
    sfs->SetMaxNumOfIterations( maxNumOfSFSIterations );

    for( unsigned int surrogate = 0; surrogate < regressorMatrices.size(); surrogate++ )
    {
      sfs->AddRegressorMatrix( regressorMatrices[surrogate] );
    }

    sfs->PerformSFS();


    //regressorTrainingMatrix = *sfs->GetOptimalSubset();
    ridgeParameter = sfs->GetRidgeParameter();

    std::cout<<"Optimal number of indices: "<<sfs->GetOptimalIndices().size()<<std::endl;
    std::cout << "Optimal indices: ";
    //regressorTrainingMatrix.set_size(sfs->GetOptimalSubset()->rows(),sfs->GetOptimalIndices().size());
    VnlMatrixType tempMatrix;
    std::set<unsigned int> optimalIndices=sfs->GetOptimalIndices();

    for( std::set<unsigned int>::iterator it = optimalIndices.begin(); it != optimalIndices.end(); ++it )
    {
      std::cout << *it << ",";
      VnlMatrixType* combined = imiMotionPredictionSFS::CombineMatrices(&tempMatrix,regressorMatrices[*it]);
      tempMatrix=*combined;
    }
    std::cout<<std::endl;


    regressorTrainingMatrix = tempMatrix;

    std::cout << "Optimal ridge parameter: " << ridgeParameter << std::endl;
    std::cout << "Optimal signal dimension: " << regressorTrainingMatrix.rows() << std::endl;

    predictionFilter->SetMulticollinearityCheck( false );


  }


  //ridgeParameter=1e-05;

  std::cout<<"Regressor:"<<regressorTrainingMatrix<<std::endl;

  predictionFilter->SetRegressorTrainingMatrix( regressorTrainingMatrix );
  predictionFilter->SetRegressandTrainingMatrix( regressandTrainingMatrix );
  predictionFilter->SetTikhonovRegularizationParameter( ridgeParameter );
  predictionFilter->TrainLSEstimator();

  /* ***************************************************************
   *    MLR PREDICTION PART:
   ****************************************************************/

  for( unsigned int i = 0; i < inputFilenames[0].size(); i++ )
  {

    vcl_ifstream fid;

    VnlMatrixType measurement;
    VnlMatrixType prediction;
    std::string fileTypePattern;

    if( useSFS )
    {
      std::set<unsigned int> optimalIndices;

      if (combineInputs)
      {
        for (unsigned int i=0; i<regressorObservationsFilenames.size();i++)
        {
          optimalIndices.insert(i);
        }
      } else
      {
        optimalIndices=sfs->GetOptimalIndices();
      }

      for( std::set<unsigned int>::iterator it = optimalIndices.begin(); it != optimalIndices.end(); ++it )
      {
        VnlMatrixType* tempInput = new VnlMatrixType();
        imiINFO( "  Reading measurement (Z_hat) "<<i+1<<" ..." );
        fid.open( inputFilenames[*it][i].c_str() );
        vnl_matlab_read_or_die( fid, *tempInput );
        fid.close();

        VnlMatrixType* combined = imiMotionPredictionSFS::CombineMatrices(&measurement,tempInput);

        measurement=*combined;

        delete combined;

      }
    }
    else
    {
      imiINFO( "  Reading measurement (Z_hat) "<<i+1<<" ..." );
      fid.open( inputFilenames[0][i].c_str() );
      vnl_matlab_read_or_die( fid, measurement );
      fid.close();
    }

    std::cout << "Measurement: " << measurement << std::endl;

    predictionFilter->PredictOutput( measurement, prediction );

    imiINFO( "Saving prediction output "<<i+1<<" ... " );

    fileTypePattern = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( predictionFilenames[i] );

    if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      vnl_matlab_filewrite predictionFile( predictionFilenames[i].c_str() );
      predictionFile.write( prediction, "Prediction" );
    }
    else if( strcmp( fileTypePattern.c_str(), "mha" ) == 0 )
    {
      imiSurrogateBasedMotionPredictionHelpers::SaveVnlVectorAsVectorField( prediction, predictionFilenames[i], regressandObservationsFilenames[0], maskFilename );
    }
    else
    {
      imiERROR( "Unsupported file type: " << fileTypePattern );
      return -1;
    }

  }

  predictionFilter->Delete();

  imiINFO( "\n--------------------------------------------" );
  imiINFO( "imiGenericMLRTrainingAndPrediction finished." );
  imiINFO( "============================================\n" );

  return 0;
} // end of main


