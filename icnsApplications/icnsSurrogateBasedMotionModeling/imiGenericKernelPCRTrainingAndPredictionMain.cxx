/** \file imiGenericKernelPCRTrainingAndPredictionMain.cxx
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
#include "imiKernelPCRMotionPrediction.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"
#include "imiGaussianKernel.h"
#include "imiLinearKernel.h"
#include "imiKernel.h"

using namespace imi;


void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenericKernelPCRTrainingAndPrediction ... \n\n";

  std::cout << "-V            List of regressand observations (.mat|.mha) | \n";
  std::cout << "-Z            List of regressor observations (.mat|.mha)\n";
  std::cout << "-I            Input measurement (.mat).\n";
  std::cout << "-O            Prediction output (.mat)|(.mha).\n";
  std::cout << "-M            Mask file to restrict the computations, but only usable with -O *.mha.";
  std::cout << "-c <modes>    Number of modes determined by the user\n";
  std::cout << "-v <percentage> Number of modes determined by the explained variation criterion (NYI)\n";
  std::cout << "-k 0|1        0: Gaussian kernel 1: Linear Kernel (=Linear PCA)\n";
  std::cout << "-s <sigma>    Sigma for Gaussian kernel (<0 == automatic choice)\n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{
  imiINFO( "===========================================================" );
  imiINFO( "====    imiGenericKernelPCRTrainingAndPrediction     ====" );
  imiINFO( "===========================================================" );
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

  int noOfComponents = -1;

  std::string maskFilename;

  std::string inputFilename;

  std::string predictionFilename;

  bool majorFlagChange = false;

  int kernelChoice=1;

  MatrixValueType sigma=1.0;

  MatrixValueType variability=0.95;

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
      i++;
      inputFilename = argv[i];
      std::cout << "Measurement (Z_hat) filename:        " << inputFilename << std::endl;
      majorFlagChange = true;
      continue;
    }

    //Prediction
    if( strcmp( argv[i], "-O" ) == 0 )
    {
      i++;
      predictionFilename = argv[i];
      std::cout << "Prediction (V_hat) filename:        " << predictionFilename << std::endl;
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

    if( strcmp( argv[i], "-c" ) == 0 )
    {
      i++;
      noOfComponents = atoi(argv[i]);
      std::cout << "Number of components used for prediction:        " << noOfComponents << std::endl;
      majorFlagChange = true;
      continue;
    }

    if( strcmp( argv[i], "-v" ) == 0 )
    {
      i++;
      variability= atof(argv[i]);
      std::cout << "Explained variation criterion:        " << variability << std::endl;
      majorFlagChange = true;
      continue;
    }

    if( strcmp( argv[i], "-s" ) == 0 )
    {
      i++;
      sigma = atof(argv[i]);
      std::cout << "Sigma (Gaussian kernel):        " << sigma << std::endl;
      majorFlagChange = true;
      continue;
    }

    if( strcmp( argv[i], "-k" ) == 0 )
    {
      i++;
      kernelChoice = atoi(argv[i]);
      switch (kernelChoice)
      {
        case 0:
          std::cout << "Kernel:        Gaussian"<< std::endl;
          break;
        case 1:
          std::cout << "Kernel:        Linear"<< std::endl;
          break;
        default:
          std::cout << "ERROR: Unknown kernel!"<< std::endl;
          break;
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

    // Anything not accounted for?
    std::cout << "Unknown param at " << i << ": " << argv[i] << std::endl;
  }

  if( noOfRegressandObservations != noOfRegressorObservations )
  {
    imiERROR( "Numbers of regressand and regressor observations are not the same!" );
    return -1;
  }

  imiKernelPCRMotionPrediction* predictionFilter = imiKernelPCRMotionPrediction::New();


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
      imiERROR( "KernelPCR not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }
  }

  if(noOfComponents > (int)regressorTrainingMatrix.cols()-1)
  {
    imiWARNING("The number of eigenvectors specified exceeds the number of regressor observations-1 available. Setting number of eigenvectors to "<<regressorTrainingMatrix.cols()-1<<" (# of regressor observations-1)!");
    noOfComponents=regressorTrainingMatrix.cols()-1;
  }

  predictionFilter->SetNumberOfComponents(noOfComponents);
  predictionFilter->SetVariabilityThreshold(variability);

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
      imiERROR( "KernelPCR not implemented for regressand file ending: " << fileTypePattern );
      return -1;
    }
  }

  /****************************************************************
   *    KernelPCR TRAINING PART:
   ****************************************************************/

  std::cout << std::endl;
  imiINFO( "Training KernelPCR system matrix ... " );

  imiGaussianKernel* gaussianKernel;
  imiLinearKernel* linearKernel;

  switch (kernelChoice)
    {
      case 0:
        //gaussian kernel
        gaussianKernel=imiGaussianKernel::New();
        predictionFilter->SetKernel(gaussianKernel);
        if(sigma<0)
        {
          sigma=gaussianKernel->ComputeSigmaValueCasero(regressorTrainingMatrix);
          //sigma=gaussianKernel->ComputeSigmaValueCremers1(regressorTrainingMatrix);
        }
        imiINFO( "  sigma (gaussian kernel): "<< sigma );
        //gaussianKernel->SetSigma(sigma);
        //gaussianKernel->SetSigma(100000.0);
        gaussianKernel->SetSigma(sigma);
        break;
      case 1:
        linearKernel=imiLinearKernel::New();
        predictionFilter->SetKernel(linearKernel);
        break;
      default:
        std::cout << "ERROR: Unknown kernel!"<< std::endl;
        break;
    }



  predictionFilter->SetRegressorTrainingMatrix( regressorTrainingMatrix );
  predictionFilter->SetRegressandTrainingMatrix( regressandTrainingMatrix );
  predictionFilter->TrainKernelPCREstimator();


  /****************************************************************
   *    KernelPCR PREDICTION PART:
   ****************************************************************/


  vcl_ifstream fid;

  VnlMatrixType measurement;
  VnlMatrixType prediction;

  imiINFO( "  Reading measurement (Z_hat) ..." );
  fid.open( inputFilename.c_str() );
  vnl_matlab_read_or_die( fid, measurement );
  fid.close();

  predictionFilter->PredictOutput(measurement,prediction);

  if( !predictionFilename.empty() )
  {
    imiINFO( "Saving prediction output ... " );

    fileTypePattern = imiSurrogateBasedMotionPredictionHelpers::GetFileEnding( regressandObservationsFilenames[0] );

    if( strcmp( fileTypePattern.c_str(), "mat" ) == 0 )
    {
      vnl_matlab_filewrite predictionFile( predictionFilename.c_str() );
      predictionFile.write( prediction, "Prediction" );
    } else if ( strcmp( fileTypePattern.c_str(), "mha" ) == 0 )
    {
      imiSurrogateBasedMotionPredictionHelpers::SaveVnlVectorAsVectorField( prediction, predictionFilename, regressandObservationsFilenames[0], maskFilename );
    } else
    {
      imiERROR( "Unsupported file type: " << fileTypePattern );
            return -1;
    }

  }




  predictionFilter->Delete();

  imiINFO( "\n------------------------------------------" );
  imiINFO( "imiGenericKernelPCRTraining finished." );
  imiINFO( "==========================================\n" );

  return 0;
} // end of main
