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
#include "imiKernelMLRMotionPrediction.h"
#include "imiMotionPredictionLOOCV.h"
#include "imiMotionPredictionSFS.h"
#include "imiSurrogateBasedMotionPredictionTypeDefinitions.h"

using namespace imi;

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "imiGenericKernelMLRTrainingAndPrediction ... \n\n";

  std::cout << "-V            List of regressand observations (.mat|.mha) | \n";
  std::cout << "-Z            List of regressor observations (.mat|.mha)\n";
  std::cout << "-I            Input measurement (.mat).\n";
  std::cout << "-O            Prediction output (.mat)|(.mha).\n";
  std::cout << "-M            Mask file to restrict the computations, but only usable with -O *.mha.";
  std::cout << "-k 0|1          0: Gaussian kernel 1: Linear Kernel (=Linear MLR)\n";
  std::cout << "-s <sigma>      Sigma for Gaussian kernel (<0 == automatic choice)\n";
  std::cout << "-h            Print this help.\n";
}

int main( int argc, char *argv[] )
{
  imiINFO( "===========================================================" );
  imiINFO( "====    imiGenericKernelMLRTrainingAndPrediction     ====" );
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

  std::vector<std::vector<std::string> > regressorObservationsFilenames;

  std::vector<std::string> regressandObservationsFilenames;

  float ridgeParameter = 0.0;
  float minRidgeParameter = 0.0;
  float maxRidgeParameter = 0.0;
  float multiplicatorRidgeParameter = 0.0;
  bool loocv = false;

  std::string maskFilename;

  std::vector<std::vector<std::string> > inputFilenames;

  std::vector<std::string> predictionFilenames;

  unsigned int kernelChoice = 0;

  MatrixValueType sigma = 1.0;

  bool useSFS = false;
  unsigned int maxNumOfSFSIterations = 1;

  char c;
  int cnt = 0;

  // Reading parameters
  while( (c = getopt( argc, argv, "Z:V:I:O:M:s:l:::k:hf:::?" )) != -1 )
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
      case 's':
        optind--;
        sigma = atof( argv[optind] );
        std::cout << "Sigma (Gaussian kernel):        " << sigma << std::endl;
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
      case 'l':
        loocv = true;
        std::cout << "Performing LOOCV to determine optimal ridge parameter (only linear kernel!!!!!!!)." << std::endl;
        minRidgeParameter = atof( argv[optind] );
        std::cout << "Min ridge parameter:        " << minRidgeParameter << std::endl;
        maxRidgeParameter = atof( argv[++optind] );
        std::cout << "Max ridge parameter:        " << maxRidgeParameter << std::endl;
        multiplicatorRidgeParameter = atof( argv[++optind] );
        std::cout << "Multiplicator ridge parameter:        " << multiplicatorRidgeParameter << std::endl;
        optind++;
        break;
      case 'k':
        optind--;
        kernelChoice = atoi( argv[optind] );
        switch( kernelChoice )
        {
          case 0:
            std::cout << "Kernel:        Gaussian" << std::endl;
            break;
          case 1:
            std::cout << "Kernel:        Linear" << std::endl;
            break;
          default:
            std::cout << "ERROR: Unknown kernel!" << std::endl;
            break;
        }
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

  imiKernelMLRMotionPrediction* predictionFilter = imiKernelMLRMotionPrediction::New();

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

      regressorMatrices.push_back( regressorTrainingMatrix );
    }
    else
    {
      imiERROR( "KernelMLR not implemented for regressor file ending: " << fileTypePattern );
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
      imiERROR( "KernelMLR not implemented for regressor file ending: " << fileTypePattern );
      return -1;
    }
  }

  /* ***************************************************************
   *    KernelMLR TRAINING PART:
   ****************************************************************/

  std::cout << std::endl;
  imiINFO( "Training OLS system matrix ... " );

  imiGaussianKernel* gaussianKernel;
  imiLinearKernel* linearKernel;

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





  switch( kernelChoice )
  {
    case 0:
      //gaussian kernel
      gaussianKernel = imiGaussianKernel::New();
      predictionFilter->SetKernel( gaussianKernel );
      if( sigma < 0 )
      {
        //sigma = gaussianKernel->ComputeSigmaValueCasero( regressorTrainingMatrix );
        //sigma=gaussianKernel->ComputeSigmaValueCremers1(regressorTrainingMatrix);
        //sigma=gaussianKernel->ComputeSigmaValueCremers2(regressorTrainingMatrix);
        sigma = gaussianKernel->ComputeSigmaValueArias( regressorTrainingMatrix );
      }
      imiINFO( "  sigma (gaussian kernel): "<< sigma );
      //gaussianKernel->SetSigma(sigma);
      //gaussianKernel->SetSigma(100000.0);
      gaussianKernel->SetSigma( sigma );
      break;
    case 1:
      linearKernel = imiLinearKernel::New();
      predictionFilter->SetKernel( linearKernel );
      break;
    default:
      std::cout << "ERROR: Unknown kernel!" << std::endl;
      break;
  }

  if( loocv )
  {
    imiMotionPredictionLOOCV* loocvTest = imiMotionPredictionLOOCV::New();
    loocvTest->SetRegressandTrainingMatrix( regressandTrainingMatrix );
    loocvTest->SetRegressorTrainingMatrix( *regressorMatrices[0] );
    //ridgeParameter=loocvTest->PerformLOOCVLinearKMLR(0.001,1000000000,10);
    switch( kernelChoice )
    {
      case 0:
        //gaussian kernel
        ridgeParameter = loocvTest->PerformLOOCVGaussianKMLR( minRidgeParameter, maxRidgeParameter, multiplicatorRidgeParameter );
        break;
      case 1:
        ridgeParameter = loocvTest->PerformLOOCVLinearKMLR( minRidgeParameter, maxRidgeParameter, multiplicatorRidgeParameter );
        break;
      default:
        std::cout << "ERROR: Unknown kernel!" << std::endl;
        break;
    }

    imiINFO( "LOOCV: "<<ridgeParameter<<" selected!" );
    
      predictionFilter->SetMulticollinearityCheck( false );

  }
  else
  {
    ridgeParameter = 0;
    predictionFilter->SetMulticollinearityCheck( true );
  }

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




  }


  //ridgeParameter=1e-05;

  std::cout<<"Regressor:"<<regressorTrainingMatrix<<std::endl;

  predictionFilter->SetRegressorTrainingMatrix( regressorTrainingMatrix );
  predictionFilter->SetRegressandTrainingMatrix( regressandTrainingMatrix );
  predictionFilter->SetTikhonovRegularizationParameter( ridgeParameter );
  predictionFilter->TrainLSEstimator();

  /* ***************************************************************
   *    KernelMLR PREDICTION PART:
   ****************************************************************/

  for( unsigned int i = 0; i < inputFilenames[0].size(); i++ )
  {

    vcl_ifstream fid;

    VnlMatrixType measurement;
    VnlMatrixType prediction;
    std::string fileTypePattern;

    if( useSFS )
    {
      std::set<unsigned int> optimalIndices=sfs->GetOptimalIndices();

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

  imiINFO( "\n------------------------------------------" );
  imiINFO( "imiGenericKernelMLRTraining finished." );
  imiINFO( "==========================================\n" );

  return 0;
} // end of main

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
    catch( itk::ExceptionObject & excep )
    {
      imiERROR( " Load failed with exception:" << excep << "\n" );
    }
    useMask = true;
    maskImage = maskReader->GetOutput();
    maskImage->Update();

    maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );
    maskIt.GoToBegin();
    while( !maskIt.IsAtEnd() )
    {
      if( maskIt.Get() != 0 )
        noOfMaskVoxels++;
      ++maskIt;
    }

    imiINFO( "Using specified mask for training ..." );
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
    catch( itk::ExceptionObject & excep )
    {
      imiERROR( " Load failed with exception:" << excep << "\n" );
      return -1;
    }
    DisplacementFieldPointerType inputField = fieldReader->GetOutput();

    // If first field is loaded, allocate vnl regressand matrix:
    if( i == 0 )
    {
      if( useMask )
      {
        noOfRegressandVectorEntries = noOfMaskVoxels * Dimension;
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
        noOfRegressandVectorEntries = noOfVoxels * Dimension;
      }

      // Allocate matrix:
      imiINFO( "  Allocate regressand VNL matrix (size: " << noOfRegressandVectorEntries << " x " << noOfRegressandObservations << ") ..." );
      if( !regressandMatrix.set_size( noOfRegressandVectorEntries, noOfRegressandObservations ) )
      {
        imiERROR( "  Cannot allocate matrix of requested size !" );
        return -1;
      }
      regressandMatrix.fill( 0.0 );
    }
    // Otherwise check size. Has to be the same for all fields.
    else
    {
      if( inputField->GetLargestPossibleRegion().GetSize() != fieldSize )
      {
        imiERROR( "  Loaded fields have different size!" );
        return -1;
      }
      if( useMask && (inputField->GetLargestPossibleRegion().GetSize() != maskSize) )
      {
        imiERROR( "  Loaded field and mask have different size!" );
        return -1;
      }
    }

    // Now, write field information into vnl matrix:
    imiINFO( "  Converting itk field to vnl matrix ..." );

    typedef itk::ImageRegionConstIterator<DisplacementFieldType> ConstFieldIteratorType;
    ConstFieldIteratorType inputIt( inputField, inputField->GetLargestPossibleRegion() );
    inputIt.GoToBegin();
    unsigned int vector_pos = 0;
    DisplacementFieldType::PixelType vectorValue;

    // First version (= default version): No mask used ...
    if( !useMask )
    {
      while( !inputIt.IsAtEnd() )
      {
        vectorValue = inputIt.Get();
        for( unsigned int d = 0; d < Dimension; d++ )
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
      while( !inputIt.IsAtEnd() )
      {
        if( maskIt.Get() != 0.0 )
        {
          vectorValue = inputIt.Get();
          for( unsigned int d = 0; d < Dimension; d++ )
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
    if( vector_pos != regressandMatrix.rows() )
    {
      imiERROR( "  Matrix not correctly filled with field values!" );
      return false;
    }
  } // End for-loop over fields.
  return true;
}

bool SavePredictedVnlVectorAsVectorField( vnl_matrix<imiRealType>& predictedOutputAsVnlVector, std::string outFieldFilename, std::string refFieldname, std::string maskFilename )
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
  catch( itk::ExceptionObject & excep )
  {
    imiERROR( " Load failed with exception:" << excep << "\n" );
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
    catch( itk::ExceptionObject & excep )
    {
      imiERROR( "  Load failed with exception:" << excep << "\n" );
    }
    useMask = true;
    maskImage = maskReader->GetOutput();
    maskImage->Update();

    maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );
    maskIt.GoToBegin();
    while( !maskIt.IsAtEnd() )
    {
      if( maskIt.Get() != 0 )
        noOfMaskVoxels++;
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
  FieldIteratorType outputIt( outputField, outputField->GetLargestPossibleRegion() );
  outputIt.GoToBegin();

  unsigned int urow = 0;
  DisplacementFieldType::PixelType vectorValue;

  // Default version: no lung mask used to fill the field:
  if( !useMask )
  {
    while( !outputIt.IsAtEnd() )
    {
      vectorValue[0] = predictedOutputAsVnlVector[urow][0];
      vectorValue[1] = predictedOutputAsVnlVector[urow + 1][0];
      vectorValue[2] = predictedOutputAsVnlVector[urow + 2][0];
      urow += 3;
      outputIt.Set( vectorValue );
      ++outputIt;
    }
  }
  // Alternative version: lung mask used to fill the field:
  // If mask == 0, fill (0,0,0) into field. Otherwise
  // fill in vnl vector entries.
  else
  {
    maskIt.GoToBegin();
    while( !outputIt.IsAtEnd() )
    {
      if( maskIt.Get() != 0 )
      {
        vectorValue[0] = predictedOutputAsVnlVector[urow][0];
        vectorValue[1] = predictedOutputAsVnlVector[urow + 1][0];
        vectorValue[2] = predictedOutputAsVnlVector[urow + 2][0];
        urow += 3;
        outputIt.Set( vectorValue );
      }
      else
      {
        vectorValue[0] = 0;
        vectorValue[1] = 0;
        vectorValue[2] = 0;
        outputIt.Set( vectorValue );
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

  typedef itk::ImageFileWriter<DisplacementFieldType> FieldWriterType;
  FieldWriterType::Pointer displacementFieldWriter = FieldWriterType::New();

  displacementFieldWriter->SetInput( outputField );
  displacementFieldWriter->SetFileName( outFieldFilename.c_str() );
  try
  {
    displacementFieldWriter->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    imiERROR( " Save failed with exception:" << excep << "\n" );
    return false;
  }

  return true;
}

