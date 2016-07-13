/** \file icnsVariationalRegistration.cpp
 *
 *  Example use case of itk remote module VariationalRegistration.
 *  Here: Either register two images with standard params or select a 4DCT
 *  image folder, extract two images, and register them.
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2016 IPMI at ICNS at UKE
 *
 ****************************************************************************/

// System includes:
#include <iostream>
#include <string>
#include <fstream>
extern "C"
{
#include "getopt.h"
}

// ITK includes:
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkImageSeriesReader.h>
#include <itkGDCMImageIO.h>
#include <itkGDCMSeriesFileNames.h>

#include <itkConfigure.h>
#include <itkExponentialDisplacementFieldImageFilter.h>
#include <itkHistogramMatchingImageFilter.h>

#include <itkVariationalRegistrationMultiResolutionFilter.h>
#include <itkVariationalRegistrationFilter.h>
#include <itkVariationalDiffeomorphicRegistrationFilter.h>
#include <itkVariationalSymmetricDiffeomorphicRegistrationFilter.h>
#include <itkVariationalRegistrationFunction.h>
#include <itkVariationalRegistrationDemonsFunction.h>
#include <itkVariationalRegistrationSSDFunction.h>
#include <itkVariationalRegistrationNCCFunction.h>
#include <itkVariationalRegistrationRegularizer.h>
#include <itkVariationalRegistrationGaussianRegularizer.h>
#include <itkVariationalRegistrationDiffusionRegularizer.h>
#include <itkVariationalRegistrationElasticRegularizer.h>
#include <itkVariationalRegistrationStopCriterion.h>
#include <itkVariationalRegistrationLogger.h>

// Project includes: NONE SO FAR

// Global typedefs:

typedef itk::Image<short, 3>                        ImageType;
typedef itk::ImageFileReader<ImageType>             ImageReaderType;
typedef itk::ImageFileWriter<ImageType>             ImageWriterType;
typedef itk::ImageSeriesReader<ImageType>           ImageSeriesReader;
typedef itk::GDCMImageIO                            ImageIOType;
typedef itk::GDCMSeriesFileNames                    NamesGeneratorType;

typedef itk::Image<itk::Vector<float, 3>, 3>        DisplacementFieldType;
typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;
typedef itk::ImageFileWriter<DisplacementFieldType> DisplacementFieldWriterType;

typedef itk::VariationalRegistrationFunction<ImageType,ImageType,DisplacementFieldType>::MaskImageType MaskImageType;
typedef itk::ImageFileReader<MaskImageType>         MaskReaderType;


// ---------------------------------------------------------------
// Declaration of additional routines:
// ---------------------------------------------------------------

std::vector<ImageType::Pointer> ReadTemporalImageSequence( char* directoryFN );

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage icnsVariationalRegistration:\n";
  std::cout << "icnsVariationalRegistration -F <fixed image FN> -M <moving image FN> -O <output field FN> [...]\n\n";
  std::cout << std::endl;
  std::cout << "-F <image FN>              [IN] Filename of the FIXED image." << std::endl;
  std::cout << "-M <image FN>              [IN] Filename of the MOVING image." << std::endl;
  std::cout << "-D <4DCT directory name>   [IN] Directory of 4D data to read." << std::endl;
  std::cout << "-S <segmentation mask>     [IN] Filename of registration mask image (optional)." << std::endl;
  std::cout << "-I <initial vector field>  [IN] Filename of initial deformation field (optional)." << std::endl;
  std::cout << "-O <vector field FN>       [OUT] Filename of the computed DISPLACEMENT field." << std::endl;
  std::cout << "-V <vector field FN>       [OUT] Filename of the computed VELOCITY field." << std::endl;
  std::cout << "-W <warped image FN>       [OUT] Filename of the warped moving image." << std::endl;
  std::cout << std::endl;
  std::cout << "-y <fixed image phase>     If a 4D image data is loaded, specify FIXED image phase." << std::endl;
  std::cout << "-z <moving image phase>    If a 4D image data is loaded, specify MOVING image phase." << std::endl;
  std::cout << std::endl;
  std::cout << "Parameters for registration filter:" << std::endl;
  std::cout << "  -i <iterations>          Number of iterations." << std::endl;
  std::cout << "  -l <levels>              Number of multi-resolution levels." << std::endl;
  std::cout << "  -t <tau>                 Registration time step." << std::endl;
  std::cout << "  -s 0|1|2                 Select search space." << std::endl;
  std::cout << "                             0: Standard (default)." << std::endl;
  std::cout << "                             1: Diffeomorphic." << std::endl;
  std::cout << "                             2: Symmetric diffeomorphic." << std::endl;
  std::cout << "  -u 0|1                   Use spacing for regularization." << std::endl;
  std::cout << "                             0: false" << std::endl;
  std::cout << "                             1: true (default)" << std::endl;
  std::cout << "  -e <exp iterations>      Number of iterations for exponentiator in case of" << std::endl;
  std::cout << "                               diffeomorphic registration (search space 1 or 2)." << std::endl;
  std::cout << std::endl;
  std::cout << "Parameters for regularizer:" << std::endl;
  std::cout << "  -r 0|1|2                 Select regularizer." << std::endl;
  std::cout << "                             0: Gaussian smoother." << std::endl;
  std::cout << "                             1: Diffusive regularizer (default)." << std::endl;
  std::cout << "                             2: Elastic regularizer." << std::endl;
  std::cout << "  -a <alpha>               Alpha for the regularization (only diffusive)." << std::endl;
  std::cout << "  -v <variance>            Variance for the regularization (only gaussian)." << std::endl;
  std::cout << "  -m <mu>                  Mu for the regularization (only elastic)." << std::endl;
  std::cout << "  -b <lambda>              Lambda for the regularization (only elasic)." << std::endl;
  std::cout << std::endl;
  std::cout << "Parameters for registration function:" << std::endl;
  std::cout << "  -f 0|1|2                 Select force term." << std::endl;
  std::cout << "                             0: Demon forces (default)." << std::endl;
  std::cout << "                             1: Sum of Squared Differences." << std::endl;
  std::cout << "                             2: Normalized Cross Correlation." << std::endl;
  std::cout << "  -q <radius>              Radius of neighborhood size for Normalized Cross Correlation." << std::endl;
  std::cout << "  -d 0|1|2                 Select image domain for force calculation." << std::endl;
  std::cout << "                             0: Warped image forces (default)." << std::endl;
  std::cout << "                             1: Fixed image forces." << std::endl;
  std::cout << "                             2: Symmetric forces." << std::endl;
  std::cout << std::endl;
  std::cout << "Parameters for stop criterion:" << std::endl;
  std::cout << "  -p 0|1|2                 Select stop criterion policy for multi-resolution." << std::endl;
  std::cout << "                             0: Use default stop criterion." << std::endl;
  std::cout << "                             1: Use simple graduated policy (default)." << std::endl;
  std::cout << "                             2: Use graduated policy." << std::endl;
  std::cout << "  -g <grad slope>          Set fitted line slope for stop criterion (default 0.005)." << std::endl;
  std::cout << std::endl;
  std::cout << "Preprocessing and general parameters:" << std::endl;
  std::cout << "  -h 0|1                   Perform histogram matching." << std::endl;
  std::cout << "                             0: false (default)" << std::endl;
  std::cout << "                             1: true" << std::endl;
  std::cout << "  -x                       Print debug information during execution." << std::endl;
  std::cout << "-?                         Print this help." << std::endl;
  std::cout << std::endl;
}


// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 2 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  std::cout << " " << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsVariationalRegistration"                << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "Reading parameters and preparing data ..." << std::endl;
  
  // Input / output data
  char* fixedImageFilename               = NULL;
  char* movingImageFilename              = NULL;
  char* warpedImageFilename              = NULL;
  char* temporalImageSequenceDirectory   = NULL;
  char* maskImageFilename                = NULL;
  char* initialDisplacementFieldFilename = NULL;
  char* outputDisplacementFieldFilename  = NULL;
  char* outputVelocityFieldFilename      = NULL;
  
  // If a temporal image sequence is loaded, fixed and moving image frames
  // have to be specified:
  int fixedImagePhase                 = -1;
  int movingImagePhase                = -1;
  
  // Registration parameters
  int numberOfIterations              = 800;   // STANDARD DEFAULT: 400
  int numberOfLevels                  = 4;     // STANDARD DEFAULT: 3
  int numberOfExponentiatorIterations = 4;
  double timestep                     = 1.0;
  int searchSpace                     = 0;     // Standard
  bool useImageSpacing                = true;
  
  // Regularizer parameters
  int regularizerType                 = 1;     // = Diffusive
  float regulAlpha                    = 0.5;
  float regulVar                      = 0.5;
  float regulMu                       = 0.5;
  float regulLambda                   = 0.5;
  
  // Force parameters
  int forceType                       = 0;     // Demon
  int forceDomain                     = 0;     // Warped moving
  int nccRadius                       = 2;
  
  // Stop criterion parameters
  int stopCriterionPolicy             = 1;     // Simple graduated is default
  float stopCriterionSlope            = 1e-05; // STANDARD DEFAULT:  0.005;
  
  // Preproc and general parameters
  bool useHistogramMatching           = false;
  bool useDebugMode                   = false;
  
  int c;
  int intVal = 0;
  
  while( (c = getopt( argc, argv, "F:M:D:S:I:O:V:W:y:z:i:n:l:t:s:u:e:r:a:v:m:b:f:d:p:g:h:q:x?" )) != -1 )
  {
    switch( c )
    {
      case 'F':
        fixedImageFilename = optarg;
        std::cout << std::endl;
        std::cout << "Fixed image filename:               " << fixedImageFilename << std::endl;
        break;
      case 'M':
        movingImageFilename = optarg;
        std::cout << "Moving image filename:              " << movingImageFilename << std::endl;
        break;
      case 'D':
        temporalImageSequenceDirectory = optarg;
        std::cout << "Directory of 4D images:             " << temporalImageSequenceDirectory << std::endl;
        break;
      case 'S':
        maskImageFilename = optarg;
        std::cout << "Mask image filename:                " << maskImageFilename << std::endl;
        break;
      case 'I':
        initialDisplacementFieldFilename = optarg;
        std::cout << "Initial displacement field FN:      " << initialDisplacementFieldFilename << std::endl;
        break;
      case 'O':
        outputDisplacementFieldFilename = optarg;
        std::cout << "Output displacement field FN:       " << outputDisplacementFieldFilename << std::endl;
        break;
      case 'V':
        outputVelocityFieldFilename = optarg;
        std::cout << "Output velocity field FN:           " << outputVelocityFieldFilename << std::endl;
        break;
      case 'W':
        warpedImageFilename = optarg;
        std::cout << "Warped moving image filename:       " << warpedImageFilename << std::endl;
        break;
      case 'y':
        fixedImagePhase = atoi( optarg );
        std::cout << "[MODE: 4D image data] Fixed phase:  " << fixedImagePhase << std::endl;
        break;
      case 'z':
        movingImagePhase = atoi( optarg );
        std::cout << "[MODE: 4D image data] Moving phase: " << movingImagePhase << std::endl;
        break;
      case 'e':
        numberOfExponentiatorIterations = atoi( optarg );
        std::cout << "No. of exp. iterations:             " << numberOfExponentiatorIterations << std::endl;
        break;
      case 'i':
      case 'n':
        numberOfIterations = atoi( optarg );
        std::cout << "No. of iterations:                  " << numberOfIterations << std::endl;
        break;
      case 'l':
        numberOfLevels = atoi( optarg );
        std::cout << "No. of multi-resolution levels:     " << numberOfLevels << std::endl;
        break;
      case 't':
        timestep = atof( optarg );
        std::cout << "Registration time step:             " << timestep << std::endl;
        break;
      case 's':
        searchSpace = atoi( optarg );
        if( searchSpace == 0 )
        {
          std::cout << "Search space:                       Standard" << std::endl;
        }
        else if( searchSpace == 1 )
        {
          std::cout << "Search space:                       Diffeomorphic" << std::endl;
        }
        else if( searchSpace == 2 )
        {
          std::cout << "Search space:                       Symmetric Diffeomorphic" << std::endl;
        }
        else
        {
          std::cerr << "ERROR: Search space unknown!" << std::endl;
          return EXIT_FAILURE;
        }
        break;
      case 'u':
        intVal = atoi( optarg );
        if( intVal == 0 )
        {
          std::cout << "Use image spacing:                  false" << std::endl;
          useImageSpacing = false;
        }
        else
        {
          std::cout << "Use image spacing:                  true" << std::endl;
          useImageSpacing = true;
        }
        break;
      case 'r':
        regularizerType = atoi( optarg );
        if( regularizerType == 0 )
        {
          std::cout << "Regularizer:                        Gaussian" << std::endl;
        }
        else if( regularizerType == 1 )
        {
          std::cout << "Regularizer:                        Diffusive" << std::endl;
        }
        else if( regularizerType == 2 )
        {
          std::cout << "Regularizer:                        Elastic" << std::endl;
        }
        else
        {
          std::cerr << "ERROR: Regularizer space unknown!" << std::endl;
          return EXIT_FAILURE;
        }
        break;
      case 'a':
        regulAlpha = atof( optarg );
        std::cout << "Regularization alpha:               " << regulAlpha << std::endl;
        break;
      case 'v':
        regulVar = atof( optarg );
        std::cout << "Regularization variance:            " << regulVar << std::endl;
        break;
      case 'm':
        regulMu = atof( optarg );
        std::cout << "Regularization mu:                  " << regulMu << std::endl;
        break;
      case 'b':
        regulLambda = atof( optarg );
        std::cout << "Regularization lambda:              " << regulLambda << std::endl;
        break;
      case 'f':
        forceType = atoi( optarg );
        if( forceType == 0 )
        {
          std::cout << "Force type:                         Demons" << std::endl;
        }
        else if ( forceType == 1 )
        {
          std::cout << "Force type:                         SSD" << std::endl;
        }
        else if( forceType == 2)
        {
          std::cout << "Force type:                         NCC" << std::endl;
        }
        else
        {
          std::cerr << "ERROR: Force type unknown!" << std::endl;
          return EXIT_FAILURE;
        }
        break;
      case 'd':
        forceDomain = atoi( optarg );
        if( forceDomain == 0 )
        {
          std::cout << "Force domain:                       Warped moving image" << std::endl;
        }
        else if( forceDomain == 1 )
        {
          std::cout << "Force domain:                       Fixed image" << std::endl;
        }
        else if( forceDomain == 2 )
        {
          std::cout << "Calc. forces in:                    Symmetric" << std::endl;
        }
        else
        {
          std::cerr << "ERROR: Force domain unknown!" << std::endl;
          return EXIT_FAILURE;
        }
        break;
      case 'p':
        stopCriterionPolicy = atoi( optarg );
        if( stopCriterionPolicy == 0 )
        {
          std::cout << "StopCriterion-Policy:               Default stop criterion on all levels." << std::endl;
        }
        else if( stopCriterionPolicy == 1 )
        {
          std::cout << "StopCriterion-Policy:               Simple graduated (- increase count on coarse levels," << std::endl;
          std::cout << "                                                      - plus line fitting on finest level)." << std::endl;
        }
        else if( stopCriterionPolicy == 2 )
        {
          std::cout << "StopCriterion-Policy:               Graduated (- max iterations on coarse levels," << std::endl;
          std::cout << "                                               - increase count on second finest level," << std::endl;
          std::cout << "                                               - plus line fitting on finest level)." << std::endl;
        }
        break;
      case 'g':
        stopCriterionSlope = atof( optarg );
        std::cout << "StopCrit. Grad. Threshold:          " << stopCriterionSlope << std::endl;
        break;
      case 'h':
        intVal = atoi( optarg );
        if( intVal == 0 )
        {
          std::cout << "Use histogram matching:             false" << std::endl;
          useHistogramMatching = false;
        }
        else
        {
          std::cout << "Use histogram matching:             true" << std::endl;
          useHistogramMatching = true;
        }
        break;
      case 'q':
        nccRadius = atoi( optarg );
        std::cout << "Radius size for NCC:                " << nccRadius << std::endl;
        break;
      case 'x':
        std::cout << "Use debug mode:                     true" << std::endl;
        useDebugMode = true;
        break;
      case '?':
        PrintHelp();
        return EXIT_SUCCESS;
      default:
        std::cerr << "ERROR: Unknown argument " << (char) c << " not processed!\n" << std::endl;
        return EXIT_FAILURE;
    }
  }
  
  // Check validity of arguments:
  
  if( fixedImageFilename == NULL && temporalImageSequenceDirectory == NULL )
  {
    std::cerr << "No fixed image!" << std::endl;
    return EXIT_FAILURE;
  }
  if( movingImageFilename == NULL && temporalImageSequenceDirectory == NULL )
  {
    std::cerr << "No moving image!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputDisplacementFieldFilename == NULL )
  {
    std::cerr << "No displacement field filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  if( temporalImageSequenceDirectory != NULL && fixedImagePhase == -1 )
  {
    std::cerr << "4D image mode selected, but no fixed phase specified!" << std::endl;
    return EXIT_FAILURE;
  }
  if( temporalImageSequenceDirectory != NULL && movingImagePhase == -1 )
  {
    std::cerr << "4D image mode selected, but no moving phase specified!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Declaring image pointers and loading input images:

  ImageType::Pointer fixedImage;
  ImageType::Pointer movingImage;
  std::vector<ImageType::Pointer> temporalImageSequence;
  
  MaskImageType::Pointer maskImage;
  DisplacementFieldType::Pointer initialDisplacementField;
  
  DisplacementFieldType::Pointer outputDisplacementField;
  DisplacementFieldType::Pointer outputVelocityField;

  
  std::cout << "--------------------------------" << std::endl;
  std::cout << "Loading input data ..." << std::endl;
  
  // First: loading fixed and moving image
  
  if( temporalImageSequenceDirectory == NULL )
  {
    std::cout << "  -> Reference image ... " << std::flush;
    ImageReaderType::Pointer fixedImageReader = ImageReaderType::New();
    fixedImageReader->SetFileName( fixedImageFilename );
    try
    {
      fixedImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading reference image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    fixedImage = fixedImageReader->GetOutput();
    fixedImage->Update();
    fixedImage->DisconnectPipeline();
    std::cout << "OK." << std::endl;
  
    std::cout << "  Moving image ... " << std::flush;
    ImageReaderType::Pointer movingImageReader = ImageReaderType::New();
    movingImageReader->SetFileName( movingImageFilename );
    try
    {
      movingImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading moving image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    ImageType::Pointer movingImage = movingImageReader->GetOutput();
    movingImage->Update();
    movingImage->DisconnectPipeline();
    std::cout << "OK." << std::endl;
  }
  else
  {
    std::cout << "  Reading temporal image sequence: " << std::endl;
    
    temporalImageSequence = ReadTemporalImageSequence( temporalImageSequenceDirectory );
    
    std::cout << "  - Image sequence length: " << temporalImageSequence.size() << std::endl;
    std::cout << "  - Individual image size: " << temporalImageSequence[0]->GetLargestPossibleRegion().GetSize() << std::endl;
    
    fixedImage = temporalImageSequence[fixedImagePhase];
    fixedImage->Update();
    
    movingImage = temporalImageSequence[movingImagePhase];
    movingImage->Update();
    
    // PSEUDO DEBUG CODE PART: check wether temporal sequence is correctly loaded:
    
    /*
    std::string TEMP_outFixedImage = "./TestOutput/TEMP_" + std::to_string( fixedImagePhase ) + ".nii.gz";

    ImageWriterType::Pointer imageWriter;
    imageWriter = ImageWriterType::New();
    
    imageWriter->SetInput( fixedImage );
    imageWriter->SetFileName( TEMP_outFixedImage );
    imageWriter->Update();
    
    std::string TEMP_outMovingImage = "./TestOutput/TEMP_" + std::to_string( movingImagePhase ) + ".nii.gz";
    
    imageWriter = ImageWriterType::New();
    
    imageWriter->SetInput( movingImage );
    imageWriter->SetFileName( TEMP_outMovingImage );
    imageWriter->Update();
    */
  }
  
  // Loading mask image (if specified):
  
  if( maskImageFilename != NULL )
  {
    std::cout << "Loading mask image ... " << std::flush;
    MaskReaderType::Pointer maskImageReader = MaskReaderType::New();
    maskImageReader->SetFileName( maskImageFilename );
    try
    {
      maskImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading mask image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    maskImage = maskImageReader->GetOutput();
    maskImage->Update();
    maskImage->DisconnectPipeline();
    std::cout << "OK." << std::endl;
  }
  
  // Loading initial displacement field (if specified):
  
  if( initialDisplacementFieldFilename != NULL )
  {
    std::cout << "Loading initial field..."  << std::flush;
    DisplacementFieldReaderType::Pointer displacementFieldReader = DisplacementFieldReaderType::New();
    displacementFieldReader->SetFileName( initialDisplacementFieldFilename );
    try
    {
      displacementFieldReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading initial displacement field." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    initialDisplacementField = displacementFieldReader->GetOutput();
    initialDisplacementField->Update();
    initialDisplacementField->DisconnectPipeline();
    std::cout << "OK" << std::endl;
  }
  
  // -------------------------------------------------------------
  // Preprocessing of image data:
  // Histogram matching
  
  if( useHistogramMatching )
  {
    std::cout << "Performing histogram matching of moving image..." << std::endl;
    typedef itk::HistogramMatchingImageFilter<ImageType, ImageType> MatchingFilterType;
    MatchingFilterType::Pointer matcher;
    
    matcher = MatchingFilterType::New();
    
    matcher->SetInput( movingImage );
    matcher->SetReferenceImage( fixedImage );
    matcher->SetNumberOfHistogramLevels( 1024 );
    matcher->SetNumberOfMatchPoints( 7 );
    matcher->ThresholdAtMeanIntensityOn();
    
    try
    {
      matcher->Update();
    }
    catch( itk::ExceptionObject& err )
    {
      std::cerr << "ERROR: Could not match input images!" << std::endl;
      return EXIT_FAILURE;
    }
    
    movingImage = matcher->GetOutput();
  }
  
  // -------------------------------------------------------------
  // Registration itself (copied from example file of
  // Variational Registration):
  
  // Setup registration function
  
  typedef itk::VariationalRegistrationFunction<
  ImageType,ImageType,DisplacementFieldType>   FunctionType;
  typedef itk::VariationalRegistrationDemonsFunction<
  ImageType, ImageType, DisplacementFieldType> DemonsFunctionType;
  typedef itk::VariationalRegistrationSSDFunction<
  ImageType, ImageType, DisplacementFieldType> SSDFunctionType;
  typedef itk::VariationalRegistrationNCCFunction<
  ImageType, ImageType, DisplacementFieldType> NCCFunctionType;
  
  FunctionType::Pointer function;
  switch( forceType )
  {
    case 0:
    {
      DemonsFunctionType::Pointer demonsFunction = DemonsFunctionType::New();
      switch( forceDomain )
      {
        case 0:
          demonsFunction->SetGradientTypeToWarpedMovingImage();
          break;
        case 1:
          demonsFunction->SetGradientTypeToFixedImage();
          break;
        case 2:
          demonsFunction->SetGradientTypeToSymmetric();
          break;
      }
      function = demonsFunction;
    }
    break;
    case 1:
    {
      SSDFunctionType::Pointer ssdFunction = SSDFunctionType::New();
      switch( forceDomain )
      {
        case 0:
          ssdFunction->SetGradientTypeToWarpedMovingImage();
          break;
        case 1:
          ssdFunction->SetGradientTypeToFixedImage();
          break;
        case 2:
          ssdFunction->SetGradientTypeToSymmetric();
          break;
      }
      function = ssdFunction;
    }
    break;
    case 2:
      NCCFunctionType::Pointer nccFunction = NCCFunctionType::New();
      NCCFunctionType::RadiusType r;
      for( unsigned int dim = 0; dim < NCCFunctionType::ImageDimension; dim++ )
      {
        r[dim] = nccRadius;
      }
      nccFunction->SetRadius( r );
      
      switch( forceDomain )
      {
        case 0:
          nccFunction->SetGradientTypeToWarpedMovingImage();
          break;
        case 1:
          nccFunction->SetGradientTypeToFixedImage();
          break;
        case 2:
          nccFunction->SetGradientTypeToSymmetric();
          break;
      }
      function = nccFunction;
      break;
  }
  //function->SetMovingImageWarper( warper );
  function->SetTimeStep( timestep );
  
  // Setup regularizer
  typedef itk::VariationalRegistrationRegularizer<DisplacementFieldType>          RegularizerType;
  typedef itk::VariationalRegistrationGaussianRegularizer<DisplacementFieldType>  GaussianRegularizerType;
  typedef itk::VariationalRegistrationDiffusionRegularizer<DisplacementFieldType> DiffusionRegularizerType;
  typedef itk::VariationalRegistrationElasticRegularizer<DisplacementFieldType>   ElasticRegularizerType;
  
  RegularizerType::Pointer regularizer;
  switch( regularizerType )
  {
    case 0:
    {
      GaussianRegularizerType::Pointer gaussRegularizer = GaussianRegularizerType::New();
      gaussRegularizer->SetStandardDeviations( vcl_sqrt( regulVar ) );
      regularizer = gaussRegularizer;
    }
    break;
    case 1:
    {
      DiffusionRegularizerType::Pointer diffRegularizer = DiffusionRegularizerType::New();
      diffRegularizer->SetAlpha( regulAlpha );
      regularizer = diffRegularizer;
    }
    break;
    case 2:
    {
      ElasticRegularizerType::Pointer elasticRegularizer = ElasticRegularizerType::New();
      elasticRegularizer->SetMu( regulMu );
      elasticRegularizer->SetLambda( regulLambda );
      regularizer = elasticRegularizer;
    }
    break;
  }
  regularizer->InPlaceOff();
  regularizer->SetUseImageSpacing( useImageSpacing );
  
  // Setup registration filter
  typedef itk::VariationalRegistrationFilter<
  ImageType,ImageType,DisplacementFieldType> RegistrationFilterType;
  typedef itk::VariationalDiffeomorphicRegistrationFilter<
  ImageType,ImageType,DisplacementFieldType> DiffeomorphicRegistrationFilterType;
  typedef itk::VariationalSymmetricDiffeomorphicRegistrationFilter<
  ImageType,ImageType,DisplacementFieldType> SymmetricDiffeomorphicRegistrationFilterType;
  
  RegistrationFilterType::Pointer regFilter;
  switch( searchSpace )
  {
    case 0:
    {
      regFilter = RegistrationFilterType::New();
    }
    break;
    case 1:
    {
      DiffeomorphicRegistrationFilterType::Pointer diffeoRegFilter =
      DiffeomorphicRegistrationFilterType::New();
      diffeoRegFilter->SetNumberOfExponentiatorIterations( numberOfExponentiatorIterations );
      regFilter = diffeoRegFilter;
    }
    break;
    case 2:
    {
      SymmetricDiffeomorphicRegistrationFilterType::Pointer symmDiffeoRegFilter =
      SymmetricDiffeomorphicRegistrationFilterType::New();
      symmDiffeoRegFilter->SetNumberOfExponentiatorIterations( numberOfExponentiatorIterations );
      regFilter = symmDiffeoRegFilter;
    }
    break;
  }
  regFilter->SetRegularizer( regularizer );
  regFilter->SetDifferenceFunction( function );
  
  // Setup multi-resolution filter
  
  unsigned int its[numberOfLevels];
  its[numberOfLevels - 1] = numberOfIterations;
  for( int level = numberOfLevels - 2; level >= 0; --level )
  {
    its[level] = its[level + 1];
  }
  
  typedef itk::VariationalRegistrationMultiResolutionFilter<ImageType,ImageType,DisplacementFieldType> MRRegistrationFilterType;
  
  MRRegistrationFilterType::Pointer mrRegFilter = MRRegistrationFilterType::New();
  mrRegFilter->SetRegistrationFilter( regFilter );
  mrRegFilter->SetMovingImage( movingImage );
  mrRegFilter->SetFixedImage( fixedImage );
  mrRegFilter->SetMaskImage( maskImage );
  mrRegFilter->SetNumberOfLevels( numberOfLevels );
  mrRegFilter->SetNumberOfIterations( its );
  mrRegFilter->SetInitialField( initialDisplacementField );
  
  //
  // Setup stop criterion
  //
  typedef itk::VariationalRegistrationStopCriterion<
  RegistrationFilterType,MRRegistrationFilterType> StopCriterionType;
  StopCriterionType::Pointer stopCriterion = StopCriterionType::New();
  stopCriterion->SetRegressionLineSlopeThreshold( stopCriterionSlope );
  stopCriterion->PerformLineFittingMaxDistanceCheckOn();
  
  switch( stopCriterionPolicy )
  {
    case 1:
      stopCriterion->SetMultiResolutionPolicyToSimpleGraduated();
      break;
    case 2:
      stopCriterion->SetMultiResolutionPolicyToGraduated();
      break;
    default:
      stopCriterion->SetMultiResolutionPolicyToDefault();
      break;
  }
  
  regFilter->AddObserver( itk::IterationEvent(), stopCriterion );
  mrRegFilter->AddObserver( itk::IterationEvent(), stopCriterion );
  mrRegFilter->AddObserver( itk::InitializeEvent(), stopCriterion );
  
  // Setup logger
  
  typedef itk::VariationalRegistrationLogger<
  RegistrationFilterType,MRRegistrationFilterType> LoggerType;
  LoggerType::Pointer logger = LoggerType::New();
  
  regFilter->AddObserver( itk::IterationEvent(), logger );
  mrRegFilter->AddObserver( itk::IterationEvent(), logger );
  
  if( useDebugMode )
  {
    regularizer->DebugOn();
    regFilter->DebugOn();
    mrRegFilter->DebugOn();
    stopCriterion->DebugOn();
    logger->DebugOn();
  }
  
  // Execute registration
  
  std::cout << "Starting registration..." << std::endl;
  
  mrRegFilter->Update();
  
  std::cout << "Registration execution finished." << std::endl;
  
  outputDisplacementField = mrRegFilter->GetDisplacementField();
  if( searchSpace == 1 || searchSpace == 2 )
  {
    outputVelocityField = mrRegFilter->GetOutput();
  }
  
  
  // Write results
  
  std::cout << "==========================================" << std::endl;
  std::cout << "WRITING output data..." << std::endl;
  
  if( outputDisplacementFieldFilename != NULL && outputDisplacementField.IsNotNull() )
  {
    std::cout << "Saving deformation field..." << std::endl;
    DisplacementFieldWriterType::Pointer DisplacementFieldWriter;
    DisplacementFieldWriter = DisplacementFieldWriterType::New();
      
    DisplacementFieldWriter->SetInput( outputDisplacementField );
    DisplacementFieldWriter->SetFileName( outputDisplacementFieldFilename );
    DisplacementFieldWriter->Update();
      
    if( outputVelocityFieldFilename != NULL && outputVelocityField.IsNotNull() )
    {
      std::cout << "Saving velocity field..." << std::endl;
      DisplacementFieldWriterType::Pointer velocityFieldWriter;
      velocityFieldWriter = DisplacementFieldWriterType::New();
        
      velocityFieldWriter->SetInput( outputVelocityField );
      velocityFieldWriter->SetFileName( outputVelocityFieldFilename );
      velocityFieldWriter->Update();
    }
  }
  
  if( warpedImageFilename != NULL )
  {
    typedef FunctionType::MovingImageWarperType MovingImageWarperType;
    MovingImageWarperType::Pointer warper = MovingImageWarperType::New();
    
    warper->SetInput( movingImage );
    warper->SetOutputParametersFromImage( fixedImage );
    warper->SetDisplacementField( outputDisplacementField );
    warper->UpdateLargestPossibleRegion();
    
    ImageWriterType::Pointer imageWriter;
    imageWriter = ImageWriterType::New();
    
    imageWriter->SetInput( warper->GetOutput() );
    imageWriter->SetFileName( warpedImageFilename );
    imageWriter->Update();
  }
  
  std::cout << "icnsVariationalRegistration FINISHED!" << std::endl;
  std::cout << "==========================================\n\n" << std::endl;
  
  return EXIT_SUCCESS;
}

// ---------------------------------------------------------------
// Implementation of aforementioned routine declarations:
// ---------------------------------------------------------------

std::vector<ImageType::Pointer> ReadTemporalImageSequence( char* directoryFN )
{
  std::vector<ImageType::Pointer> temporalImageSequence;
  
  // Generate reader:
  
  ImageIOType::Pointer imageIO = ImageIOType::New();
  ImageSeriesReader::Pointer seriesReader = ImageSeriesReader::New();
  seriesReader->SetImageIO( imageIO );
  
  // Generate filenames to read:
  
  NamesGeneratorType::Pointer namesGenerator = NamesGeneratorType::New();
  namesGenerator->SetDirectory( directoryFN );
  
  typedef std::vector<std::string> SeriesIDContainer;
  const SeriesIDContainer& seriesUID = namesGenerator->GetSeriesUIDs();
  SeriesIDContainer::const_iterator seriesItr = seriesUID.begin();
  SeriesIDContainer::const_iterator seriesEnd = seriesUID.end();
  
  std::vector<std::string> filenames;
  while( seriesItr != seriesEnd )
  {
    // Define image to be read:
    
    ImageType::Pointer currentImage;
    
    // Reading data:
    
    std::cout << "  Reading series " << seriesItr->c_str() << " ... " << std::flush;
    filenames = namesGenerator->GetFileNames( seriesItr->c_str() );
    seriesReader->SetFileNames( filenames );
    
    try
    {
      seriesReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while reading series." << std::endl;
      std::cerr << excp << std::endl;
    }
    
    // Writing data into vector:
    
    currentImage = seriesReader->GetOutput();
    currentImage->Update();
    currentImage->DisconnectPipeline();
    temporalImageSequence.push_back( currentImage );
    
    std::cout << "OK" << std::endl;
    seriesItr++;
  }
  
  return temporalImageSequence;
}
