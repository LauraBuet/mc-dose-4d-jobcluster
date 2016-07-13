/** \file icnsAffineImageRegistration.cpp
 *
 *  \b Initial \b Author: Jan Ehrhardt \n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#include <string>
#include <fstream>
extern "C"
{
#include "getopt.h"
}

#include <itkAffineTransform.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkCastImageFilter.h>
#include <itkCenteredTransformInitializer.h>
#include <itkCommand.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkLinearInterpolateImageFunction.h>
#include <itkMinimumMaximumImageCalculator.h>
#include <itkMeanSquaresImageToImageMetric.h>
#include <itkMultiResolutionImageRegistrationMethod.h>
#include <itkNormalizedCorrelationImageToImageMetric.h>
//#include <itkPowellOptimizer.h>
#include <itkRecursiveMultiResolutionPyramidImageFilter.h>
#include <itkRegularStepGradientDescentOptimizer.h>
#include <itkResampleImageFilter.h>
#include <itkSimilarity3DTransform.h>
#include <itkTransformFileReader.h>
#include <itkTransformFileWriter.h>

// Project includes: NONE so far.

// GLOBAL VARIABLES, part I:

double MaxStepLengthLevel0 = 1.0;
double MinStepLengthLevel0 = 1.0;


// ---------------------------------------------------------------
// DEFINE CLASS CommandIterationUpdate
// ---------------------------------------------------------------

class CommandIterationUpdate: public itk::Command
{
public:
  typedef CommandIterationUpdate  Self;
  typedef itk::Command            Superclass;
  typedef itk::SmartPointer<Self> Pointer;
  //typedef const OptimizerType *OptimizerPointer;
  itkNewMacro( Self );

protected:
  CommandIterationUpdate() :
    m_CumulativeIterationIndex( 0 ),
    m_LevelIterationIndex( 0 ),
    m_CurrentLevel( 0 ),
    m_BestLevel( 0 ),
    m_BestValue( 32000 )
  {};

public:
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef const OptimizerType *                    OptimizerPointer;
  typedef OptimizerType::MeasureType               ValueType;
  typedef OptimizerType::ParametersType            ParametersType;

  void Execute( itk::Object *caller, const itk::EventObject & event )
  {
    Execute( (const itk::Object *) caller, event );
  }

  void Execute( const itk::Object * object, const itk::EventObject & event )
  {
    OptimizerPointer optimizer = dynamic_cast<OptimizerPointer>( object );
    if( !(itk::IterationEvent().CheckEvent( &event )) )
    {
      return;
    }

    if( m_BestValue > optimizer->GetValue() || m_CumulativeIterationIndex==0 )
    {
      m_BestValue = optimizer->GetValue();
      m_BestPosition = optimizer->GetCurrentPosition();
      m_BestLevel = m_CurrentLevel;
      std::cout <<  "New best value:" << std::endl;
    }
    std::cout << std::endl << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << "  " << m_CumulativeIterationIndex << std::endl;

    m_CumulativeIterationIndex++;
    m_LevelIterationIndex++;
  }
  
  unsigned int m_CumulativeIterationIndex;
  unsigned int m_LevelIterationIndex;
  unsigned int m_BestIteration;
  unsigned int m_CurrentLevel, m_BestLevel;
  ValueType m_BestValue;
  ParametersType m_BestPosition;
  
};

// GLOBAL VARIABLES, part II
bool bUseBestValuesPerIteration = true;
bool bUseItsMinMaxStepSchedule = true;
bool bUseBestFromLastLevel = true;
CommandIterationUpdate::Pointer globalObserverPointer = NULL;


// ---------------------------------------------------------------
// DEFINE CLASS RegistrationInterfaceCommand
// ---------------------------------------------------------------

template<typename TRegistration>
class RegistrationInterfaceCommand: public itk::Command
{
public:
  typedef RegistrationInterfaceCommand Self;
  typedef itk::Command                 Superclass;
  typedef itk::SmartPointer<Self>      Pointer;
  itkNewMacro( Self );

protected:
  RegistrationInterfaceCommand()
  {};

public:
  typedef TRegistration                            RegistrationType;
  typedef RegistrationType*                        RegistrationPointer;
  typedef itk::RegularStepGradientDescentOptimizer OptimizerType;
  typedef OptimizerType*                           OptimizerPointer;

  void Execute( itk::Object * object, const itk::EventObject & event )
  {
    if( !(itk::IterationEvent().CheckEvent( &event )) )
    {
      return;
    }
    RegistrationPointer registration = dynamic_cast<RegistrationPointer>( object );
    OptimizerPointer optimizer = dynamic_cast<OptimizerPointer>( registration->GetOptimizer() );

    std::cout << "Start Level " << registration->GetCurrentLevel() << std::endl;

    if( registration->GetCurrentLevel() == 0 )
    {
      optimizer->SetMaximumStepLength( MaxStepLengthLevel0 );
      optimizer->SetMinimumStepLength( MinStepLengthLevel0 );
    }
    else
    {
      if (bUseItsMinMaxStepSchedule)
      {
        optimizer->SetMaximumStepLength( optimizer->GetCurrentStepLength() * 0.5 );
        optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 10.0 );
        optimizer->SetNumberOfIterations( (int) (optimizer->GetNumberOfIterations() * 0.75) + 1 );
      }

      if( bUseBestValuesPerIteration && globalObserverPointer.IsNotNull() )
      {
        // generate a new optimizer as workaround
        std::cout << "               best last level: " << globalObserverPointer->m_BestLevel << std::endl;
        std::cout << " best value in the last levels: " << globalObserverPointer->m_BestValue << std::endl;
        std::cout << "               best parameters: " << globalObserverPointer->m_BestPosition << std::endl;

        registration->SetInitialTransformParametersOfNextLevel( globalObserverPointer->m_BestPosition );
        globalObserverPointer->m_LevelIterationIndex = 0;
        globalObserverPointer->m_CurrentLevel = registration->GetCurrentLevel();
      }
    }
    if ( globalObserverPointer.IsNotNull() && bUseBestFromLastLevel ) globalObserverPointer->m_BestValue *= 2;

    std::cout << " set number of iterations = " << optimizer->GetNumberOfIterations() << std::endl;
    std::cout << "  set maximum step length = " << optimizer->GetMaximumStepLength() << std::endl;
    std::cout << "  set minimum step length = " << optimizer->GetMinimumStepLength() << "\n" << std::endl;
  }
  void Execute( const itk::Object *, const itk::EventObject & )
  {
    return;
  }
};


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage icnsAffineRegistration:\n";
  std::cout << "icnsAffineRegistration -f <fixed image file> -m <moving image file> -o <output image file>\n\n";
    
  std::cout << "-F <image>                 Filename of the fixed image.\n";
  std::cout << "-M <image>                 Filename of the moving image.\n";
  std::cout << "-A <transform>             Filename of the output transform.\n";
  std::cout << "-O <image>                 Filename of the output image.\n";
  std::cout << "-n <iter>                  Number of iterations (affine and similarity, each level).\n";
  std::cout << "-l <levels>                Number of multi-resolution levels.\n";
  std::cout << "-s <length>                Maximum step length of optimization.\n";
  std::cout << "-p 0|1                     Use similarity pre-registration (default on).\n";
  std::cout << "-q <iter>                  Number of similarity iterations (use after -n).\n";
  std::cout << "-c                         Initialize by center of mass matching.\n";
  std::cout << "-u                         Use NCC metric.\n";
  std::cout << "-m 0|1                     Use number of its, Min-/MaxStepLength level schedule (default on).\n";
  std::cout << "-b 0|1|2                   Use best value from all levels and iterations (1), best value from last level (default 2).\n";
    
  std::cout << "-h                         Print this help.\n";
  std::cout << "\n";
}


// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 4 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  // Required typedefs:
  typedef itk::Image<short, 3>                 ImageType;
  typedef itk::ImageFileReader<ImageType>      ImageReaderType;
  typedef itk::ImageFileWriter<ImageType>      ImageWriterType;

  std::cout << "==========================================" << std::endl;
  std::cout << "icnsAffineRegistration"                     << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "Reading parameters ..." << std::endl;

  char* fixedImageFilename       = NULL;
  char* movingImageFilename      = NULL;
  char* outputImageFilename      = NULL;
  char* initialTransformFilename = NULL;
  char* outputTransformFilename  = NULL;
  
  int numberOfIterations = 10;
  int numberOfSimilarityIterations = numberOfIterations;
  int numberOfLevels = 3;
  int intVal = 1;
  bool bUseSimilarityReg = true;
  bool bUseNCCMetric = false;
  bool bUseCenterOfMassForInitialization = false;
  double maxStepLength = 0.1;
  double minStepLength = 0.0001;

  int c;
  while( (c = getopt( argc, argv, "A:F:M:O:n:l:s:p:q:t:m:b:cuh" )) != -1 )
  {
    switch( c )
    {
      case 'F':
        fixedImageFilename = optarg;
        std::cout << std::endl;
        std::cout << " Fixed image filename:            " << fixedImageFilename << std::endl;
        break;
      case 'M':
        movingImageFilename = optarg;
        std::cout << " Moving image filename:           " << movingImageFilename << std::endl;
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << " Output image filename:           " << outputImageFilename << std::endl;
        break;
      case 'A':
        outputTransformFilename = optarg;
        std::cout << " Output transform filename:       " << outputTransformFilename << std::endl;
        break;
      case 'n':
        numberOfIterations = atoi( optarg );
        numberOfIterations = (numberOfIterations < 1) ? 1 : numberOfIterations;
        std::cout << " Number of iterations:            " << numberOfIterations << std::endl;
        numberOfSimilarityIterations = numberOfIterations;
        break;
      case 'l':
        numberOfLevels = atoi( optarg );
        numberOfLevels = (numberOfLevels < 1) ? 1 : numberOfLevels;
        std::cout << " Number of levels:                " << numberOfLevels << std::endl;
        break;
      case 's':
        maxStepLength = atof( optarg );
        std::cout << " Max step length:                 " << maxStepLength << std::endl;
        break;
      case 'c':
        bUseCenterOfMassForInitialization = true;
        std::cout << " Initialize with center of mass." << std::endl;
        break;
      case 'p':
        intVal = atoi( optarg );
        if( intVal == 0 )
        {
          std::cout << " Similarity pre-registration:     false" << std::endl;
          bUseSimilarityReg = false;
        }
        else
        {
          std::cout << " Similarity pre-registration:     true" << std::endl;
          bUseSimilarityReg = true;
        }
        break;
      case 'm':
        intVal = atoi( optarg );
        if( intVal == 0 )
        {
          std::cout << " Scheduling:                      false" << std::endl;
          bUseItsMinMaxStepSchedule = false;
        }
        else
        {
          std::cout << " Scheduling:                      true" << std::endl;
          bUseItsMinMaxStepSchedule = true;
        }
        break;
      case 'b':
        intVal = atoi( optarg );
        if( intVal == 0 )
        {
          std::cout << " Best iteration value:            false" << std::endl;
          bUseBestValuesPerIteration = false;
        }
        else if( intVal == 1 )
        {
          std::cout << " Best all levels iteration value: true" << std::endl;
          bUseBestValuesPerIteration = true;
          bUseBestFromLastLevel = false;
        }
        else if( intVal == 2 )
        {
          std::cout << " Best last level value:           true" << std::endl;
          bUseBestValuesPerIteration = true;
          bUseBestFromLastLevel = true;
        }
        break;
      case 't':
        initialTransformFilename = optarg;
        std::cout << " Initial transform filename:      " << initialTransformFilename << std::endl;
        break;
      case 'u':
        bUseNCCMetric = true;
        std::cout << " Use NCC metric for registration." << std::endl;
        break;
      case 'q':
        numberOfSimilarityIterations = atoi( optarg );
        numberOfSimilarityIterations = (numberOfSimilarityIterations < 1) ? 1 : numberOfSimilarityIterations;
        std::cout << " Number of similarity iterations: " << numberOfSimilarityIterations << std::endl;
        break;
      case 'h':
      case '?':
        PrintHelp();
        exit( 1 );
        break;
      default:
        std::cerr << "ERROR: Argument " << (char) c << " not processed !\n" << std::endl;
    }
  }

  // Check validity of arguments:
  
  if( fixedImageFilename == NULL )
  {
    std::cerr << "No fixed image!" << std::endl;
    return EXIT_FAILURE;
  }
  if( movingImageFilename == NULL )
  {
    std::cerr << "No moving image!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputImageFilename == NULL )
  {
    std::cerr << "No output image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputTransformFilename == NULL )
  {
    std::cerr << "No output transform filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( bUseSimilarityReg && (initialTransformFilename != NULL) )
  {
    std::cerr << "WARNING: Similarity pre-registration and initial affine Transformation not supported!" << std::endl;
    std::cerr << "WARNING: We will NOT use the similarity pre-Registration." << std::endl;
    bUseSimilarityReg = false;
  }
  if( bUseCenterOfMassForInitialization && (initialTransformFilename != NULL) )
  {
    std::cerr << "WARNING: Center of mass initialization and initial affine Transformation not supported!" << std::endl;
    std::cerr << "We will NOT use center of mass initialization." << std::endl;
    bUseCenterOfMassForInitialization = false;
  }

  // -------------------------------------------------------------
  // Loading input images:

  std::cout << "  --------------------------------" << std::endl;
  std::cout << "  Loading input data ..." << std::endl;
  
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
  ImageType::Pointer fixedImage = fixedImageReader->GetOutput();
  fixedImage->Update();
  fixedImage->DisconnectPipeline();
  std::cout << "OK." << std::endl;
  
  std::cout << "  -> Moving image ... " << std::flush;
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

  // Defining registration types:

  typedef float                                                                                 InternalPixelType;
  typedef itk::Image< InternalPixelType, 3 >                                                    InternalImageType;
  typedef itk::AffineTransform< double, 3 >                                                     AffineTransformType;
  typedef itk::Similarity3DTransform< double >                                                  SimilarityTransformType;
  typedef itk::LinearInterpolateImageFunction< InternalImageType, double >                      InterpolatorType;
  typedef itk::NormalizedCorrelationImageToImageMetric< InternalImageType, InternalImageType >  NCCMetricType;
  typedef itk::MeanSquaresImageToImageMetric< InternalImageType, InternalImageType >            MSDMetricType;
  typedef itk::RegularStepGradientDescentOptimizer                                              OptimizerType;
  //typedef OptimizerType::Pointer                                                                OptimizerPointer;
  typedef OptimizerType::ScalesType                                                             OptimizerScalesType;
  typedef itk::MultiResolutionImageRegistrationMethod< InternalImageType, InternalImageType >   RegistrationType;

  // Initializing objects:
  
  InterpolatorType::Pointer interpolator        = InterpolatorType::New();
  RegistrationType::Pointer registration        = RegistrationType::New();
  AffineTransformType::Pointer affineTransform  = AffineTransformType::New();
  SimilarityTransformType::Pointer simTransform = SimilarityTransformType::New();

  // Cast input images to internal image type
  
  typedef itk::CastImageFilter<ImageType, InternalImageType> CastImageFilterType;

  std::cout << "  Casting input images to real valued images ... " << std::flush;
  
  CastImageFilterType::Pointer fixedCaster = CastImageFilterType::New();
  fixedCaster->SetInput( fixedImage );
  fixedCaster->Update();

  CastImageFilterType::Pointer movingCaster = CastImageFilterType::New();
  movingCaster->SetInput( movingImage );
  movingCaster->Update();
  
  std::cout << "OK." << std::endl;

  //  std::cout<<"\nFixed:\n";
  //  fixedCaster->GetOutput()->Print(std::cout);
  //  std::cout<<"\nMoving:\n";
  //  movingCaster->GetOutput()->Print(std::cout);

  // -------------------------------------------------------------
  // Setting up initial parameters

  // Set initial transform parameters (default):
  
  affineTransform->SetIdentity();
  simTransform->SetIdentity();

  if( bUseCenterOfMassForInitialization )
  {
    std::cout << "Initialize with center of mass ..." << std::endl;
    if( bUseSimilarityReg )
    {
      typedef itk::CenteredTransformInitializer<SimilarityTransformType, InternalImageType, InternalImageType> SimTransformInitializerType;
      SimTransformInitializerType::Pointer initializer = SimTransformInitializerType::New();
      initializer->SetTransform( simTransform );
      initializer->SetFixedImage( fixedCaster->GetOutput() );
      initializer->SetMovingImage( movingCaster->GetOutput() );
      initializer->MomentsOn();
      initializer->InitializeTransform();
      std::cout << "  similarity Parameters: [";
      for( unsigned int i = 0; i < simTransform->GetNumberOfParameters(); i++ )
        std::cout << simTransform->GetParameters()[i] << ",";
      std::cout << " ]\n";
    }
    else
    {
      typedef itk::CenteredTransformInitializer<AffineTransformType, InternalImageType, InternalImageType> AffTransformInitializerType;
      AffTransformInitializerType::Pointer initializer = AffTransformInitializerType::New();
      initializer->SetTransform( affineTransform );
      initializer->SetFixedImage( fixedCaster->GetOutput() );
      initializer->SetMovingImage( movingCaster->GetOutput() );
      initializer->MomentsOn();
      initializer->InitializeTransform();
      std::cout << std::endl;
      std::cout << "  affine Parameters: [";
      for( unsigned int i = 0; i < affineTransform->GetNumberOfParameters(); i++ )
        std::cout << affineTransform->GetParameters()[i] << ",";
      std::cout << " ]\n";
    }
  }
  
  // Loading input transforms:
  
  if( initialTransformFilename != NULL)
  {
    std::cout << "  Loading initial transform file ... "  << std::flush;
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
    itk::TransformFileReaderTemplate<double>::Pointer transformReader = itk::TransformFileReaderTemplate<double>::New();
#else
    itk::TransformFileReader::Pointer transformReader = itk::TransformFileReader::New();
#endif
    transformReader->SetFileName( initialTransformFilename );
    try
    {
      transformReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading transform file." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "OK." << std::endl;

    // Check if read transform is an affine transform.
    // If so, assign it to the affineTransform pointer.
    
    typedef itk::TransformFileReaderTemplate< double > TransformReaderType;
    typedef TransformReaderType::TransformListType* TransformListType;
    
    TransformListType transforms = transformReader->GetTransformList();
    std::cout << "  Number of transforms = " << transforms->size() << std::endl;
    
    TransformReaderType::TransformListType::const_iterator it = transforms->begin();
    if( !strcmp((*it)->GetNameOfClass(), "AffineTransform") )
    {
      affineTransform = static_cast< AffineTransformType* >( (*it).GetPointer() );
      affineTransform->Print(std::cout);
    }
  }

  // -------------------------------------------------------------
  // Initialize optimizer and registration

  // Set number of iterations and levels:
  
  registration->SetNumberOfLevels( numberOfLevels );

  if( bUseItsMinMaxStepSchedule )
  {
    MaxStepLengthLevel0 = maxStepLength * pow( 2.0, numberOfLevels - 1 );
    MinStepLengthLevel0 = minStepLength * pow( 10.0, numberOfLevels - 1 );
    if( MinStepLengthLevel0 >= MaxStepLengthLevel0 )
      MinStepLengthLevel0 = MaxStepLengthLevel0 * 0.1;
  }
  else
  {
    MaxStepLengthLevel0 = maxStepLength;
    MinStepLengthLevel0 = minStepLength;
  }

  // Setup registration filter:
  
  registration->SetFixedImage( fixedCaster->GetOutput() );
  registration->SetMovingImage( movingCaster->GetOutput() );
  registration->SetFixedImageRegion( fixedCaster->GetOutput()->GetBufferedRegion() );

  registration->SetInterpolator( interpolator );
  if( bUseNCCMetric )
  {
    NCCMetricType::Pointer metric = NCCMetricType::New();
    registration->SetMetric( metric );
  }
  else
  {
    MSDMetricType::Pointer metric = MSDMetricType::New();
    registration->SetMetric( metric );
  }

  typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
  CommandType::Pointer command = CommandType::New();
  registration->AddObserver( itk::IterationEvent(), command );

  std::cout << "  Registration parameters:" << std::endl;
  std::cout << "  -> numberOfLevels: "<<registration->GetNumberOfLevels() << std::endl;
  std::cout << "  -> numberOfIterations: "<<numberOfIterations << std::endl;
  if( bUseSimilarityReg )
    std::cout << "  -> numberOfSimilarityIterations: "<<numberOfSimilarityIterations << std::endl;

  // -------------------------------------------------------------
  // START SIMILARITY REGISTRATION
 
  if( bUseSimilarityReg )
  {
    std::cout << "  Starting similarity registration...\n" << std::endl;

    // Set optimizer scales.
    OptimizerScalesType optimizerScales( simTransform->GetNumberOfParameters() );

    optimizerScales[0] = 1.0;
    optimizerScales[1] = 1.0;
    optimizerScales[2] = 1.0;
    optimizerScales[3] = 1.0 / 800.0;
    optimizerScales[4] = 1.0 / 800.0;
    optimizerScales[5] = 1.0 / 800.0;
    optimizerScales[6] = 1.0;

    OptimizerType::Pointer optimizer = OptimizerType::New();

    optimizer->SetScales( optimizerScales );
    optimizer->SetNumberOfIterations( numberOfSimilarityIterations );
    optimizer->SetMaximumStepLength( MaxStepLengthLevel0 );
    optimizer->SetMinimumStepLength( MinStepLengthLevel0 );
    optimizer->MinimizeOn();

    // Set Observers
    globalObserverPointer = CommandIterationUpdate::New();
    optimizer->AddObserver( itk::IterationEvent(), globalObserverPointer );

    registration->SetOptimizer( optimizer );

    // Set registration filter.
    registration->SetTransform( simTransform );
    registration->SetInitialTransformParameters( simTransform->GetParameters() );
    try
    {
      registration->Update();
    }
    catch( itk::ExceptionObject & err )
    {
      std::cout << "ExceptionObject caught !" << std::endl;
      std::cout << err << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "Similarity Registration finished with condition " << optimizer->GetStopCondition() << "." << std::endl;

    // Get the final registration parameters
    OptimizerType::ParametersType finalSimParameters = registration->GetLastTransformParameters();
    if( bUseBestValuesPerIteration && globalObserverPointer.IsNotNull() )
    {
      finalSimParameters = globalObserverPointer->m_BestPosition;
      std::cout << "                   best level: " << globalObserverPointer->m_BestLevel << std::endl;
      std::cout << "best value in the last levels: " << globalObserverPointer->m_BestValue << std::endl;
      std::cout << "              best parameters: " << globalObserverPointer->m_BestPosition << std::endl;
    }

    // Generate a final transformation
    SimilarityTransformType::Pointer finalSimTransform = SimilarityTransformType::New();
    finalSimTransform->SetCenter( simTransform->GetCenter() );
    finalSimTransform->SetParameters( finalSimParameters );

    // Copy the similarity transform to the affine transform
    affineTransform->SetCenter( finalSimTransform->GetCenter() );
    affineTransform->SetMatrix( finalSimTransform->GetMatrix() );
    affineTransform->SetTranslation( finalSimTransform->GetTranslation() );
    //affineTransform->SetOffset( finalSimTransform->GetOffset() );

    //
    // Debug Output
    //
    std::cout << "Final parameters:\n" << std::endl;
    for( unsigned int i = 0; i < finalSimTransform->GetNumberOfParameters(); i++ )
    {
      std::cout << finalSimParameters[i] << ",";
    }
    std::cout << "Similarity Transform Matrix:\n" << std::endl;
    finalSimTransform->Print( std::cout );
  }

  // -------------------------------------------------------------
  // START AFFINE REGISTRATION
  
  // Set optimization scales:
  
  OptimizerScalesType optimizerScales( affineTransform->GetNumberOfParameters() );

  optimizerScales[0]  = 1.0;
  optimizerScales[1]  = 1.0;
  optimizerScales[2]  = 1.0;
  optimizerScales[3]  = 1.0;
  optimizerScales[4]  = 1.0;
  optimizerScales[5]  = 1.0;
  optimizerScales[6]  = 1.0;
  optimizerScales[7]  = 1.0;
  optimizerScales[8]  = 1.0;
  optimizerScales[9]  = 1.0 / 800.0;
  optimizerScales[10] = 1.0 / 800.0;
  optimizerScales[11] = 1.0 / 800.0;

  OptimizerType::Pointer optimizer = OptimizerType::New();
  optimizer->SetScales( optimizerScales );
  optimizer->SetNumberOfIterations( numberOfIterations );
  optimizer->SetMaximumStepLength( MaxStepLengthLevel0 );
  optimizer->SetMinimumStepLength( MinStepLengthLevel0 );
  optimizer->MinimizeOn();

  // Set observers:
  
  globalObserverPointer = CommandIterationUpdate::New();
  optimizer->AddObserver( itk::IterationEvent(), globalObserverPointer );

  registration->SetOptimizer( optimizer );

  
  // Start registration:
  
  std::cout << "  Starting affine registration ...\n" << std::endl;
  registration->SetTransform( affineTransform );
  registration->SetInitialTransformParameters( affineTransform->GetParameters() );
  try
  {
    registration->Update();
  }
  catch( itk::ExceptionObject & err )
  {
    std::cout << "ExceptionObject caught !" << std::endl;
    std::cout << err << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "  Registration finished with condition " << optimizer->GetStopCondition() << "." << std::endl;

  // Get the final registration parameters:
  
  OptimizerType::ParametersType finalAffineParameters = registration->GetLastTransformParameters();
  if( bUseBestValuesPerIteration && globalObserverPointer.IsNotNull() )
  {
    finalAffineParameters = globalObserverPointer->m_BestPosition;
    std::cout << "                   best level: " << globalObserverPointer->m_BestLevel << std::endl;
    std::cout << "best value in the last levels: " << globalObserverPointer->m_BestValue << std::endl;
    std::cout << "              best parameters: " << globalObserverPointer->m_BestPosition << std::endl;
  }

  // Generate a final transformation:
  
  AffineTransformType::Pointer finalAffineTransform = AffineTransformType::New();
  finalAffineTransform->SetCenter( affineTransform->GetCenter() );
  finalAffineTransform->SetParameters( finalAffineParameters );

  // Debug Output:
  
  std::cout << "Final affine parameters:\n" << std::endl;
  for( unsigned int i = 0; i < finalAffineTransform->GetNumberOfParameters(); i++ )
  {
    std::cout << finalAffineParameters[i] << ",";
  }
  std::cout << "Affine Transformation:\n" << std::endl;
  finalAffineTransform->Print( std::cout );

  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "WRITING output data ..."                  << std::endl;
  
  // TRANSFORM file:
  
  std::cout << "  Transform file ... "  << std::flush;
  OptimizerType::Pointer lastOptimizer = dynamic_cast< OptimizerType* >( registration->GetOptimizer() );
  
  // Writing data:
  
#if (ITK_VERSION_MAJOR == 4 && ITK_VERSION_MINOR >= 5) || ITK_VERSION_MAJOR > 4
  itk::TransformFileWriterTemplate<double>::Pointer transformWriter = itk::TransformFileWriterTemplate<double>::New();
#else
  itk::TransformFileWriter::Pointer transformWriter = itk::TransformFileWriter::New();
#endif
  //transformWriter->SetInput( affineTransform );
  transformWriter->SetInput( finalAffineTransform );
  transformWriter->SetFileName( outputTransformFilename );
  try
  {
    transformWriter->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "   ERROR while writing transformation to file." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;
  
  // WARPED MOVING image data:
  
  std::cout << "  Warped moving image file ... "  << std::endl;

  // First step: resample image:

  std::cout << "  -> Resampling output image ... " << std::endl;
  std::cout << "  -> Testing if moving image is binary image ... " << std::flush;
  
  bool movingImageIsBinaryImage = false;
  
  typedef itk::MinimumMaximumImageCalculator< ImageType > MinMaxCalulationFilterType;
  MinMaxCalulationFilterType::Pointer minMaxCalculator = MinMaxCalulationFilterType::New();
  minMaxCalculator->SetImage( movingImage );
  minMaxCalculator->Compute();
  
  typedef itk::BinaryThresholdImageFilter<ImageType, ImageType> BinaryThresholdFilterType;
  BinaryThresholdFilterType::Pointer thresholdFilter;
  if( ( minMaxCalculator->GetMinimum() == 0 ) && ( minMaxCalculator->GetMinimum() == 1 ) )
  {
    std::cout << "Yes." << std::endl;
    std::cout << "  -> Rescaling image dynamics to [0;128] ... " << std::endl;
    
    movingImageIsBinaryImage = true;
    
    thresholdFilter = BinaryThresholdFilterType::New();
    thresholdFilter->SetOutsideValue( 0 );
    thresholdFilter->SetInsideValue( 128 );
    thresholdFilter->SetLowerThreshold( 1 );
    thresholdFilter->SetUpperThreshold( 255 );
    thresholdFilter->SetInput( movingImage );
      
    try
    {
      thresholdFilter->Update();
      std::cout << "OK." << std::endl;
    }
    catch( itk::ExceptionObject &err )
    {
      std::cerr << "Error during thresholding binary input image! Exception error: " << err << std::endl;
      return EXIT_FAILURE;
    }
  }
  else
  {
    std::cout << "No. Continuing by standard resampling." << std::endl;
  }

  typedef itk::ResampleImageFilter<ImageType, ImageType> ResampleFilterType;
  ResampleFilterType::Pointer resample = ResampleFilterType::New();

  resample->SetTransform( finalAffineTransform );
  //resample->SetTransform( registration->GetOutput()->Get() );
  if( movingImageIsBinaryImage )
  {
    resample->SetInput( thresholdFilter->GetOutput() );
  }
  else
  {
    resample->SetInput( movingImage );
  }

  ImageType::PixelType backgroundGrayLevel = 0;
  resample->SetDefaultPixelValue( backgroundGrayLevel );
  /*resample->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
  resample->SetOutputOrigin( fixedImage->GetOrigin() );
  resample->SetOutputSpacing( fixedImage->GetSpacing() );
  resample->SetOutputDirection( fixedImage->GetDirection() );*/
  resample->UseReferenceImageOn();
  resample->SetReferenceImage(fixedImage);
  resample->Update();
  
  BinaryThresholdFilterType::Pointer rescalingFilter;
  if( movingImageIsBinaryImage )
  {
    std::cout << "Post-processing resampled image: rescaling dynamics to [0,1] ... "  << std::flush;
      
    rescalingFilter = BinaryThresholdFilterType::New();
    rescalingFilter->SetOutsideValue( 0 );
    rescalingFilter->SetInsideValue( 1 );
    rescalingFilter->SetLowerThreshold( 64 );
    rescalingFilter->SetUpperThreshold( 255 );
    rescalingFilter->SetInput( resample->GetOutput() );
      
    try
    {
      rescalingFilter->Update();
    }
    catch( itk::ExceptionObject &err )
    {
      std::cerr << "ERROR: Cannot execute rescaling! Exception error: " << err << std::endl;
      return EXIT_FAILURE;
    }
      
    std::cout << "OK." << std::endl;
  }

  // Writing image:
  
  if( outputImageFilename != NULL )
  {
    std::cout << "  -> Writing file ... "  << std::flush;
    ImageWriterType::Pointer imageWriter = ImageWriterType::New();
    if( movingImageIsBinaryImage )
    {
      imageWriter->SetInput( rescalingFilter->GetOutput() );
    }
    else
    {
      imageWriter->SetInput( resample->GetOutput() );
    }
    imageWriter->SetFileName( outputImageFilename );
    try
    {
      imageWriter->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "   ERROR while writing warped target image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    std::cout << "OK." << std::endl;
  }
  
  return EXIT_SUCCESS;
}
