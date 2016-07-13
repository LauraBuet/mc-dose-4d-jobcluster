/** \file icnsRobustPCADecompositionEvaluation.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2016 ICNS, UKE
 *
 ****************************************************************************/

// System includes:
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <numeric>

extern "C"
{
#include "getopt.h"
}

// ITK includes
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryErodeImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkConnectedComponentImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkLabelOverlapMeasuresImageFilter.h>
#include <itkLabelShapeKeepNObjectsImageFilter.h>
#include <itkMaskImageFilter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkThresholdImageFilter.h>

// Project includes: NONE SO FAR.

// Global typedefs:
typedef short                                                    ImagePixelType;
typedef itk::Image<ImagePixelType, 3>                            ImageType;
typedef itk::ImageFileReader<ImageType>                          ImageReaderType;
typedef itk::ImageFileWriter<ImageType>                          ImageWriterType;
typedef itk::BinaryBallStructuringElement<ImagePixelType, 3 >    StructuringElementType;
typedef itk::BinaryErodeImageFilter<ImageType, ImageType, StructuringElementType> BinaryErodeImageFilterType;
typedef itk::BinaryThresholdImageFilter<ImageType, ImageType>    BinaryThresholdFilterType;
typedef itk::ConnectedComponentImageFilter<ImageType, ImageType> ConnectedComponentImageFilterType;
typedef itk::LabelOverlapMeasuresImageFilter<ImageType>          LabelOverlapMeasuresFilterType;
typedef itk::LabelShapeKeepNObjectsImageFilter< ImageType >      KeepNObjectsImageFilterType;
typedef itk::MaskImageFilter<ImageType, ImageType>               MaskImageFilterType;
typedef itk::RescaleIntensityImageFilter<ImageType, ImageType>   RescaleFilterType;
typedef itk::ThresholdImageFilter<ImageType>                     ThresholdFilterType;


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage :\n";
  std::cout << "icnsRobustPCADecompositionEvaluation -I <list of low-rank images> -S <list of sparse images> "
               "-B <list of lesion masks> [-M <list of brain masks>] [-O <list of output files> -s <slices output images>] [...] \n";
  std::cout << "----------------------------------------------------------------------------------" << std::endl;
  std::cout << "INPUT: " << std::endl;
  //std::cout << "-I  <list of input images>           Filenames of original images.\n";
  //std::cout << "-L  <list of input images>           Filenames of low rank images.\n";
  std::cout << "-S  <list of input images>           Filenames of sparse images.\n";
  std::cout << "-G  <list of mask images>            Filenames of lesion masks (G = ground truth).\n";
  std::cout << "-M  <list of mask images>            Filenames of brain masks.\n";
  std::cout << "-t  <lower threshold>                Lower threshold applied for low-rank image binarization.\n";
  std::cout << "-r  <erosion radius>                 Radius of (currently isotropic) brain mask erosion filter.\n";
  std::cout << "----------------------------------------------------------------------------------" <<std::endl;
  std::cout << "Output: " << std::endl;
  std::cout << "-O  <list of output images>          Filenames of output lesion segmentation images.\n";
  //std::cout << "-l  <slices of output images>        Slices for 2D output images.\n ";
  std::cout << "----------------------------------------------------------------------------------" << std::endl;
  std::cout << "-h                                   Print this help.\n\n";
}


// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 1 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }

  std::cout << "================================================" << std::endl;
  std::cout << "icnsImageRPCADecomposition: Evaluation part" << std::endl;
  std::cout << "================================================" << std::endl;
  std::cout << "Reading parameters...\n" << std::endl;

  // Initializing parameters with default values:
  
  std::vector<std::string> inputSparseImageFilenames;
  std::vector<std::string> inputLesionSegmentationFilenames;
  std::vector<std::string> inputBrainSegmentationFilenames;
  std::vector<std::string> outputLesionSegmentationFilenames;
  
  std::string logFile;

  std::vector<ImageType::Pointer> sparseImages;
  std::vector<ImageType::Pointer> sparseImagesTemp;
  std::vector<ImageType::Pointer> lesionSegmentationImages;
  std::vector<ImageType::Pointer> brainSegmentationImages;
  
  ImagePixelType lowerThreshold = 100;
  ImagePixelType erosionRadius = 1;
  bool useMaskingSparseImage = false;
  bool iterateThreshold = false;
  
  // Reading parameters:
  
  char c;
  int cnt = 0;

  while( (c = getopt( argc, argv, "S:G:M:O:L:t:r:h?" )) != -1 )
  {
    switch( c )
    {
      case 'S':
        // Reading sparse rank input filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputSparseImageFilenames.push_back( argv[optind] );
          std::cout << "  SR-image [" << cnt << "]:                   " << argv[optind] << std::endl;
        }
        break;
      case 'G':
        // Reading the ground truth lesion segm filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputLesionSegmentationFilenames.push_back( argv[optind] );
          std::cout << "  Ground truth lesion segm [" << cnt << "]:   " << argv[optind] << std::endl;
        }
        break;
      case 'M':
        // Reading the barin mask filenames, as long as the next sign isn't '-'
        useMaskingSparseImage = true;
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputBrainSegmentationFilenames.push_back( argv[optind] );
          std::cout << "  Brain segmentation data [" << cnt << "]:    " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        // Reading output lesion segmentation filenames, as long as the next sign isn't '-'
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputLesionSegmentationFilenames.push_back( argv[optind] );
          std::cout << "  Out lesion segm [" << cnt << "]:            " << argv[optind] << std::endl;
        }
        break;
      case 'L':
        logFile = optarg;
        std::cout << "  Logfile:                                    " << logFile << std::endl;
        break;
      case 't':
        // Reading threshold value for segmentation of sparse rank image
        lowerThreshold = atoi( optarg );
        if( lowerThreshold == -1 ) iterateThreshold = true;
        std::cout << "  Lower threshold for SR image segmentation: " << lowerThreshold << std::endl;
        break;
      case 'r':
        erosionRadius = atoi( optarg );
        std::cout << "  Erosion radius for brain mask erosion:     " << erosionRadius << std::endl;
        break;
      case 'h':
      case '?':
        PrintHelp();
        exit( 1 );
        break;
      default:
        std::cerr << "Argument "<<(char)c<<" not processed !\n" << std::endl;
        break;
    }
  }
  
  // Check validity of arguments:
  
  if( inputSparseImageFilenames.size() < 1 )
  {
    std::cerr << "ERROR: At least one input image required!\n" << std::endl;
    return EXIT_FAILURE;
  }
  if( inputSparseImageFilenames.size() != inputLesionSegmentationFilenames.size() )
  {
    std::cerr << "ERROR: Number of sparse images and lesion segm data do not correspond!\n" << std::endl;
    return EXIT_FAILURE;
  }
  if( useMaskingSparseImage && (inputSparseImageFilenames.size() != inputBrainSegmentationFilenames.size() ) )
  {
    std::cerr << "ERROR: useMaskingSparseImage selected and number of sparse images and brain segm data do not correspond!\n" << std::endl;
    return EXIT_FAILURE;
  }

  // -------------------------------------------------------------
  // Loading input image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading sparse rank images ... " << std::endl;
  
  ImageReaderType::Pointer inputImageReader;
  
  for( unsigned int i = 0; i < inputSparseImageFilenames.size(); i++ )
  {
    std::cout << "  Reading: " << inputSparseImageFilenames[i] << std::flush;
    ImageType::Pointer inputImage;
    
    inputImageReader = ImageReaderType::New();
    inputImageReader->SetFileName( inputSparseImageFilenames[i].c_str() );
    inputImage = inputImageReader->GetOutput();
    inputImage->Update();
    
    inputImage->DisconnectPipeline();
    sparseImages.push_back( inputImage );
    sparseImagesTemp.push_back( inputImage );
    
    std::cout << " -> Ok." << std::endl;
  }
  
  std::cout << "Loading lesion segmentation images ... " << std::endl;
  for( unsigned int i = 0; i < inputLesionSegmentationFilenames.size(); i++ )
  {
    std::cout << "  Reading: " << inputLesionSegmentationFilenames[i] << std::flush;
    ImageType::Pointer inputImage;
    
    inputImageReader = ImageReaderType::New();
    inputImageReader->SetFileName( inputLesionSegmentationFilenames[i].c_str() );
    inputImage = inputImageReader->GetOutput();
    inputImage->Update();
    
    inputImage->DisconnectPipeline();
    lesionSegmentationImages.push_back( inputImage );
    
    std::cout << " -> Ok." << std::endl;
  }
  
  if( useMaskingSparseImage )
  {
    std::cout << "Loading brain mask images ... " << std::endl;
    
    for( unsigned int i = 0; i < inputBrainSegmentationFilenames.size(); i++ )
    {
      std::cout << "  Reading: " << inputBrainSegmentationFilenames[i] << std::flush;
      ImageType::Pointer inputImage;
      
      inputImageReader = ImageReaderType::New();
      inputImageReader->SetFileName( inputBrainSegmentationFilenames[i].c_str() );
      inputImage = inputImageReader->GetOutput();
      inputImage->Update();
      
      inputImage->DisconnectPipeline();
      brainSegmentationImages.push_back( inputImage );
      
      std::cout << " -> Ok." << std::endl;
    }
  }

  // -------------------------------------------------------------
  // For each sparse data set: threshold data set, compute
  // largest connected component.
  // -------------------------------------------------------------

  std::cout << "----------------------------------------" << std::endl;
  
  MaskImageFilterType::Pointer maskImageFilter;
  BinaryErodeImageFilterType::Pointer erodeFilter;
  
  if( useMaskingSparseImage )
  {
    std::cout << "Masking sparse images ... " << std::endl;
  
    for( unsigned int i = 0; i < sparseImages.size(); i++ )
    {
      std::cout << "  Image " << i << std::flush;
      
      ImageType::Pointer erodedBrainMask;
      ImageType::Pointer maskedImage;
      
      StructuringElementType structuringElement;
      structuringElement.SetRadius( erosionRadius );
      structuringElement.CreateStructuringElement();
      
      erodeFilter = BinaryErodeImageFilterType::New();
      erodeFilter->SetInput( brainSegmentationImages[i] );
      erodeFilter->SetKernel( structuringElement );
      erodeFilter->SetErodeValue( 1 );
      erodeFilter->SetBackgroundValue( 0 );
      try
      {
        erodeFilter->Update();
      }
      catch( itk::ExceptionObject &err )
      {
        std::cerr << "Cannot execute itkBinaryErodeImageFilter! Exception error: " << err << std::endl;
        return EXIT_FAILURE;
      }
      
      erodedBrainMask = erodeFilter->GetOutput();
      erodedBrainMask->Update();
      erodedBrainMask->DisconnectPipeline();
      brainSegmentationImages[i] = erodedBrainMask;
      
      maskImageFilter = MaskImageFilterType::New();
      maskImageFilter->SetInput( sparseImages[i] );
      maskImageFilter->SetMaskImage( brainSegmentationImages[i] );
      try
      {
        maskImageFilter->Update();
      }
      catch( itk::ExceptionObject &err )
      {
        std::cerr << "Cannot execute itkMaskImageFilter! Exception error: " << err << std::endl;
        return EXIT_FAILURE;
      }
      
      maskedImage = maskImageFilter->GetOutput();
      maskedImage->Update();
      maskedImage->DisconnectPipeline();
      sparseImages[i] = maskedImage;
      
      std::cout << " OK." << std::endl;
      
      // Kick out negative values:
      ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
      ImageType::Pointer thresholdedImage;
        
      thresholdFilter->ThresholdOutside( 0, itk::NumericTraits<ImageType::PixelType>::max() );
      thresholdFilter->SetOutsideValue( 0 );
      thresholdFilter->SetInput( sparseImages[i] );
      try
      {
        thresholdFilter->Update();
      }
      catch( itk::ExceptionObject &err )
      {
        std::cerr << "Cannot execute itkBinaryThresholdFilter! Exception error: " << err << std::endl;
        return EXIT_FAILURE;
      }
        
      thresholdedImage = thresholdFilter->GetOutput();
      thresholdedImage->Update();
      thresholdedImage->DisconnectPipeline();
      sparseImages[i] = thresholdedImage;
        
      std::cout << " OK." << std::endl;
      
      // Rescale image dynamics to 0 to 255:
      
      RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
      ImageType::Pointer ccImage;

      rescaleFilter->SetOutputMinimum( 0 );
      rescaleFilter->SetOutputMaximum( 255 );
      rescaleFilter->SetInput( sparseImages[i] );
        
      try
      {
        rescaleFilter->Update();
      }
      catch( itk::ExceptionObject &err )
      {
        std::cerr << "Cannot execute largest component extraction! Exception error: " << err << std::endl;
        return EXIT_FAILURE;
      }
        
      ccImage = rescaleFilter->GetOutput();
      ccImage->Update();
      ccImage->DisconnectPipeline();
      sparseImages[i] = ccImage;
    }
  }
  
  // Iterate over threshold range:
  
  std::vector< float > optimalDiceVals;
  std::vector< float > optimalFPR;
  std::vector< float > optimalFNR;
  float optimalMeanDiceVal        = 0.0;
  float optimalStdevDiceVal       = 0.0;

  ImagePixelType optimalThreshold = -1;
  unsigned int n_ThresholdValues  = 1;
  unsigned int initialThreshold   = 10;
  unsigned int incrThreshold      = 2;
  if( iterateThreshold ) n_ThresholdValues = 90;
  
  for( unsigned int currentThresholdIndex = 1; currentThresholdIndex <= n_ThresholdValues; currentThresholdIndex++ )
  {
    if( n_ThresholdValues > 1 ) lowerThreshold = (currentThresholdIndex-1)*incrThreshold + initialThreshold;
    
    std::cout << "Thresholding sparse image by theta = " << lowerThreshold << " ... " << std::endl;
    
    BinaryThresholdFilterType::Pointer thresholdFilter;
    for( unsigned int i = 0; i < sparseImages.size(); i++ )
    {
      std::cout << "  Image " << i << std::flush;
      ImageType::Pointer thresholdedImage;
      
      thresholdFilter = BinaryThresholdFilterType::New();
      thresholdFilter->SetOutsideValue( 0 );
      thresholdFilter->SetInsideValue( 1 );
      thresholdFilter->SetLowerThreshold( lowerThreshold );
      thresholdFilter->SetUpperThreshold( itk::NumericTraits<ImageType::PixelType>::max() );
      thresholdFilter->SetInput( sparseImages[i] );
      try
      {
        thresholdFilter->Update();
      }
      catch( itk::ExceptionObject &err )
      {
        std::cerr << "Cannot execute itkBinaryThresholdFilter! Exception error: " << err << std::endl;
        return EXIT_FAILURE;
      }
      
      thresholdedImage = thresholdFilter->GetOutput();
      thresholdedImage->Update();
      thresholdedImage->DisconnectPipeline();
      sparseImagesTemp[i] = thresholdedImage;
      
      std::cout << " OK." << std::endl;
    }
    
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Extracting largest connected component ... " << std::endl;
    
    ConnectedComponentImageFilterType::Pointer connected;
    KeepNObjectsImageFilterType::Pointer keepNObjectsImageFilter;
    RescaleFilterType::Pointer rescaleFilter;
    
    for( unsigned int i = 0; i < sparseImages.size(); i++ )
    {
      std::cout << "  Image " << i << ": " << std::flush;
      ImageType::Pointer ccImage;
      
      connected = ConnectedComponentImageFilterType::New ();
      connected->SetInput( sparseImagesTemp[i] );
      connected->Update();
      
      std::cout << "[no of components = " << connected->GetObjectCount() << "] " << std::flush;
      
      keepNObjectsImageFilter = KeepNObjectsImageFilterType::New();
      keepNObjectsImageFilter->SetInput( connected->GetOutput() );
      keepNObjectsImageFilter->SetBackgroundValue( 0 );
      keepNObjectsImageFilter->SetNumberOfObjects( 1 );
      keepNObjectsImageFilter->SetAttribute( KeepNObjectsImageFilterType::LabelObjectType::NUMBER_OF_PIXELS );
      
      rescaleFilter = RescaleFilterType::New();
      rescaleFilter->SetOutputMinimum( 0 );
      //rescaleFilter->SetOutputMaximum( itk::NumericTraits<ImagePixelType>::max() );
      rescaleFilter->SetOutputMaximum( 1 );
      rescaleFilter->SetInput( keepNObjectsImageFilter->GetOutput() );
      
      try
      {
        rescaleFilter->Update();
      }
      catch( itk::ExceptionObject &err )
      {
        std::cerr << "Cannot execute largest component extraction! Exception error: " << err << std::endl;
        return EXIT_FAILURE;
      }
      
      ccImage = rescaleFilter->GetOutput();
      ccImage->Update();
      ccImage->DisconnectPipeline();
      sparseImagesTemp[i] = ccImage;
      
      std::cout << "OK. " << std::endl;
    }
    
    // -------------------------------------------------------------
    // Compute Dice between ccImage and ground truth.
    // -------------------------------------------------------------
    
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Computing DICE coefs ... " << std::endl;
    
    std::vector<float> diceCoefs;
    std::vector<float> falsePositiveRate;
    std::vector<float> falseNegativeRate;
    LabelOverlapMeasuresFilterType::Pointer labelOverlapMeasuresFilter;
    
    for( unsigned int i = 0; i < sparseImagesTemp.size(); i++ )
    {
      std::cout << "  Image [" << i << "]: " << std::endl;
      
      labelOverlapMeasuresFilter = LabelOverlapMeasuresFilterType::New();
      labelOverlapMeasuresFilter->SetSourceImage( sparseImagesTemp[i] );
      labelOverlapMeasuresFilter->SetTargetImage( lesionSegmentationImages[i] );
      labelOverlapMeasuresFilter->Update();
      
      LabelOverlapMeasuresFilterType::MapType labelMap = labelOverlapMeasuresFilter->GetLabelSetMeasures();
      LabelOverlapMeasuresFilterType::MapType::const_iterator it;
      
      for( it = labelMap.begin(); it != labelMap.end(); ++it )
      {
        if( (*it).first == 0 ) continue;
        
        int label = (*it).first;
        
        std::cout << "  -> LABEL " << label << ": " << std::flush;
        //std::cout << labelOverlapMeasuresFilter->GetTargetOverlap( label ) << std::endl; //
        std::cout << labelOverlapMeasuresFilter->GetDiceCoefficient( label ) << std::endl;
      }
      diceCoefs.push_back( labelOverlapMeasuresFilter->GetDiceCoefficient( 1 ) );
      falsePositiveRate.push_back( labelOverlapMeasuresFilter->GetFalsePositiveError(1) );
      falseNegativeRate.push_back( labelOverlapMeasuresFilter->GetFalseNegativeError(1) );
    }
    
    float sum = std::accumulate( diceCoefs.begin(), diceCoefs.end(), 0.0 );
    float mean = sum/diceCoefs.size();
    
    float sq_sum = std::inner_product( diceCoefs.begin(), diceCoefs.end(), diceCoefs.begin(), 0.0 );
    float stdev = std::sqrt( sq_sum / diceCoefs.size() - mean * mean );
    
    std::cout << "  TOTAL: " << mean << " +/- " << stdev << std::endl;
    
    if( mean > optimalMeanDiceVal )
    {
      optimalMeanDiceVal  = mean;
      optimalStdevDiceVal = stdev;
      optimalDiceVals     = diceCoefs;
      optimalThreshold    = lowerThreshold;
      optimalFPR          = falsePositiveRate;
      optimalFNR          = falseNegativeRate;
    }
  }
  
  std::cout << "+++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "BEST TOTAL: " << optimalMeanDiceVal << " +/- " << optimalStdevDiceVal << std::endl;
  std::cout << "THRESHOLD:  " << optimalThreshold << std::endl;
  std::cout << "+++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "DICE" << std::endl;
  
  for( unsigned int i = 0; i < sparseImagesTemp.size(); i++ )
  {
    std::cout << optimalDiceVals[i] << std::endl;
  }
  
  std::cout << "+++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "SPC" << std::endl;
  
  for( unsigned int i = 0; i < sparseImagesTemp.size(); i++ )
  {
    std::cout << 1-optimalFPR[i] << std::endl;
  }
  
  std::cout << "+++++++++++++++++++++++++++++++" << std::endl;
  std::cout << "SENS" << std::endl;
  
  for( unsigned int i = 0; i < sparseImagesTemp.size(); i++ )
  {
    std::cout << 1-optimalFNR[i] << std::endl;
  }

  // -------------------------------------------------------------
  // Writing output data:
  // -------------------------------------------------------------
  
  //std::cout << "----------------------------------------" << std::endl;
  //std::cout << "Writing output image ... "                << std::endl;
  
  ImageWriterType::Pointer imageWriter;
  for( unsigned int i = 0; i < sparseImages.size(); i++ )
  {
    //std::cout << "  Image [" << i << "] ... " << std::flush;
    
    imageWriter = ImageWriterType::New();
    imageWriter->SetInput( sparseImagesTemp[i] );
    //imageWriter->SetInput( brainSegmentationImages[i] );
    imageWriter->SetFileName( outputLesionSegmentationFilenames[i] );
    try
    {
      imageWriter->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "ERROR while writing output lesion segm image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    //std::cout << "OK." << std::endl;
  }
  
  

  // TODO: write logfile
  
  return EXIT_SUCCESS;
}
