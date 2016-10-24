/** \file icnsVectorFieldStatistics.cxx
 *
 *  \b Initial \b Author: Rene Werner\n\n
 *  \b Copyright (C) 2016 IPMI at ICNS at UKE
 *
 ****************************************************************************/

// System includes:
#include <string>
#include <fstream>

// ITK includes:
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"

// Project includes: None so far

// Global typedefs:

typedef itk::Image<short, 3>                        ImageType;
typedef itk::ImageFileReader<ImageType>             ImageReaderType;
typedef itk::Image<itk::Vector<float, 3>, 3>        DisplacementFieldType;
typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;

// ---------------------------------------------------------------
// Declaration of additional routines:
// ---------------------------------------------------------------

void ComputeVectorFieldMagnitudeStatistics( const DisplacementFieldType::Pointer field,
                                            std::ofstream& logfile,
                                            const ImageType::Pointer mask,
                                            const bool useMask );
void ComputeVectorFieldComponentStatistics( const DisplacementFieldType::Pointer field,
                                            std::ofstream& logfile,
                                            const ImageType::Pointer mask,
                                            const bool useMask,
                                            int component );

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "icnsVectorFieldStatistics -I <vector field FN> -S <mask image FN> -L <logfile FN> [...]\n";
  std::cout << std::endl;
  std::cout << "-I <image FN>              [IN] Vector field to analyze." << std::endl;
  std::cout << "-S <mask image FN>         [IN] Mask image to constrain analysis." << std::endl;
  std::cout << "-F <logfile FN>            [OUT] Filename of logfile to write." << std::endl;
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
  std::cout << "icnsVectorFieldStatistics"                  << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Reading parameters and preparing data ..." << std::endl;

  // Input / output filenames:
  char* vectorFieldFilename = NULL;
  char* maskImageFilename   = NULL;
  char* logfileFilename     = NULL;
  
  // Other params:
  bool useMask              = false;
  
  // Reading data:
  
  int c;
  while( (c = getopt( argc, argv, "I:S:L:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        vectorFieldFilename = optarg;
        std::cout << "  Vector field filename:     " << vectorFieldFilename << std::endl;
        break;
      case 'M':
      case 'S':
        maskImageFilename = optarg;
        useMask = true;
        std::cout << "  Mask image filename:       " << maskImageFilename << std::endl;
        break;
      case 'L':
        logfileFilename = optarg;
        std::cout << "  Logfile filename:          " << logfileFilename << std::endl;
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_SUCCESS;
      default:
        std::cout << "ERROR: Argument " << (char) c << " not processed!\n" << std::endl;
        return EXIT_FAILURE;
    }
  }

  // Check validity of arguments:
  
  if( vectorFieldFilename == NULL )
  {
    std::cerr << "No vector field given!" << std::endl;
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  if( logfileFilename == NULL )
  {
    std::cerr << "No logfile given!" << std::endl;
    return EXIT_FAILURE;
  }

  // -------------------------------------------------------------
  // Declaring image pointers and loading input images:
  
  DisplacementFieldType::Pointer vectorField;
  ImageType::Pointer maskImage;
  
   std::cout << "------------------------------------------" << std::endl;
  std::cout << "Loading input data ..." << std::endl;
  
  // Loading input vector field:
  
  std::cout << "  Loading vector field ..." << std::flush;
  DisplacementFieldReaderType::Pointer vectorFieldReader = DisplacementFieldReaderType::New();
  vectorFieldReader->SetFileName( vectorFieldFilename );
  try
  {
    vectorFieldReader->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "ERROR while loading vector field." << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  vectorField = vectorFieldReader->GetOutput();
  vectorField->Update();
  std::cout << "OK." << std::endl;
  
  // Loading mask image:
  
  if( useMask )
  {
    std::cout << "  Loading mask image ..." << std::flush;
    ImageReaderType::Pointer maskImageReader = ImageReaderType::New();
    maskImageReader->SetFileName( maskImageFilename );
    try
    {
      maskImageReader->Update();
    }
    catch( itk::ExceptionObject& excp )
    {
      std::cerr << "  ERROR while loading reference image." << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    maskImage = maskImageReader->GetOutput();
    maskImage->Update();
    std::cout << "OK." << std::endl;
  }

  // Generate logfile object:
  
  std::ofstream logfile;
  logfile.open( logfileFilename, std::ios::app );
      
  logfile << "ICNS VECTOR FIELD STATISTICS" << std::endl;
  logfile << "Vector field: " << vectorFieldFilename << std::endl;
  if( useMask ) logfile << "Mask image:   " << maskImageFilename << std::endl;
  
  // Compute displacement field magnitude statistics:
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Computing vector field magnitude statistics ..." << std::flush;
  ComputeVectorFieldMagnitudeStatistics( vectorField, logfile, maskImage, useMask );
  std::cout << "OK." << std::endl;
  
  // Compute displacement field components statistics:
  
  std::cout << "Computing vector field component 0 statistics ..." << std::flush;
  ComputeVectorFieldComponentStatistics( vectorField, logfile, maskImage, useMask, 0 );
  std::cout << "OK." << std::endl;
  
  std::cout << "Computing vector field component 1 statistics ..." << std::flush;
  ComputeVectorFieldComponentStatistics( vectorField, logfile, maskImage, useMask, 1 );
  std::cout << "OK." << std::endl;
  
  std::cout << "Computing vector field component 2 statistics ..." << std::flush;
  ComputeVectorFieldComponentStatistics( vectorField, logfile, maskImage, useMask, 2 );
  std::cout << "OK." << std::endl;
  
  // Show logfile content and close logfile:
  
  logfile.close();
  
  std::ifstream writtenLogfile( logfileFilename );
  std::string str;
  std::string fileContents;
  while (std::getline( writtenLogfile, str ))
  {
    fileContents += str;
    fileContents.push_back('\n');
  }
  std::cout << fileContents;

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Finished." << std::endl;
  std::cout << "==========================================" << std::endl;

  return EXIT_SUCCESS;
}


// Implementation of additional routines:

void ComputeVectorFieldMagnitudeStatistics( const DisplacementFieldType::Pointer field,
                                           std::ofstream& logfile,
                                           const ImageType::Pointer mask,
                                           const bool useMask )
{
  // Prepare result variables:
  unsigned int nVoxels          = 0;
  double averageVectorMagnitude;// = 0.0;
  double stdVectorMagnitude;//     = 0.0;
  double minVectorMagnitude;//    = 0.0;
  double maxVectorMagnitude;//     = 0.0;
  
  // Generate field and image iterators:
  typedef itk::ImageRegionConstIterator< DisplacementFieldType > ConstFieldIteratorType;
  typedef itk::ImageRegionConstIterator< ImageType > ConstMaskIteratorType;
  
  ConstFieldIteratorType fieldIt( field, field->GetLargestPossibleRegion() );
  ConstMaskIteratorType maskIt;
  if( useMask )
  {
    maskIt = ConstMaskIteratorType( mask, mask->GetLargestPossibleRegion() );
    maskIt.GoToBegin();
  }

  DisplacementFieldType::PixelType voxelValue;
  for( fieldIt.GoToBegin(); !fieldIt.IsAtEnd(); ++fieldIt)
  {
    // If mask is given and pixel is outside mask, continue:
    if( useMask && maskIt.Get() == 0 )
    {
      ++maskIt;
      continue;
    }
    
    // Otherwise update statistics:
    nVoxels++;
    voxelValue             = fieldIt.Get();
    double vectorMagnitude = sqrt( voxelValue[0]*voxelValue[0] + voxelValue[1]*voxelValue[1] + voxelValue[2]*voxelValue[2] );
    averageVectorMagnitude += vectorMagnitude;
    stdVectorMagnitude     += vectorMagnitude*vectorMagnitude;
    minVectorMagnitude     = (minVectorMagnitude < vectorMagnitude) ? minVectorMagnitude : vectorMagnitude;
    maxVectorMagnitude     = (maxVectorMagnitude > vectorMagnitude) ? maxVectorMagnitude : vectorMagnitude;
    
    // And update mask iterator:
    if( useMask ) ++maskIt;
  }
  
  // Calculate average values:
  averageVectorMagnitude /= nVoxels; // E(X)
  stdVectorMagnitude     /= nVoxels; // E(X^2)
  stdVectorMagnitude     = sqrt( stdVectorMagnitude - averageVectorMagnitude*averageVectorMagnitude ); // V(X) = E(X^2) - E(X)^2
  
  // Write to logfile:
  logfile << "Vector field magnitude [mm]: " << averageVectorMagnitude << " pm " << stdVectorMagnitude;
  logfile << " [min: " << minVectorMagnitude << "; max: " << maxVectorMagnitude << "]" << std::endl;
}

void ComputeVectorFieldComponentStatistics( const DisplacementFieldType::Pointer field,
                                            std::ofstream& logfile,
                                            const ImageType::Pointer mask,
                                            const bool useMask,
                                            int component )
{
  // Prepare result variables:
  unsigned int nVoxels = 0;
  double averageValue;//  = 0.0;
  double stdValue;//      = 0.0;
  double minValue;//      = 0.0;
  double maxValue;//      = 0.0;
  
  // Generate field and image iterators:
  typedef itk::ImageRegionConstIterator< DisplacementFieldType > ConstFieldIteratorType;
  typedef itk::ImageRegionConstIterator< ImageType > ConstMaskIteratorType;
  
  ConstFieldIteratorType fieldIt( field, field->GetLargestPossibleRegion() );
  ConstMaskIteratorType maskIt;
  if( useMask )
  {
    maskIt = ConstMaskIteratorType( mask, mask->GetLargestPossibleRegion() );
    maskIt.GoToBegin();
  }
  
  double voxelValue = 0.0;
  for( fieldIt.GoToBegin(); !fieldIt.IsAtEnd(); ++fieldIt)
  {
    // If mask is given and pixel is outside mask, continue:
    if( useMask && maskIt.Get() == 0 )
    {
      ++maskIt;
      continue;
    }
    
    // Otherwise update statistics:
    nVoxels++;
    voxelValue   = fieldIt.Get()[component];
    averageValue += voxelValue;
    stdValue     += voxelValue*voxelValue;
    minValue     = (minValue < voxelValue) ? minValue : voxelValue;
    maxValue     = (maxValue > voxelValue) ? maxValue : voxelValue;
    
    // And update mask iterator:
    if( useMask ) ++maskIt;
  }
  
  // Calculate average values:
  averageValue /= nVoxels; // E(X)
  stdValue     /= nVoxels; // E(X^2)
  stdValue     = sqrt( stdValue - averageValue*averageValue ); // V(X) = E(X^2) - E(X)^2
  
  // Write to logfile:
  logfile << "Vector field comp. " << component << " [mm]:   " << averageValue << " pm " << stdValue;
  logfile << " [min: " << minValue << "; max: " << maxValue << "]" << std::endl;
}
