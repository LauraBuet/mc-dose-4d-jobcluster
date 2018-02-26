/** \file icnsComposeVectorFieldS.cxx
 *
 *  \b Initial \b Author: Rene Werner\n\n
 *  \b Copyright (C) 2016 IPMI at ICNS at UKE
 *
 ****************************************************************************/

// System includes:
#include <string>
#include <fstream>

// ITK includes:
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkComposeDisplacementFieldsImageFilter.h>

// Project includes: None so far

// Global typedefs:

typedef itk::Image<short, 3>                        ImageType;
typedef itk::ImageFileReader<ImageType>             ImageReaderType;
typedef itk::Image<itk::Vector<double, 3>, 3>       DisplacementFieldType;
typedef itk::ImageFileReader<DisplacementFieldType> DisplacementFieldReaderType;
typedef itk::ImageFileWriter<DisplacementFieldType> DisplacementFieldWriterType;

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "icnsComposeVectorFields -I <vector field FN> -S <vector field FN> -O <vector field FN> [...]\n";
  std::cout << std::endl;
  std::cout << "-I <image FN>              [IN] Vector field I." << std::endl;
  std::cout << "-S <image FN>              [IN] Vector field II." << std::endl;
  std::cout << "-O <image FN>              [OUT] Filename of output vector field." << std::endl;
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
  
  std::cout << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsComposeVectorFields"                  << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Reading parameters and preparing data ..." << std::endl;

  // Input / output filenames:
  char* vectorFieldFilename = NULL;
  char* vectorFieldFilenameII   = NULL;
  char* outputVectorFieldFilename     = NULL;
  
  // Reading data:
  
  int c;
  while( (c = getopt( argc, argv, "I:S:O:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        vectorFieldFilename = optarg;
        std::cout << "  Vector field filename:     " << vectorFieldFilename << std::endl;
        break;
      case 'S':
        vectorFieldFilenameII = optarg;
        std::cout << "  Vector field filename:     " << vectorFieldFilenameII << std::endl;
        break;
      case 'O':
        outputVectorFieldFilename = optarg;
        std::cout << "  Output vector field filename:          " << outputVectorFieldFilename << std::endl;
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
    std::cerr << "No vector field I given!" << std::endl;
    PrintHelp();
    return EXIT_FAILURE;
  }

  if( vectorFieldFilenameII == NULL )
  {
    std::cerr << "No vector field II given!" << std::endl;
    PrintHelp();
    return EXIT_FAILURE;
  }

  // -------------------------------------------------------------
  // Declaring image pointers and loading input images:
  
  DisplacementFieldType::Pointer vectorField;
  DisplacementFieldType::Pointer vectorFieldII;
  DisplacementFieldType::Pointer composedDisplacementField;

  
   std::cout << "------------------------------------------" << std::endl;
  std::cout << "Loading input data ..." << std::endl;
  
  // Loading input vector field:
  
  std::cout << "  Loading vector field I ..." << std::flush;
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
  vectorFieldReader->Update();
  vectorField = vectorFieldReader->GetOutput();
  vectorField->Update();
  std::cout << "OK." << std::endl;
  
  // Loading input vector field:

  std::cout << "  Loading vector field II ..." << std::flush;
  DisplacementFieldReaderType::Pointer vectorFieldReaderII = DisplacementFieldReaderType::New();
  vectorFieldReaderII->SetFileName( vectorFieldFilenameII );
  try
  {
    vectorFieldReaderII->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "ERROR while loading vector field." << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  vectorFieldReaderII->Update();
  vectorFieldII = vectorFieldReaderII->GetOutput();
  vectorFieldII->Update();
  std::cout << "OK." << std::endl;

  std::cout << "  Compose vector fields ..." << std::flush;

  typedef itk::ComposeDisplacementFieldsImageFilter<DisplacementFieldType,DisplacementFieldType> ComposeDisplacementFieldsImageFilterType;

  ComposeDisplacementFieldsImageFilterType::Pointer composeFilter = ComposeDisplacementFieldsImageFilterType::New();

  composeFilter->SetDisplacementField(vectorField);
  composeFilter->SetWarpingField(vectorFieldII);
  composeFilter->Update();

  composedDisplacementField = composeFilter->GetOutput();
  composedDisplacementField->Update();
  std::cout << "OK" << std::endl;


  if( outputVectorFieldFilename != NULL)
  {
    std::cout << "Saving composed displacement field..." << std::endl;
    DisplacementFieldWriterType::Pointer DisplacementFieldWriter;
    DisplacementFieldWriter = DisplacementFieldWriterType::New();

    DisplacementFieldWriter->SetInput( composedDisplacementField );
    DisplacementFieldWriter->SetFileName( outputVectorFieldFilename );
    DisplacementFieldWriter->Update();
  }

  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Finished." << std::endl;
  std::cout << "==========================================" << std::endl;

  return EXIT_SUCCESS;

}
