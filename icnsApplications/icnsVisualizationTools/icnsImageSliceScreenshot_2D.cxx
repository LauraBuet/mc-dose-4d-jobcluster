/** \file icnsImageSliceScreenshot_2D.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2016 Department of Computational Neuroscience,
 *     University Medical Center Hamburg-Eppendorf
 *
 ****************************************************************************/

// System includes:
#include <iostream>
extern "C"
{
#include "getopt.h"
}

// ITK includes:
#include <itkExtractImageFilter.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkChangeInformationImageFilter.h>

// VTK includes:
#include <itkImageToVTKImageFilter.h>

#include "vtkCamera.h"
#include "vtkVersion.h"
#include "vtkImageViewer.h"
#include "vtkImageMapper3D.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkImageActor.h"
#include "vtkInteractorStyleImage.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkWindowLevelLookupTable.h"
#include "vtkImageMapToColors.h"
#include "vtkImageMapToWindowLevelColors.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"

// Project includes: NONE so far.

// Global typedefs:
typedef itk::Image<short, 2>                          ImageType;
typedef itk::ImageFileReader<ImageType>               ImageReaderType;
typedef itk::ImageFileWriter<ImageType>               ImageWriterType;
typedef itk::ExtractImageFilter<ImageType, ImageType> ExtractImageFilterType;
typedef itk::ChangeInformationImageFilter<ImageType>  ChangeImageInformationFilterType;
typedef itk::ImageToVTKImageFilter<ImageType>         ITKToVTKConnectorType;


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage:\n";
  std::cout << "icnsImageSliceScreenshot_2D -I <input image> -O <output image> "
              "-w <vis. window> -l <vis. level>\n\n";
  //            "-w <vis. window> -l <vis. level> [-i <x-index> <y-index> <z-index>]\n\n";
    
  std::cout << "-I <input image>               Filename of input image to be visualized.\n";
  std::cout << "-O <output image>              Filename of output screenshot image.\n";
  std::cout << "-l <vis. level>                Level (= window center) used for visualization.\n";
  std::cout << "-w <vis. window>               window (= window width) used for visualization.\n";
  // std::cout << "-i <x> <y> <z>                 Starting index of cropping region.\n";
  // std::cout << "-s <x> <y> <z>                 Size of cropping region.\n";
  
  std::cout << "-h                             Print this help.\n";
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
  
  std::cout << "========================================" << std::endl;
  std::cout << "icnsImageSliceScreenshot_2D" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "Reading parameters ... " << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  
  char* inputImageFilename  = NULL;
  char* outputImageFilename = NULL;
  
  float visualizationLevel = 0.0;
  float visualizationWindow = 0.0;
  //ImageType::IndexType croppingRegionStartingIndex;
  //croppingRegionStartingIndex.Fill( 0 );
  //ImageType::SizeType croppingRegionSize;
  //croppingRegionSize.Fill( 0 );
  
  // Reading parameters: 
  while( (c = getopt( argc, argv, "I:O:l:w:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputImageFilename = optarg;
        std::cout << "  Image to visualize:                " << inputImageFilename << std::endl;
        break;
      case 'O':
        outputImageFilename = optarg;
        std::cout << "  Screenshot filename:               " << outputImageFilename << std::endl;
        break;
      case 'l':
        visualizationLevel = atof( optarg );
        std::cout << "  Viualization level:                " << visualizationLevel << std::endl;
        break;
      case 'w':
        visualizationWindow = atof( optarg );
        std::cout << "  Viualization window:               " << visualizationWindow << std::endl;
        break;
/*      case 'i':
#ifdef __APPLE__
        std::cout << "  -- getopt with APPLE system --" << std::endl;
        croppingRegionStartingIndex[0] = atoi(argv[optind-1]);
        croppingRegionStartingIndex[1] = atoi(argv[optind]);
        croppingRegionStartingIndex[2] = atoi(argv[optind+1]);
        optind = optind+2;
#elif __unix
        std::cout << "  -- getopt with unix system --" << std::endl;
        croppingRegionStartingIndex[0] = atof(argv[optind]);
        croppingRegionStartingIndex[1] = atof(argv[optind+1]);
        croppingRegionStartingIndex[2] = atof(argv[optind+2]);
        optind = optind+3;
#endif
        std::cout << "  Cropping region starting index:    " << croppingRegionStartingIndex << std::endl;
        break;
      case 's':
#ifdef __APPLE__
        std::cout << "  -- getopt with APPLE system --" << std::endl;
        croppingRegionSize[0] = atoi(argv[optind-1]);
        croppingRegionSize[1] = atoi(argv[optind]);
        croppingRegionSize[2] = atoi(argv[optind+1]);
        optind = optind+2;
#elif __unix
        std::cout << "  -- getopt with unix system --" << std::endl;
        croppingRegionSize[0] = atof(argv[optind]);
        croppingRegionSize[1] = atof(argv[optind+1]);
        croppingRegionSize[2] = atof(argv[optind+2]);
        optind = optind+3;
#endif
        std::cout << "  Cropping region size:              " << croppingRegionSize << std::endl;
        break;
 */
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_FAILURE;
        break;
      default:
        std::cout << "  Argument " << (char)c << " not processed!\n" << std::endl;
        return EXIT_FAILURE;
    }
  }
  
  // Check validity of arguments:
  
  if( inputImageFilename == NULL )
  {
    std::cerr << "No input image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputImageFilename == NULL )
  {
    std::cerr << "No output image filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input image ... " << std::flush;

  ImageReaderType::Pointer inputImageReader = ImageReaderType::New();
  inputImageReader->SetFileName( inputImageFilename );
  try
  {
    inputImageReader->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "  ERROR while loading input image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Visualizing image:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Visualizing image:" << std::endl;
  
  // Converting ITK to VTK image:
  
  ITKToVTKConnectorType::Pointer originalConnector = ITKToVTKConnectorType::New();
  originalConnector->SetInput( inputImageReader->GetOutput());
  
  // Generate LUT for Window/Level-acc. visualization:
  
  double lowerVisualizationRangeBound = -1;
  double upperVisualizationRangeBound = -1;
  if( visualizationLevel == -1 || visualizationWindow == -1 )
  {
    std::cout << "  No specific or invalid window/level values specified.\n  Display image with max. intensity range." << std::endl;
    lowerVisualizationRangeBound = originalConnector->GetOutput()->GetScalarRange()[0];
    upperVisualizationRangeBound = originalConnector->GetOutput()->GetScalarRange()[1];
  }
  else
  {
    lowerVisualizationRangeBound = visualizationLevel - 0.5*visualizationWindow;
    upperVisualizationRangeBound = visualizationLevel + 0.5*visualizationWindow;
  }
  
  std::cout << "Range: " << originalConnector->GetOutput()->GetScalarRange() << std::endl;
  
  vtkSmartPointer<vtkLookupTable> currentLUT = vtkSmartPointer<vtkLookupTable>::New();
  currentLUT->SetRange( lowerVisualizationRangeBound, upperVisualizationRangeBound ); // image intensity range
  currentLUT->SetValueRange(0.0, 1.0);                                                // from black to white
  currentLUT->SetSaturationRange(0.0, 0.0);                                           // no color saturation
  currentLUT->SetRampToLinear();
  currentLUT->Build();
  
  // Create mapper and actor:
  
  vtkSmartPointer<vtkImageMapToWindowLevelColors> imageMapper = vtkImageMapToWindowLevelColors::New();
  imageMapper->SetLookupTable( currentLUT );
  imageMapper->SetInputData( originalConnector->GetOutput() );
  imageMapper->Update();
  
  vtkSmartPointer<vtkImageActor> originalActor = vtkSmartPointer<vtkImageActor>::New();
  originalConnector->Update();
  originalActor->GetMapper()->SetInputConnection( imageMapper->GetOutputPort() );
  originalActor->InterpolateOff();
  
  // Visualize data:
  
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);
  
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(renderer);
  renderer->AddActor(originalActor);
  
  renderer->ResetCamera();
  
  // Set up camera to fill the renderer with the image
  double origin[3];
  double spacing[3];
  int extent[6];
  originalConnector->GetOutput()->GetOrigin( origin );
  originalConnector->GetOutput()->GetSpacing( spacing );
  originalConnector->GetOutput()->GetExtent( extent );
  
  vtkCamera* camera = renderer->GetActiveCamera();
  camera->ParallelProjectionOn();
  
  double xc = origin[0] + 0.5*(extent[0] + extent[1])*spacing[0];
  double yc = origin[1] + 0.5*(extent[2] + extent[3])*spacing[1];
  //double xd = (extent[1] - extent[0] + 1)*spacing[0];
  double yd = (extent[3] - extent[2] + 1)*spacing[1];
  double d = camera->GetDistance();
  camera->SetParallelScale(0.5*yd-0.5);
  camera->SetFocalPoint(xc,yc,0.0);
  camera->SetPosition(xc,yc,d);
  
  // Render again to set the correct view
  
  renderWindow->Render();
  
  vtkSmartPointer<vtkInteractorStyleImage> style =
  vtkSmartPointer<vtkInteractorStyleImage>::New();
  interactor->SetInteractorStyle(style);
  
  // Screenshot
  
  vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter =
  vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowToImageFilter->SetInput(renderWindow);
  //windowToImageFilter->SetMagnification(3); //set the resolution of the output image (3 times the current resolution of vtk render window)
  //windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
  windowToImageFilter->ReadFrontBufferOff(); // read from the back buffer
  windowToImageFilter->Update();
  
  vtkSmartPointer<vtkPNGWriter> writer =
  vtkSmartPointer<vtkPNGWriter>::New();
  writer->SetFileName( outputImageFilename );
  writer->SetInputConnection(windowToImageFilter->GetOutputPort());
  writer->Write();
  
  //interactor->Start();
  
/*  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output image ... "                << std::flush;
  
  ImageWriterType::Pointer imageWriter = ImageWriterType::New();
  imageWriter->SetInput( changeImageInformationFilter->GetOutput() );
  imageWriter->SetFileName( outputImageFilename );
  try
  {
    imageWriter->Update();
  }
  catch( itk::ExceptionObject& excp )
  {
    std::cerr << "ERROR while writing warped target image." << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "OK." << std::endl;*/
  
  return EXIT_SUCCESS;
}
