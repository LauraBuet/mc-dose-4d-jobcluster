/** \file icnsICPBasedSurfaceRegistration.cxx
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

// ITK includes: None so far.

// VTK includes:


// XXX
#include <vtkActor.h>
#include <vtkLandmarkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkVector.h>

#include <vtkInteractorStyleSwitch.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkInteractorStyleTrackball.h>

#include <vtkCallbackCommand.h>

#include <vtkObjectFactory.h>
#include <vtkLookupTable.h>
#include <vtkOBJReader.h>
#include <vtkPoints.h>
#include <vtkPointPicker.h>
#include <vtkProperty.h>
#include <vtkRendererCollection.h>
#include <vtkSphereSource.h>


// Project includes:
#include "vtkTexturingHelper.h"

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage :\n";
  std::cout << "icnsLandmarkBasedSurfaceRegistration -F <fixed vtk mesh> -M <moving vtk mesh> [...]\n\n";
    
  std::cout << "-F <fixed vtk mesh>            Filename of fixed mesh.\n";
  std::cout << "-M <moving vtk mesh>           Filename of moving mesh.\n";
  std::cout << "-L <fixed mesh landmarks>      Filename of fixed landmarks.\n";
  std::cout << "-N <moving mesh landmarks>     Filename of fixed landmarks.\n";
  
  std::cout << "-t <texture file prefix>       Prefix of texture files.\n";
  std::cout << "-e <texture file postfix>      Postfix of texture files (e.g. .png).\n";
  std::cout << "-n <number of texture files>   Number of texture files.\n";
  
  std::cout << "-h                             Print this help.\n";
  std::cout << "\n";
}

// ---------------------------------------------------------------
// Additional routine declarations:
// ---------------------------------------------------------------

static void ClickCallbackFunction ( vtkObject* caller,
                                   long unsigned int eventId,
                                   void* clientData,
                                   void* callData );

void Write4x4Matrix( vtkMatrix4x4* matrix, std::string filename );
vtkMatrix4x4* Load4x4Matrix( std::string filename );

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
  
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsLandmarkBasedSurfaceRegistration" << std::endl;
  std::cout << "==========================================" << std::endl;
  
  // -------------------------------------------------------------
  // Initializing parameters with default values:
  
  char* inputFixedMeshLandmarksFilename  = NULL;
  char* inputMovingMeshLandmarksFilename = NULL;
  char* inputTEMPMeshFilename          = NULL;
  
  char* inputTexturePrefix               = NULL;
  char* inputTexturePostfix              = NULL;
  unsigned int nTextureFiles             = 0;
  
  std::vector< std::string > inputMeshFilenames;
  std::vector< std::string > inputTransformFilenames;
  std::vector< std::string > outputMeshFilenames;
  std::vector< std::string > outputTransformFilenames;
  
  // Reading parameters:
  
  int c;
  int cnt = 0;
  
  while( (c = getopt( argc, argv, "I:J:O:T:L:N:Z:t:e:n:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          if( cnt == 0 )
          {
            inputMeshFilenames.push_back( argv[optind] );
            std::cout << "Input fixed mesh:  " << inputMeshFilenames[0] << std::endl;
          }
          if( cnt == 1 )
          {
            inputMeshFilenames.push_back( argv[optind] );
            std::cout << "Input moving mesh: " << inputMeshFilenames[1] << std::endl;
          }
          if( cnt > 1 ) std::cout << "WARNING: >2 input mesh filenames specified. Ignoring entries." << std::endl;
        }
        break;
      case 'J':
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          if( cnt == 0 )
          {
            inputTransformFilenames.push_back( argv[optind] );
            std::cout << "Initial input fixed transform:  " << inputTransformFilenames[0] << std::endl;
          }
          if( cnt == 1 )
          {
            inputTransformFilenames.push_back( argv[optind] );
            std::cout << "Initial input moving mesh: " << inputTransformFilenames[1] << std::endl;
          }
          if( cnt > 1 ) std::cout << "WARNING: >2 input transform filenames specified. Ignoring entries." << std::endl;
        }
        break;
      case 'Z':
        inputTEMPMeshFilename = optarg;
        std::cout << "TEMP mesh filename:             " << inputTEMPMeshFilename << std::endl;
        break;
      case 'O':
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          if( cnt == 0 )
          {
            outputMeshFilenames.push_back( argv[optind] );
            std::cout << "Output fixed mesh: " << outputMeshFilenames[0] << std::endl;
          }
          if( cnt == 1 )
          {
            outputMeshFilenames.push_back( argv[optind] );
            std::cout << "Output moving mesh: " << outputMeshFilenames[1] << std::endl;
          }
          if( cnt > 1 ) std::cout << "WARNING: >2 output mesh filenames specified. Ignoring entries." << std::endl;
        }
        break;
      case 'T':
        optind--;
        cnt = 0;
        std::cout << std::endl;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          if( cnt == 0 )
          {
            outputTransformFilenames.push_back( argv[optind] );
            std::cout << "Output fixed transform: " << outputTransformFilenames[0] << std::endl;
          }
          if( cnt == 1 )
          {
            outputTransformFilenames.push_back( argv[optind] );
            std::cout << "Output moving transform: " << outputTransformFilenames[1] << std::endl;
          }
          if( cnt > 1 ) std::cout << "WARNING: >2 output transform filenames specified. Ignoring entries." << std::endl;
        }
        break;
      case 'L':
        inputFixedMeshLandmarksFilename = optarg;
        std::cout << "Fixed mesh landmarks:           " << inputFixedMeshLandmarksFilename << std::endl;
        break;
      case 'N':
        inputMovingMeshLandmarksFilename = optarg;
        std::cout << "Moving mesh landmarks:          " << inputMovingMeshLandmarksFilename << std::endl;
        break;
      case 't':
        inputTexturePrefix = optarg;
        std::cout << "  Input texture file prefix:  " << inputTexturePrefix << std::endl;
        break;
      case 'e':
        inputTexturePostfix = optarg;
        std::cout << "  Input texture file postfix: " << inputTexturePostfix << std::endl;
        break;
      case 'n':
        nTextureFiles = atoi( optarg );
        std::cout << "  Number of texture files:    " << nTextureFiles << std::endl;
        break;
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
  
  if( inputMeshFilenames.empty() )
  {
    std::cerr << "ERROR: No input meshes specified!" << std::endl;
    return EXIT_FAILURE;
  }
  if( inputMeshFilenames.size() != 2 )
  {
    std::cerr << "ERROR. Exactly two input filenames required (fixed and moving mesh)!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input meshes:
  
  vtkSmartPointer<vtkPolyData> fixedMesh;
  vtkSmartPointer<vtkPolyData> movingMesh;
  vtkSmartPointer<vtkPolyData> tempMesh;
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading fixed mesh ... " << std::flush;
  
  vtkSmartPointer<vtkPolyDataReader> fixedMeshReader;
  vtkTexturingHelper fixedMeshHelper;
  
  std::string filename = inputMeshFilenames[0];
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "Reading VTK data ..." << std::endl;
    fixedMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    fixedMeshReader->SetFileName( filename.c_str() );
    fixedMeshReader->Update();
    fixedMesh = fixedMeshReader->GetOutput();
  }
  else
  {
    std::cout << "Reading OBJ data ..." << std::endl;
    fixedMeshHelper.ReadGeometryFile( filename.c_str() );
    fixedMesh = fixedMeshHelper.GetPolyData();
  }

  std::cout << "OK." << std::endl;
  std::cout << "Loading moving mesh ... " << std::flush;
  
  vtkSmartPointer<vtkPolyDataReader> movingMeshReader;
  vtkTexturingHelper movingMeshHelper;
  
  filename = inputMeshFilenames[1];
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "Reading VTK data ..." << std::endl;
    movingMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    movingMeshReader->SetFileName( filename.c_str() );
    movingMeshReader->Update();
    movingMesh = movingMeshReader->GetOutput();
  }
  else
  {
    std::cout << "Reading OBJ data ..." << std::endl;
    movingMeshHelper.ReadGeometryFile( filename.c_str() );
    movingMesh = movingMeshHelper.GetPolyData();
  }
  
  std::cout << "OK." << std::endl;
  
  // TEMP
  vtkSmartPointer<vtkMatrix4x4> matrix1 = Load4x4Matrix( inputTransformFilenames[0] );
  
  
  
  std::cout << "Loading TEMP mesh ... " << std::flush;
  
  vtkSmartPointer<vtkPolyDataReader> tempMeshReader;
  vtkTexturingHelper tempMeshHelper;
  
  filename = inputTEMPMeshFilename;
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "Reading VTK data ..." << std::endl;
    tempMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    tempMeshReader->SetFileName( inputTEMPMeshFilename );
    tempMeshReader->Update();
    tempMesh = tempMeshReader->GetOutput();
  }
  else
  {
    std::cout << "Reading OBJ data ..." << std::endl;
    tempMeshHelper.ReadGeometryFile( inputTEMPMeshFilename );
    tempMesh = tempMeshHelper.GetPolyData();
  }
  
  std::cout << "OK." << std::endl;
  
  // Loading input landmarks (optional; if not specified, pseudo landmarks are generated):
  // TEMP SOLUTION:
  
  // FIXED LANDMARKS: ECO 06 (i.e. top part):
  
  vtkSmartPointer<vtkPoints> fixedMeshLandmarks = vtkSmartPointer<vtkPoints>::New();
  fixedMeshLandmarks->InsertNextPoint( 58.897, -41.0685, 7.2383 );
  fixedMeshLandmarks->InsertNextPoint( 58.7651, -20.1566, 15.8083 );
  fixedMeshLandmarks->InsertNextPoint( 44.8672, -10.762, 30.1874 );
  fixedMeshLandmarks->InsertNextPoint( -17.8375, 32.5014, 22.5529 );
  fixedMeshLandmarks->InsertNextPoint( -37.6899, 38.714, 13.7837 );
  
  // TEMP LANDMARKS: ECO 06 (i.e. top part):
  //YEAH!
  vtkSmartPointer<vtkPoints> tempFixedMeshLandmarks = vtkSmartPointer<vtkPoints>::New();
  tempFixedMeshLandmarks->InsertNextPoint( -35.1626, 46.8751, -18.2815 );
  tempFixedMeshLandmarks->InsertNextPoint( -13.6249, 53.9925, -20.7532 );
  tempFixedMeshLandmarks->InsertNextPoint( 10.5811, 44.9248, -16.1256 );
  tempFixedMeshLandmarks->InsertNextPoint( 33.0445, 25.6715, 3.91884 );
  tempFixedMeshLandmarks->InsertNextPoint( 30.9051, 13.487, 22.6511 );
  tempFixedMeshLandmarks->InsertNextPoint( 36.4909, -11.4603, 28.7534 );
  tempFixedMeshLandmarks->InsertNextPoint( 20.5124, 22.6622, -30.0448 );
  
  
  /*
  fixedMeshLandmarks->InsertNextPoint( -25.5974, -10.7835, -12.6974 );
  fixedMeshLandmarks->InsertNextPoint( 8.53846, 10.9905, -24.5173 );
  fixedMeshLandmarks->InsertNextPoint( 8.84473, 11.3597, -24.8401 );
  fixedMeshLandmarks->InsertNextPoint( -23.396, 6.57009, -26.3628 );
  fixedMeshLandmarks->InsertNextPoint( -22.2657, 6.98938, -26.3504 );
  /*
  fixedMeshLandmarks->InsertNextPoint( 9.21098, 11.3476, -24.7576 );
  fixedMeshLandmarks->InsertNextPoint( 8.07267, 12.3653, -33.7108 );
  fixedMeshLandmarks->InsertNextPoint( 2.26376, 6.36441, -27.0665 );
  fixedMeshLandmarks->InsertNextPoint( 16.0132, 9.60645, -27.4497 );
  fixedMeshLandmarks->InsertNextPoint( 10.9212, 3.7784, -21.174 );
  fixedMeshLandmarks->InsertNextPoint( -22.1363, 6.89309, -26.223 );
  fixedMeshLandmarks->InsertNextPoint( -29.9017, 2.19442, -28.8828 );
  fixedMeshLandmarks->InsertNextPoint( -23.3821, 8.50089, -35.1005 );
  fixedMeshLandmarks->InsertNextPoint( -14.1338, 5.39627, -28.4619 );
  fixedMeshLandmarks->InsertNextPoint( -19.5956, -1.10621, -22.6882 );
   */
  
  // MOVING LANDMARKS: ECO 06 ... NO; ENTIRE!
  
  
  vtkSmartPointer<vtkPoints> movingMeshLandmarks = vtkSmartPointer<vtkPoints>::New();
  // NO!
  movingMeshLandmarks->InsertNextPoint( -61.5773, -5.41659, 61.7419 );
  movingMeshLandmarks->InsertNextPoint( -46.0431, -21.476, 63.197 );
  movingMeshLandmarks->InsertNextPoint( -24.8988, -20.083, 69.0768 );
  movingMeshLandmarks->InsertNextPoint( 40.5567, -15.5864, 27.0243 );
  movingMeshLandmarks->InsertNextPoint( 52.2918, -9.44844, 9.61071 );
  
  // YEAH!
  vtkSmartPointer<vtkPoints> tempMovingMeshLandmarks = vtkSmartPointer<vtkPoints>::New();
  tempMovingMeshLandmarks->InsertNextPoint( 47.9933, 58.3786, -17.7889 );
  tempMovingMeshLandmarks->InsertNextPoint( 64.5093, 49.7574, -4.00626 );
  tempMovingMeshLandmarks->InsertNextPoint( 67.2316, 33.5068, 14.8721 );
  tempMovingMeshLandmarks->InsertNextPoint( 54.5685, 18.6301, 45.0772 );
  tempMovingMeshLandmarks->InsertNextPoint( 35.4663, 22.3284, 57.8411 );
  tempMovingMeshLandmarks->InsertNextPoint( 17.1097, 5.31267, 66.4572 );
  tempMovingMeshLandmarks->InsertNextPoint( 61.0975, 5.49275, 13.1036 );
  
  /*
  movingMeshLandmarks->InsertNextPoint( -21.941, -4.86134, -21.2955 );
  movingMeshLandmarks->InsertNextPoint( -41.2271, 21.4081, 7.39269 );
  movingMeshLandmarks->InsertNextPoint( -41.8412, 22.2552, 7.87222 );
  movingMeshLandmarks->InsertNextPoint( -46.507, 8.63798, -21.3082 );
  movingMeshLandmarks->InsertNextPoint( -45.8503, 8.64068, -20.4694 );
  /*
  movingMeshLandmarks->InsertNextPoint( -41.856, 21.983, 7.37117 );
  movingMeshLandmarks->InsertNextPoint( -44.2849, 29.1029, 4.47734 );
  movingMeshLandmarks->InsertNextPoint( -38.7051, 21.1185, 0.459296 );
  movingMeshLandmarks->InsertNextPoint( -38.7655, 28.5365, 12.8549 );
  movingMeshLandmarks->InsertNextPoint( -33.3208, 20.7873, 8.63892 );
  movingMeshLandmarks->InsertNextPoint( -46.2924, 8.8483, -20.6806 );
  movingMeshLandmarks->InsertNextPoint( -43.3967, 7.3167, -28.232 );
  movingMeshLandmarks->InsertNextPoint( -49.2672, 15.7597, -24.2769 );
  movingMeshLandmarks->InsertNextPoint( -42.6067, 15.7738, -14.7766 );
  movingMeshLandmarks->InsertNextPoint( -36.7485, 7.60761, -18.8141 );
   */
  
  // -------------------------------------------------------------
  // Computing transformation:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Computing transformation matrix ... " << std::flush;
  
  vtkSmartPointer<vtkLandmarkTransform> landmarkBasedTransform = vtkSmartPointer<vtkLandmarkTransform>::New();
  landmarkBasedTransform->SetSourceLandmarks( fixedMeshLandmarks );
  landmarkBasedTransform->SetTargetLandmarks( movingMeshLandmarks );
  landmarkBasedTransform->SetModeToRigidBody();
  landmarkBasedTransform->Update();
  
  // Display computed transformation matrix:
  
  vtkMatrix4x4* mat = landmarkBasedTransform->GetMatrix();
  std::cout << "Matrix: " << *mat;
  
  if( !outputTransformFilenames.empty() ) Write4x4Matrix( mat, outputTransformFilenames[0] );
  
  
  // -------------------------------------------------------------
  // Applying transformation:
  
  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  //transformFilter->SetInputData( movingMesh );
  transformFilter->SetInputData( fixedMesh );
  transformFilter->SetTransform( landmarkBasedTransform );
  transformFilter->Update();
  
  // -------------------------------------------------------------
  // Computing transformation:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Computing transformation matrix (2) ... " << std::flush;
  
  vtkSmartPointer<vtkLandmarkTransform> landmarkBasedTransformTEMP = vtkSmartPointer<vtkLandmarkTransform>::New();
  landmarkBasedTransformTEMP->SetSourceLandmarks( tempFixedMeshLandmarks );
  landmarkBasedTransformTEMP->SetTargetLandmarks( tempMovingMeshLandmarks );
  landmarkBasedTransformTEMP->SetModeToRigidBody();
  landmarkBasedTransformTEMP->Update();
  
  // Display computed transformation matrix:
  
  vtkMatrix4x4* matTEMP = landmarkBasedTransformTEMP->GetMatrix();
  std::cout << "Matrix: " << *matTEMP;
  
  if( !outputTransformFilenames.empty() ) Write4x4Matrix( mat, outputTransformFilenames[1] );
  
  // -------------------------------------------------------------
  // Applying transformation:
  
  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilterTEMP = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  //transformFilter->SetInputData( movingMesh );
  transformFilterTEMP->SetInputData( tempMesh );
  transformFilterTEMP->SetTransform( landmarkBasedTransformTEMP );
  transformFilterTEMP->Update();
  
  // -------------------------------------------------------------
  // Rendering data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Rendering data ... " << std::flush;
  
  // Create mapper and actor for fixed mesh:
  
  vtkSmartPointer<vtkPolyDataMapper> fixedMeshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  //fixedMeshMapper->SetInputConnection( fixedMeshReader->GetOutputPort() );
  fixedMeshMapper->SetInputConnection( transformFilter->GetOutputPort() );
  //fixedMeshMapper->ScalarVisibilityOff();
  
  vtkSmartPointer<vtkActor> fixedMeshActor = vtkSmartPointer<vtkActor>::New();
  fixedMeshActor->SetMapper( fixedMeshMapper );
  
  // Create mapper and actor for moving mesh:
  
  vtkSmartPointer<vtkPolyDataMapper> movingMeshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  //movingMeshMapper->SetInputConnection( transformFilter->GetOutputPort() );
  movingMeshMapper->SetInputData( movingMesh );
  movingMeshMapper->ScalarVisibilityOff();
  
  vtkSmartPointer<vtkActor> movingMeshActor = vtkSmartPointer<vtkActor>::New();
  movingMeshActor->SetMapper( movingMeshMapper );
  
  // Create mapper and actor for TEMP mesh:
  
  vtkSmartPointer<vtkPolyDataMapper> movingMeshMapperTEMP = vtkSmartPointer<vtkPolyDataMapper>::New();
  //movingMeshMapper->SetInputConnection( transformFilter->GetOutputPort() );
  movingMeshMapperTEMP->SetInputConnection( transformFilterTEMP->GetOutputPort() );
  movingMeshMapperTEMP->ScalarVisibilityOff();
  
  vtkSmartPointer<vtkActor> movingMeshActorTEMP = vtkSmartPointer<vtkActor>::New();
  movingMeshActorTEMP->SetMapper( movingMeshMapperTEMP );
  movingMeshActorTEMP->GetProperty()->SetColor(0,0,1);
 
  // Create a renderer, render window, and interactor:
  
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer( renderer );
 
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );
  
  vtkSmartPointer<vtkInteractorStyleSwitch> style = vtkSmartPointer<vtkInteractorStyleSwitch>::New();
  renderWindowInteractor->SetInteractorStyle( style );
  //vtkSmartPointer<MouseInteractorStylePP> style = vtkSmartPointer<MouseInteractorStylePP>::New();
  //renderWindowInteractor->SetInteractorStyle( style );
  //style->SetCurrentRenderer( renderer );
  
  vtkSmartPointer<vtkCallbackCommand> clickCallback;
  if( !outputMeshFilenames.empty() )
  {
    clickCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    clickCallback->SetCallback( ClickCallbackFunction );
    clickCallback->SetClientData( &outputMeshFilenames );
    renderWindowInteractor->AddObserver( vtkCommand::KeyPressEvent, clickCallback );
  }
 
  // Add the actor to the scene:
  
  renderer->AddActor( fixedMeshActor );
  //renderer->AddActor( movingMeshActor );
  renderer->AddActor( movingMeshActorTEMP );
  renderer->SetBackground( 1.0, 1.0, 1.0 );
 
  // Render and interact:
  
  renderWindow->SetWindowName( "icnsSurfaceRegistration" );
  renderWindow->Render();
  renderWindowInteractor->Start();
  
  // After closing:
  
  std::cout << "TEST" << fixedMeshActor->GetUserTransform() << std::endl;
 
  return EXIT_SUCCESS;
}


// ---------------------------------------------------------------
// Additional routines definitions:
// ---------------------------------------------------------------

void ClickCallbackFunction (vtkObject *caller,
                            long unsigned int eventID,
                            void * clientData,
                            void * callData )
{
  // Get the calling interactor:
  vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
  std::string key = iren->GetKeySym();
  
  // If return and the output filename is set, save transformed data to file:
  
  if(key == "Return")
  {
    std::cout << "Click callback: RETURN." << std::endl;
    std::cout << "Saving transformed polydata to file: " << std::endl;
    
    std::vector< std::string >* names = static_cast<std::vector< std::string > * >( clientData );
    std::cout << "SIZE OF POINTER: " << names->size() << std::endl;
    std::cout << (*names)[0] << std::endl;
    std::cout << (*names)[1] << std::endl;
    
    // Extract renderer from caller data:
    
    vtkSmartPointer<vtkRenderer> currentRenderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    
    for( unsigned int iActorID = 0; iActorID < 2; iActorID++ )
    {
      vtkSmartPointer<vtkActor> currentActor       = static_cast<vtkActor*>( currentRenderer->GetActors()->GetItemAsObject( iActorID ) );
      vtkMatrix4x4* currentTransformMatrix         = currentActor->GetMatrix();
      vtkSmartPointer<vtkPolyData> currentPolyData = static_cast<vtkPolyData*>( currentActor->GetMapper()->GetInput() );
      
      vtkSmartPointer<vtkPolyDataWriter> currentWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
      currentWriter->SetFileName( (*names)[iActorID].c_str() );
      currentWriter->SetInputData( currentPolyData );
      currentWriter->Write();
    }
  }
} // end of ClickCallbackFunction()

void Write4x4Matrix( vtkMatrix4x4* matrix, std::string filename )
{
  ofstream txtFile;
  txtFile.open( filename );
  
  for( unsigned int i = 0; i < 4; i++ )
  {
    for( unsigned int j = 0; j < 4; j++ )
    {
      txtFile << matrix->GetElement( i, j );
      if( j < 3 )  txtFile << " ";
      if( j == 3 ) txtFile << "\n";
    }
  }
  
  txtFile.close();
} // end of Write4x4Matrix()

vtkMatrix4x4* Load4x4Matrix( std::string filename )
{
  std::cout << "Reading matrix from file: " << filename << "..." << std::endl;
  
  vtkMatrix4x4* mat = vtkMatrix4x4::New();
  
  std::string line;
  ifstream txtFile( filename );
  
  if( txtFile.is_open() )
  {
    for( unsigned int i = 0; i < 4; i++ )
    {
      for( unsigned int j = 0; j < 4; j++ )
      {
        txtFile >> mat->Element[i][j];
      }
    }
    txtFile.close();
  }
  
  std::cout << *mat << std::endl;

  return mat;
}
