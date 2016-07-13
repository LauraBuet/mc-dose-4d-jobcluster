/** \file icnsLandmarkBasedSurfaceRegistration.cpp
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

#include <vtkMatrixToLinearTransform.h>

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
  
  // -------------------------------------------------------------
  // Loading initial transforms and transforming meshes
  // (if filenames are specified):
  
  if( !inputTransformFilenames.empty() )
  {
    // Fixed mesh:
    
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Loading input transforms ... " << std::flush;
    
    vtkSmartPointer<vtkMatrix4x4> fixedMeshTransformMatrix = Load4x4Matrix( inputTransformFilenames[0] );
    
    vtkSmartPointer<vtkMatrixToLinearTransform> fixedMeshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
    fixedMeshTransform->SetInput( fixedMeshTransformMatrix );
    
    vtkSmartPointer<vtkTransformPolyDataFilter> fixedMeshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    fixedMeshTransformFilter->SetInputData( fixedMesh );
    fixedMeshTransformFilter->SetTransform( fixedMeshTransform );
    fixedMeshTransformFilter->Update();
    
    fixedMesh = fixedMeshTransformFilter->GetOutput();
    
    // Moving mesh:
    
    vtkSmartPointer<vtkMatrix4x4> movingMeshTransformMatrix = Load4x4Matrix( inputTransformFilenames[1] );
    
    vtkSmartPointer<vtkMatrixToLinearTransform> movingMeshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
    movingMeshTransform->SetInput( movingMeshTransformMatrix );
    
    vtkSmartPointer<vtkTransformPolyDataFilter> movingMeshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    movingMeshTransformFilter->SetInputData( movingMesh );
    movingMeshTransformFilter->SetTransform( movingMeshTransform );
    movingMeshTransformFilter->Update();
    
    movingMesh = movingMeshTransformFilter->GetOutput();
  }
  
 
 
  
  // -------------------------------------------------------------
  // Loading landmarks and computing transformation: NYI
  
  /*
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
  */
  
  // -------------------------------------------------------------
  // Rendering data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Rendering data ... " << std::flush;
  
  // Create mapper and actor for fixed mesh:
  
  vtkSmartPointer<vtkPolyDataMapper> fixedMeshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  fixedMeshMapper->SetInputData( fixedMesh );
  //fixedMeshMapper->SetInputConnection( transformFilter->GetOutputPort() );
  //fixedMeshMapper->ScalarVisibilityOff();
  
  vtkSmartPointer<vtkActor> fixedMeshActor = vtkSmartPointer<vtkActor>::New();
  fixedMeshActor->SetMapper( fixedMeshMapper );
  
  // Create mapper and actor for moving mesh:
  
  vtkSmartPointer<vtkPolyDataMapper> movingMeshMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  movingMeshMapper->SetInputData( movingMesh );
  movingMeshMapper->ScalarVisibilityOff();
  
  vtkSmartPointer<vtkActor> movingMeshActor = vtkSmartPointer<vtkActor>::New();
  movingMeshActor->SetMapper( movingMeshMapper );
  //movingMeshActorTEMP->GetProperty()->SetColor(0,0,1);
 
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
  
  // Generating callback and callback clientdata:
  
  vtkSmartPointer<vtkCallbackCommand> clickCallback;
  std::vector< std::string > callbackClientData = outputMeshFilenames;
  if( !callbackClientData.empty() )
  {
    // If output transform filenames are specified, also pass
    // them to callback:
    
    if( !outputTransformFilenames.empty() )
    {
      callbackClientData.push_back( outputTransformFilenames[0] );
      callbackClientData.push_back( outputTransformFilenames[1] );
    }
    
    // Finally generate callback and observer:
    
    clickCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    clickCallback->SetCallback( ClickCallbackFunction );
    clickCallback->SetClientData( &callbackClientData );
    renderWindowInteractor->AddObserver( vtkCommand::KeyPressEvent, clickCallback );
  }
 
  // Add the actor to the scene:
  
  renderer->AddActor( fixedMeshActor );
  renderer->AddActor( movingMeshActor );
  renderer->SetBackground( 1.0, 1.0, 1.0 );
 
  // Render and interact:
  
  renderWindow->SetWindowName( "icnsSurfaceRegistration" );
  renderWindow->Render();
  renderWindowInteractor->Start();
  
  // Cleaning up: TODO
  
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
      
      if( (names->size() == 4) && (iActorID == 0) )
      {
        Write4x4Matrix( currentTransformMatrix, (*names)[2] );
      }
      if( (names->size() == 4) && (iActorID == 1) )
      {
        Write4x4Matrix( currentTransformMatrix, (*names)[3] );
      }
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
  else
  {
    std::cerr << "ERROR: cannot open filename. Returning identity matrix." << std::endl;
  }
  
  std::cout << *mat << std::endl;
  return mat;
}
