/** \file icnsManualSurfaceRegistration.cpp
 *
 *  Purpose:
 *  Load two or more surface models, manually adjust their position and
 *  save resulting models (and/or transformations).
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
#include <vtkMatrix4x4.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkPolyDataCollection.h>
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
#include <vtkProperty.h>
#include <vtkRendererCollection.h>

#include <vtkMatrixToLinearTransform.h>

// Project includes:
#include "vtkTexturingHelper.h" // required to load multitexture OBJs

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsManualSurfaceRegistration -F <fixed vtk mesh> -M <moving vtk mesh> [...]\n" << std::endl;
  std::cout << std::endl;
  std::cout << "INPUT DATA" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "-I <list of surface meshes>         IN: Filename list of meshes to be displayed and transformed." << std::endl;
  std::cout << "                                        Possible endings: VTK or OBJ." << std::endl;
  std::cout << "-T <list of transformations>        IN: Filename list of initial transformations (4x4 matrix)" << std::endl;
  std::cout << "                                        to be applied to meshes." << std::endl;
  std::cout << "-t <list of texture file prefixes>  IN: List of texture prefixes (e.g. ..._eco06_hell)." << std::endl;
  std::cout << "-e <list of texture file prefixes>  IN: List of texture postfixes (e.g. .png)."  << std::endl;
  std::cout << "-n <number of texture files>        IN: Number of texture files."  << std::endl;
  std::cout << std::endl;
  std::cout << "OUTPUT DATA" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "-O <list of surface meshes>         OUT: Filename list of output meshes (only VTK)." << std::endl;
  std::cout << "-W <list of transformations>        OUT: Filename lust of transformations."  << std::endl;
  std::cout << "                                         NB: If initial transformations are given, the output" << std::endl;
  std::cout << "                                         will be a concat of the initial and the actual tranform." << std::endl ;
  std::cout << std::endl;
  std::cout << "-h                                  Print this help." << std::endl;
  std::cout << std::endl;
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
  if( argc < 2 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  std::cout << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "icnsManualSurfaceRegistration" << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  
  // -------------------------------------------------------------
  // Initializing parameters with default values:
  
  std::vector<std::string>   inputMeshFilenames;
  std::vector<std::string>   inputTextureFilePrefixes;
  std::vector<std::string>   inputTextureFilePostfixes;
  std::vector< unsigned int> inputNumberOfTextureFiles;
  std::vector<std::string>   inputTransformFilenames;
  
  std::vector< std::string > outputMeshFilenames;
  std::vector< std::string > outputTransformFilenames;
  
  bool readingOBJdata = false;
  
  // Initializing renderer:
  
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  
  // Reading parameters.
  // Note that the filnames of the lists are considered to *not* contain
  // a '-', because this is assumed to be reserved and to start the next
  // parameter readout.
  
  int c;
  int cnt = 0;
  
  while( (c = getopt( argc, argv, "I:T:t:e:n:O:W:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputMeshFilenames.push_back( argv[optind] );
          std::cout << "  IN-Mesh [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'T':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTransformFilenames.push_back( argv[optind] );
          std::cout << "  IN-Transform [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 't':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTextureFilePrefixes.push_back( argv[optind] );
          std::cout << "  Texture file prefix [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'e':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTextureFilePostfixes.push_back( argv[optind] );
          std::cout << "  Texture file postfix [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'n':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputNumberOfTextureFiles.push_back( atoi(argv[optind]) );
          std::cout << "  Number of texture files [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputMeshFilenames.push_back( argv[optind] );
          std::cout << "  OUT-Mesh [" << cnt << "] = " << argv[optind] << std::endl;
        }
        break;
      case 'W':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          outputTransformFilenames.push_back( argv[optind] );
          std::cout << "  OUT-Transform [" << cnt << "] = " << argv[optind] << std::endl;
        }
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
  if( inputMeshFilenames.size() < 2 )
  {
    std::cerr << "ERROR. Exactly two input filenames required (fixed and moving mesh)!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input meshes:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input meshes ... " << std::endl;
  
  vtkSmartPointer<vtkPolyDataCollection> polyDataCollection = vtkSmartPointer<vtkPolyDataCollection>::New();
  vtkSmartPointer<vtkPolyDataReader> polyDataReader;
  
  for( unsigned iFilenames = 0; iFilenames < inputMeshFilenames.size(); iFilenames++ )
  {
    std::string currentFilename = inputMeshFilenames[iFilenames];
    vtkSmartPointer<vtkPolyData> currentPolyData = vtkSmartPointer<vtkPolyData>::New();
    
    // Extracting file ending and choosing approproate reader:
    
    // First case: Standard VTK
    if( currentFilename.substr( currentFilename.find_last_of(".") + 1) == "vtk" )
    {
      std::cout << "  Reading VTK data: " << currentFilename << std::endl;
      
      polyDataReader = vtkSmartPointer<vtkPolyDataReader>::New();
      polyDataReader->SetFileName( currentFilename.c_str() );
      polyDataReader->Update();
      
      currentPolyData->DeepCopy( polyDataReader->GetOutput() );
      polyDataCollection->AddItem( currentPolyData );
      
      // Reading input matrix (if specified):
      if( ( !inputTransformFilenames.empty() ) &&
          ( std::strcmp( inputTransformFilenames[iFilenames].c_str(), "NONE" ) != 0 ) )
      {
        std::cout << "  Reading corresponding input transform: " << inputTransformFilenames[iFilenames] << std::endl;
        vtkSmartPointer<vtkMatrix4x4> currentTransformMatrix = Load4x4Matrix( inputTransformFilenames[iFilenames] );
        
        // Printing matrix:
        for( unsigned int i = 0; i < 4; i++ )
        {
          std::cout << "    ";
          for( unsigned int j = 0; j < 4; j++ )
          {
            std::cout << currentTransformMatrix->Element[i][j] << " ";
          }
          std::cout << std::endl;
        }
        
        vtkSmartPointer<vtkMatrixToLinearTransform> meshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
        meshTransform->SetInput( currentTransformMatrix );
        
        vtkSmartPointer<vtkTransformPolyDataFilter> meshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        meshTransformFilter->SetInputData( static_cast<vtkPolyData*>( polyDataCollection->GetItemAsObject( iFilenames ) ) );
        meshTransformFilter->SetTransform( meshTransform );
        meshTransformFilter->Update();
        
        polyDataCollection->ReplaceItem( iFilenames, meshTransformFilter->GetOutput() );
      }
    }
    // Second case: OBJ
    else
    {
      std::cout << "  Reading OBJ data ..." << currentFilename << std::endl;
      vtkTexturingHelper objReader;
      readingOBJdata = true;
      
      std::cout << "  Reading geometry file ... " << currentFilename << std::endl;
      objReader.ReadGeometryFile( currentFilename );
      
      // Read texture files:
      if( !inputNumberOfTextureFiles.empty() )
      {
        std::cout << "  Reading textures ... " << std::endl;
        objReader.ReadTextureFiles( inputTextureFilePrefixes[iFilenames],
                                   inputTextureFilePostfixes[iFilenames],
                                   inputNumberOfTextureFiles[iFilenames] );
        std::cout << "  Apply textures ... " << std::endl;
        objReader.ApplyTextures();
      }
      
      // Reading input matrix (if specified):
      if( ( !inputTransformFilenames.empty() ) &&
          ( std::strcmp( inputTransformFilenames[iFilenames].c_str(), "NONE" ) != 0 ) )
      {
        vtkSmartPointer<vtkMatrix4x4> currentTransformMatrix = Load4x4Matrix( inputTransformFilenames[iFilenames] );
        vtkSmartPointer<vtkMatrixToLinearTransform> meshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
        
        // Printing matrix:
        for( unsigned int i = 0; i < 4; i++ )
        {
          std::cout << "  ";
          for( unsigned int j = 0; j < 4; j++ )
          {
            std::cout << currentTransformMatrix->Element[i][j] << " ";
          }
          if( i < 3 ) std::cout << std::endl;
        }

        meshTransform->SetInput( currentTransformMatrix );
        objReader.GetActor()->SetUserTransform( meshTransform );
      }
      
      // Directly add actor to renderer:
      //currentPolyData->DeepCopy( objReader.GetPolyData() );
      renderer->AddActor( objReader.GetActor() );
    }
  }
  std::cout << "OK." << std::endl;
  

  // -------------------------------------------------------------
  // Setting up data pipeline:
  
  // Generating mapper and actor for each input file and
  // adding each actor to renderer:
  
  if( !readingOBJdata )
  {
    for( unsigned iFilenames = 0; iFilenames < inputMeshFilenames.size(); iFilenames++ )
    {
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputData( static_cast< vtkPolyData* >( polyDataCollection->GetItemAsObject( iFilenames ) ) );
      //mapper->SetInputConnection( inputMeshReader->GetOutputPort() );
      mapper->ScalarVisibilityOff();
      
      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper( mapper );
      
      renderer->AddActor( actor );
    }
  }
  
  // Generating renderWindow and adding renderer to window:
  
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetWindowName( "icnsManualRegistration" );
  renderWindow->AddRenderer( renderer );
  
  // An interactor:
  
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );
  
  // Create appropriate interactor/interaction style:
  //vtkSmartPointer<MouseInteractorStylePP> style =
  //  vtkSmartPointer<MouseInteractorStylePP>::New();
  //vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  //renderWindowInteractor->SetInteractorStyle(style);
  //style->SetCurrentRenderer(renderer);
  
  // Generating callback and callback clientdata:
  
  vtkSmartPointer<vtkCallbackCommand> clickCallback;
  std::vector< std::string > callbackClientData = outputMeshFilenames;
  if( !callbackClientData.empty() )
  {
    // If output transform filenames are specified, also pass
    // them to callback:
    
    if( !outputTransformFilenames.empty() )
    {
      callbackClientData.insert( callbackClientData.end(), outputTransformFilenames.begin(), outputTransformFilenames.end() );
      
      // And, if input transforms are given, also pass them:
      if( !inputTransformFilenames.empty() )
      {
        callbackClientData.insert( callbackClientData.end(), inputTransformFilenames.begin(), inputTransformFilenames.end() );
      }
    }
    
    // Finally generate callback and observer:
    
    clickCallback = vtkSmartPointer<vtkCallbackCommand>::New();
    clickCallback->SetCallback( ClickCallbackFunction );
    clickCallback->SetClientData( &callbackClientData );
    renderWindowInteractor->AddObserver( vtkCommand::KeyPressEvent, clickCallback );
  }
  
  //renderer->AddActor( actor );
  renderer->SetBackground( 1.0, 1.0, 1.0 ); // Background color white
  renderer->GradientBackgroundOn();
  renderer->SetBackground( 1, 1, 1 );
  renderer->SetBackground2( 0, 0, 0 );
  renderer->ResetCamera();
  
  renderWindow->SetSize( 800, 600 );
  renderWindow->Render();
  
  // -------------------------------------------------------------
  // Visualizing data (rendering):
  
  // This starts the event loop and as a side effect causes an initial render:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Starting visualization. Key bindings:"    << std::endl;
  std::cout << "  [A]      Changing to individual actor." << std::endl;
  std::cout << "  [M]      Changing back to mapper (no change of actor transformation)" << std::endl;
  std::cout << "  [RETURN] Saving output data."           << std::endl;
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}


// ---------------------------------------------------------------
// Additional routines definitions:
// ---------------------------------------------------------------

void ClickCallbackFunction( vtkObject *caller,
                            long unsigned int eventID,
                            void * clientData,
                            void * callData )
{
  // Get the calling interactor:
  
  vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
  std::string key = iren->GetKeySym();
  
  // If return and the output filename is set, save transformed data to file:
  
  if( key == "Return" )
  {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Saving current state/data to file ... " << std::endl;
    
    // Get current renderer from caller data to extract desired actors:
    
    vtkSmartPointer<vtkRenderer> currentRenderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    unsigned int nActors = currentRenderer->GetActors()->GetNumberOfItems();
    
    // Cast clientdata to desired output filenames.
    // If more filenames are present than actors, these are interpreted as
    // transform filenames. At this, the question is whether initial transforms
    // where specified. Thus, structurally, if the size of the list of filenames
    // is 3x the number of actors, it is assumed that
    //   1. the polydata names are given,
    //   2. the outputtranform filenames are given,
    //   3. the inputtransform filenames are given.
    // These are transferred to the corresponding vectors.
    
    std::vector< std::string >* entireListOfFilenames = static_cast<std::vector< std::string > * >( clientData );
    std::vector< std::string >  outputPolyDataFilenames;
    std::vector< std::string >  outputTransformFilenames;
    std::vector< std::string >  inputTransformFilenames;
    
    unsigned int nFilenames = entireListOfFilenames->size();
    
    // Plausibility check: The number of filenames should be either
    // identical or 2x or 3x the number of actors:
    
    if( !( ( nFilenames == nActors )   ||
           ( nFilenames == 2*nActors ) ||
           ( nFilenames == 3*nActors ) ) )
    {
      std::cerr << "ERROR: number of filenames in clientdata does not match number of actors in caller." << std::endl;
      return;
    }
    
    // If passed plausibility check, generate filename vectors:
    
    for( unsigned int iFilenames = 0; iFilenames < nFilenames; iFilenames++ )
    {
      std::cout << (*entireListOfFilenames)[iFilenames] << std::endl;
      
      // 1st case: polydata filenames (mandatory):
      if( iFilenames < nActors )
      {
        outputPolyDataFilenames.push_back( (*entireListOfFilenames)[iFilenames] );
        continue;
      }
      // 2nd case: output transform filenames (optional):
      if( ( nActors <= iFilenames ) && ( iFilenames < 2*nActors ) )
      {
        outputTransformFilenames.push_back( (*entireListOfFilenames)[iFilenames] );
        continue;
      }
      // 3rd case: input transform filenames (optional):
      if( ( 2*nActors <= iFilenames ) && ( iFilenames < 3*nActors ) )
      {
        inputTransformFilenames.push_back( (*entireListOfFilenames)[iFilenames] );
      }
    }
    
    // Writing data: Polydata
    
    std::cout << "Saving transformed polydata to file ... " << std::endl;
    for( unsigned int iFiles = 0; iFiles < outputPolyDataFilenames.size(); iFiles++ )
    {
      std::cout << "  " << outputPolyDataFilenames[iFiles] << std::endl;
      
      vtkSmartPointer<vtkActor> currentActor       = static_cast<vtkActor*>( currentRenderer->GetActors()->GetItemAsObject( iFiles ) );
      vtkSmartPointer<vtkPolyData> currentPolyData = static_cast<vtkPolyData*>( currentActor->GetMapper()->GetInput() );
      
      vtkSmartPointer<vtkPolyDataWriter> currentWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
      currentWriter->SetFileName( (*entireListOfFilenames)[iFiles].c_str() );
      currentWriter->SetInputData( currentPolyData );
      currentWriter->Write();
    }
    std::cout << "OK." << std::endl;
    
    // Writing data: Transformation files
    
    if( !outputTransformFilenames.empty() )
    {
      std::cout << "Saving output transformation to file ... " << std::endl;
      for( unsigned int iFiles = 0; iFiles < outputTransformFilenames.size(); iFiles++ )
      {
        std::cout << "  " << outputTransformFilenames[iFiles] << std::endl;
        
        vtkSmartPointer<vtkActor> currentActor       = static_cast<vtkActor*>( currentRenderer->GetActors()->GetItemAsObject( iFiles ) );
        vtkMatrix4x4* currentTransformMatrix         = currentActor->GetMatrix();
        
        // Printing matrix:
        for( unsigned int i = 0; i < 4; i++ )
        {
          std::cout << "    ";
          for( unsigned int j = 0; j < 4; j++ )
          {
            std::cout << currentTransformMatrix->Element[i][j] << " ";
          }
          std::cout << std::endl;
        }
        
        // If a corresponding input transform is specified, load it
        // premultiply it to current transform:
        
        if( ( !inputTransformFilenames.empty() ) &&
            ( std::strcmp( inputTransformFilenames[iFiles].c_str(), "NONE" ) != 0 ) )
        {
          std::cout << "  Preloading input transform " << inputTransformFilenames[iFiles] << " ... " << std::endl;
          
          vtkMatrix4x4* currentInputTransformMatrix = Load4x4Matrix( inputTransformFilenames[iFiles] );
          
          // Printing matrix:
          
          for( unsigned int i = 0; i < 4; i++ )
          {
            std::cout << "    ";
            for( unsigned int j = 0; j < 4; j++ )
            {
              std::cout << currentInputTransformMatrix->Element[i][j] << " ";
            }
            std::cout << std::endl;
          }
          
          vtkMatrix4x4::Multiply4x4( currentTransformMatrix, currentInputTransformMatrix, currentTransformMatrix );
          std::cout << "  Combined matrix:" << std::endl;
          
          for( unsigned int i = 0; i < 4; i++ )
          {
            std::cout << "    ";
            for( unsigned int j = 0; j < 4; j++ )
            {
              std::cout << currentTransformMatrix->Element[i][j] << " ";
            }
            std::cout << std::endl;
          }
        }
        
        // Finally writing the matrix:
        
        Write4x4Matrix( currentTransformMatrix, outputTransformFilenames[iFiles] );
      }
      std::cout << "OK." << std::endl;
    }
    
  } // end of callback "RETURN"
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
  return mat;
}
