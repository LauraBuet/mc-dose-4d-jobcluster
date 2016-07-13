/** \file icnsViewOBJMultiTexture.cpp
 *
 *  Nomen est omen: shows an OBJ file with multiple texture files.
 *  Note that the standard VTK pipeline is not able to do so, i.e. to 
 *  sufficiently read the mtl-file. 
 *
 *  As a workaround, here the helper classes from 
 *  https://github.com/Scylardor/vtkTexturingHelper
 *  are used. This means that the mtl is not read, but the individual 
 *  texture files are read.
 *
 *  For more information see
 *  http://scylardor.fr/2013/05/06/making-multi-texturing-work-with-vtk/
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
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>

// Project includes:
#include "vtkTexturingHelper.h"

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage :\n";
  std::cout << "icnsViewOBJMultiTexture -I <OBJ file> -T <texture file prefix> -E <texture file ending> -N <number of texture files> \n\n";
    
  std::cout << "-I <OBJ file>                 OBJ file to view.\n";
  std::cout << "-T <texture file prefix>      Texture file prefix.\n";
  std::cout << "-E <texture file ending>      Texture file ending (e.g. .png).\n";
  std::cout << "-N <number of texture files>  Texture file ending.\n";
  
  std::cout << "-h                            Print this help.\n";
  std::cout << "\n";
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
  
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsViewOBJ (multiple texture files)" << std::endl;
  std::cout << "==========================================" << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  char* inputOBJFile        = NULL;
  char* inputTexturePrefix  = NULL;
  char* inputTexturePostfix = NULL;
  unsigned int nTextureFiles = 0;
  
  // Reading parameters:
  
  while( (c = getopt( argc, argv, "I:T:E:N:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputOBJFile = optarg;
        std::cout << "  Input OBJ file:             " << inputOBJFile << std::endl;
        break;
      case 'T':
        inputTexturePrefix = optarg;
        std::cout << "  Input texture file prefix:  " << inputTexturePrefix << std::endl;
        break;
      case 'E':
        inputTexturePostfix = optarg;
        std::cout << "  Input texture file postfix: " << inputTexturePostfix << std::endl;
        break;
      case 'N':
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
  
  // -------------------------------------------------------------
  // Check input:
  
  if( inputOBJFile == NULL )
  {
    std::cerr << "ERROR: No OBJ filename!" << std::endl;
    return EXIT_FAILURE;
  }

  // -------------------------------------------------------------
  // Setting up helper class:
  
  vtkTexturingHelper helper;
  try
  {
    // Read the geometry file and extract texture coordinates data from OBJ file:
    helper.ReadGeometryFile( inputOBJFile );
    
    // Read texture files:
    if( nTextureFiles > 0 )
    {
      helper.ReadTextureFiles( inputTexturePrefix, inputTexturePostfix, nTextureFiles );
      helper.ApplyTextures();
    }
  }
  catch (const vtkTexturingHelperException & helperExc)
  {
    std::cout << "Whoops! ERRROR during texture import by helper class! " << helperExc.what() << std::endl;
  }
  
  // -------------------------------------------------------------
  // Generating visualization pipeline:
  
  vtkSmartPointer<vtkActor> actor = helper.GetActor();
  actor->GetProperty()->SetRepresentationToWireframe();
  
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer( renderer );
  
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );
  //vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  //renderWindowInteractor->SetRenderWindow( renderWindow );
  
  // Create appropriate interactor/interaction style:
  
  //vtkSmartPointer<vtkInteractorStyleSwitch> style = vtkSmartPointer<vtkInteractorStyleSwitch>::New();
  //iren->SetInteractorStyle( style );
  
  // Add the actors to the renderer, set the background and size:
  
  renderer->AddActor( actor );
  renderer->SetBackground( 1.0, 1.0, 1.0 );
  renderer->ResetCamera();
  renderWindow->SetSize( 800, 600 );
  renderWindow->Render();

  // -------------------------------------------------------------
  // Visualizing data (rendering):

  // This starts the event loop and as a side effect causes an initial render:
  renderWindowInteractor->Start();
    
  // Cleaning up:
  renderer->Delete();
  renderWindow->Delete();
  renderWindowInteractor->Delete();
  
  return EXIT_SUCCESS;
}
