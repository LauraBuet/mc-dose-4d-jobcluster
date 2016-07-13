/** \file icnsVE3D_defineLandmarks.cxx
 *
 *  Project: VonEconomo3D
 *  Purpose: Write a list of landmarks that can be used for subsequent 
 *  landmarmark-based registration.
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2016 Department of Computational Neuroscience,
 *     University Medical Center Hamburg-Eppendorf
 *
 ****************************************************************************/

// System includes:
#include <iostream>
#include <fstream>
#include <string>
extern "C"
{
#include "getopt.h"
}

// ITK includes: None so far.

// VTK includes:
#include <vtkActor.h>
#include <vtkObjectFactory.h>
#include <vtkLandmarkTransform.h>
#include <vtkLookupTable.h>
#include <vtkMatrix4x4.h>
#include <vtkOBJReader.h>
#include <vtkPoints.h>
#include <vtkPointPicker.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransformPolyDataFilter.h>

#include <vtkVector.h>

#include <vtkInteractorStyleTrackballCamera.h>

// Project includes:
#include "vtkTexturingHelper.h"

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage :\n";
  std::cout << "icnsVonEconomo -I <input vtk mesh> -O <output landmarks>\n\n";
    
  std::cout << "-I <input vtk mesh>            Filename of the mesh to be displayed.\n";
  std::cout << "-O <output landmarks>          Filename of txt-file to write selected landmarks at..\n";
  
  std::cout << "-h                             Print this help.\n";
  std::cout << "\n";
}

// ---------------------------------------------------------------
// Helper class for desired interaction between scene and user:
// ---------------------------------------------------------------

class MouseInteractorStylePP : public vtkInteractorStyleTrackballCamera
{
public:
  static MouseInteractorStylePP* New();
  vtkTypeMacro(MouseInteractorStylePP, vtkInteractorStyleTrackballCamera);
    
  char* landmarkListFilename = NULL;
  std::vector<vtkVector3f> pickedPoints;
  vtkVector3f currentlyPickedPoint;
  
  virtual void OnKeyPress()
  {
    // Get the keypress
    vtkRenderWindowInteractor *rwi = this->Interactor;
    std::string key = rwi->GetKeySym();
    
    // Output the key that was pressed
    std::cout << "Pressed " << key << std::endl;
    
    // Handle an arrow key
    if(key == "Up")
    {
      std::cout << "The up arrow was pressed." << std::endl;
    }
    
    // Handle a "normal" key
    if(key == "a")
    {
      std::cout << "The a key was pressed." << std::endl;
    }
   
    vtkInteractorStyleTrackballCamera::OnKeyPress();
  }
  
  //void OnLeftButtonDown()
  void OnRightButtonDown()
  {
    // Just getting and outputting position of picked point:
      
    std::cout << "Picking pixel: " << this->Interactor->GetEventPosition()[0] << " " << this->Interactor->GetEventPosition()[1] << std::endl;
    
    this->Interactor->GetPicker()->Pick( this->Interactor->GetEventPosition()[0],
                                         this->Interactor->GetEventPosition()[1],
                                         0,  // always zero.
                                         this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer() );
    double pickedPosition[3];
    this->Interactor->GetPicker()->GetPickPosition( pickedPosition );
    
    //std::cout << "Picked value: " << pickedPosition[0] << " " << pickedPosition[1] << " " << pickedPosition[2] << std::endl;
    currentlyPickedPoint.Set( pickedPosition[0], pickedPosition[1], pickedPosition[2] );
    pickedPoints.push_back( currentlyPickedPoint );
    
    ofstream myfile( landmarkListFilename );
    std::cout << "List of picked points: " << std::endl;
    for( unsigned int i=0; i < pickedPoints.size(); i++ )
    {
      std::cout << "Point " << i << ": [" << pickedPoints[i].GetX() << ";" << pickedPoints[i].GetY() << ";" << pickedPoints[i].GetZ() << "]" << std::endl;
      if( ( landmarkListFilename != NULL ) && ( myfile.is_open() ) )
      {
        myfile << pickedPoints[i].GetX() << " " << pickedPoints[i].GetY() << " " << pickedPoints[i].GetZ() << "\n";
      }
    }
    myfile.close();
    
    // Drawing spheres at picked positions:
    
    vtkSmartPointer<vtkSphereSource> newSphereSource = vtkSmartPointer<vtkSphereSource>::New();
    newSphereSource->SetRadius( 2 );
    newSphereSource->SetCenter( pickedPosition );
    
    vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    newMapper->SetInputConnection( newSphereSource->GetOutputPort() );
    
    vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
    newActor->SetMapper( newMapper );
    //newActor->GetProperty()->SetColor(1,0,0);
    
    this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddViewProp( newActor );
    this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->Render();
    
    // Forward events:
      
    //vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    vtkInteractorStyleTrackballCamera::OnRightButtonDown();
  }
};

vtkStandardNewMacro(MouseInteractorStylePP);


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
  std::cout << "icnsVE3D: define landmarks" << std::endl;
  std::cout << "==========================================" << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  char* inputMeshFilename = NULL;
  char* outputLandmarkListFilename = NULL;
  
  // Reading parameters:
  
  while( (c = getopt( argc, argv, "I:O:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputMeshFilename = optarg;
        std::cout << "Input mesh:            " << inputMeshFilename << std::endl;
        break;
      case 'O':
        outputLandmarkListFilename = optarg;
        std::cout << "Output landmark list:  " << outputLandmarkListFilename << std::endl;
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
  
  if( inputMeshFilename == NULL )
  {
    std::cerr << "ERROR: No input mesh filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // Reading data:
  
  vtkSmartPointer<vtkPolyData> inputMesh;
  vtkSmartPointer<vtkPolyDataReader> inputMeshReader;
  vtkTexturingHelper helper;
  
  // Extracting file ending:
  
  std::string filename = inputMeshFilename;
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "Reading VTK data ..." << std::endl;
    inputMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    inputMeshReader->SetFileName( inputMeshFilename );
    inputMeshReader->Update();
    inputMesh = inputMeshReader->GetOutput();
    
  }
  else
  {
    std::cout << "Reading OBJ data ..." << std::endl;
    helper.ReadGeometryFile( inputMeshFilename );
    inputMesh = helper.GetPolyData();
  }

  // Assembling rendering pipeline:
  
  vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputData( inputMesh );
  mapper->ScalarVisibilityOff();
  
  vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper( mapper );
  
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  renderer->AddActor(actor);
  
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  renderer->SetBackground(1,1,1); // Background color white
  
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow) ;

  // Create appropriate interactor/interaction style:
  vtkSmartPointer<MouseInteractorStylePP> style = vtkSmartPointer<MouseInteractorStylePP>::New();
  if( outputLandmarkListFilename != NULL ) style->landmarkListFilename = outputLandmarkListFilename;
  //vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle( style );
  style->SetCurrentRenderer( renderer );
  
  
  // Starting visualization:
  
  renderWindow->Render();
  renderWindowInteractor->Start();
 
  return EXIT_SUCCESS;
}
