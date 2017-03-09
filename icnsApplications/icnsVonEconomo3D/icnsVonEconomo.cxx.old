/** \file icnsVonEconomo.cpp
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

#include <vtkInteractorStyleTrackballCamera.h>

// Project includes: None so far.

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage :\n";
  std::cout << "icnsVonEconomo -I <input vtk mesh> \n\n";
    
  std::cout << "-I <input vtk mesh>            Filename of the mesh to be displayed.\n";
  
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
    
  //std::vector<vtkActor*> Numbers;
  
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
    
    std::cout << "Picked value: " << pickedPosition[0] << " " << pickedPosition[1] << " " << pickedPosition[2] << std::endl;
    
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
  std::cout << "icnsVonEconomo" << std::endl;
  std::cout << "==========================================" << std::endl;
  
  // Initializing parameters with default values:
  
  int c;
  char* inputMeshFilename_1 = NULL;
  char* inputMeshFilename_2 = NULL;
  
  // Reading parameters:
  
  while( (c = getopt( argc, argv, "I:J:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputMeshFilename_1 = optarg;
        std::cout << "  Input image filename 1:            " << inputMeshFilename_1 << std::endl;
        break;
      case 'J':
        inputMeshFilename_2 = optarg;
        std::cout << "  Input image filename 2:            " << inputMeshFilename_2 << std::endl;
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
  
  /*if( inputMeshFilename_1 == NULL )
  {
    std::cerr << "ERROR: No input mesh filename!" << std::endl;
    return EXIT_FAILURE;
  }*/
    
  // DUMMY FUNCTIONALITY:
  
  std::string filename_1 = "/Users/rwerner/Dropbox/ICNS_vonEconomo___private/DatenMichel_NEU/eco01-Julia_M_2016-4-18-16-20-7.vtk";
  std::string filename_2 = "/Users/rwerner/Dropbox/ICNS_vonEconomo___private/DatenMichel_NEU/eco02-Julia_M_2016-5-1-18-47-19.vtk";
  //std::string filename = "/Users/rwerner/Dropbox/ICNS_vonEconomo___private/DatenMichel_NEU/eco01-Julia_M_2016-4-18-16-20-7.vtk";
  //std::string filename = "/Users/rwerner/Dropbox/ICNS_vonEconomo___private/DatenMichel_NEU/eco02-Julia_M_2015-9-26-21-57-6.vtk";
  //std::string filename = "/Users/rwerner/Dropbox/ICNS_vonEconomo___private/DatenMichel_NEU/eco02-Julia_M_2016-4-28-13-42-8.vtk";
  
  std::cout << "FILENAME 1: " << filename_1 << std::endl;
  std::cout << "FILENAME 2: " << filename_2 << std::endl;
  
  vtkSmartPointer<vtkPolyDataReader> inputMeshReader_1 = vtkSmartPointer<vtkPolyDataReader>::New();
  inputMeshReader_1->SetFileName( filename_1.c_str() );
  inputMeshReader_1->Update();
  vtkSmartPointer<vtkPolyDataMapper> mapper_1 = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper_1->SetInputConnection( inputMeshReader_1->GetOutputPort() );
  
  vtkSmartPointer<vtkPolyDataReader> inputMeshReader_2 = vtkSmartPointer<vtkPolyDataReader>::New();
  inputMeshReader_2->SetFileName( filename_2.c_str() );
  inputMeshReader_2->Update();
  vtkSmartPointer<vtkPolyDataMapper> mapper_2 = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper_2->SetInputConnection( inputMeshReader_2->GetOutputPort() );
  
  vtkSmartPointer<vtkLookupTable> lut = vtkLookupTable::New();
  lut->SetNumberOfTableValues(58);
		lut->SetTableRange(0,57);
		lut->Build();
  
		lut->SetTableValue(0, 1, 1, 1, 1); // undefined
		lut->SetTableValue(1, 158.0/255.0, 204/255.0, 46.0, 1); //FA 603
		lut->SetTableValue(2, 1.0/255.0, 173.0/255.0, 75.0/255.0, 1); //FB 610
		lut->SetTableValue(3, 0.0/255.0, 136.0/255.0, 90.0/255.0, 1); //FC 612
		lut->SetTableValue(4, 102.0/255.0, 158.0/255.0, 71.0/255.0, 1); //FCBm 613
		lut->SetTableValue(5, 51.0/255.0, 74.0/255.0, 58.0/255.0, 1); //FD 615
		lut->SetTableValue(6, 34.0/255.0, 84.0/255.0, 59.0/255.0, 1); //FDt 616
		lut->SetTableValue(7, 46.0/255.0, 137.0/255.0, 70.0/255.0, 1); //FDd 617
		lut->SetTableValue(8, 115.0/255.0, 239.0/255.0, 211.0/255.0, 1); //FE 618
		lut->SetTableValue(9, 0.0/255.0, 115.0/255.0, 102.0/255.0, 1); //FF 620
		lut->SetTableValue(10, 75.0/255.0, 80.0/255.0, 58.0/255.0, 1); //FG 645
		lut->SetTableValue(11, 97.0/255.0, 97.0/255.0, 63.0/255.0, 1); //FH 655
		lut->SetTableValue(12, 42.0/255.0, 70.0/255.0, 58.0/255.0, 1); //FJ 663
		lut->SetTableValue(13, 128.0/255.0, 129.0/255.0, 69.0/255.0, 1); //FK 672
		lut->SetTableValue(14, 64.0/255.0, 87.0/255.0, 93.0/255.0, 1); //FL 106
		lut->SetTableValue(15, 81.0/255.0, 90.0/255.0, 97.0/255.0, 1); //FM 149
		lut->SetTableValue(16, 140.0/255.0, 159.0/255.0, 174.0/255.0, 1); //FN 155
		lut->SetTableValue(17, 28.0/255.0, 57.0/255.0, 89.0/255.0, 1); //HA 201
		lut->SetTableValue(18, 36.0/255.0, 50.0/255.0, 79.0/255.0, 1); //HB 209
		lut->SetTableValue(19, 43.0/255.0, 202.0/255.0, 244.0/255.0, 1); //HC 210
		lut->SetTableValue(20, 49.0/255.0, 52.0/255.0, 157.0/255.0, 1); //HD 219
		lut->SetTableValue(21, 0.0/255.0, 123.0/255.0, 154.0/255.0, 1); //HE 221
		lut->SetTableValue(22, 20.0/255.0, 53.0/255.0, 72.0/255.0, 1); //HF 246
		lut->SetTableValue(23, 141.0/255.0, 41.0/255.0, 127.0/255.0, 1); //IA 901
		lut->SetTableValue(24, 91.0/255.0, 51.0/255.0, 124.0/255.0, 1); //IB 902
		lut->SetTableValue(25, 91.0/255.0, 44.0/255.0, 200.0/255.0, 1); //ICD 903
		lut->SetTableValue(26, 0.0/255.0, 221.0/255.0, 240.0/255.0, 1); //LA1 211
		lut->SetTableValue(27, 122.0/255.0, 200.0/255.0, 249.0/255.0, 1); //LA2 212
		lut->SetTableValue(28, 10.0/255.0, 78.0/255.0, 161.0/255.0, 1); //LB 205
		lut->SetTableValue(29, 1.0/255.0, 81.0/255.0, 150.0/255.0, 1); //LC1 265
		lut->SetTableValue(30, 43.0/255.0, 62.0/255.0, 167.0/255.0, 1); //LC2 256
		lut->SetTableValue(31, 38.0/255.0, 42.0/255.0, 51.0/255.0, 1); //LC3 250
		lut->SetTableValue(32, 67.0/255.0, 84.0/255.0, 100.0/255.0, 1); //LD 164
		lut->SetTableValue(33, 199.0/255.0, 195.0/255.0, 245.0/255.0, 1); //LE1 941
		lut->SetTableValue(34, 91.0/255.0, 44.0/255.0, 200.0/255.0, 1); //LE2 903
		lut->SetTableValue(35, 255.0/255.0, 232.0/255.0, 152.0/255.0, 1); //OA 701
		lut->SetTableValue(36, 25.0/255.0, 21.0/255.0, 5.0/255.0, 1); //OB (715: 255,210,5) new
		lut->SetTableValue(37, 255.0/255.0, 236.0/255.0, 1.0/255.0, 1); //OC 717
		//lut->SetTableValue(38, 211.0/255.0, 33.0/255.0, 81.0/255.0, 1); //PA 488
		//lut->SetTableValue(39, 182.0/255.0, 32.0/255.0, 70.0/255.0, 1); //PB1 474
		//lut->SetTableValue(40, 179.0/255.0, 36.0/255.0, 92.0/255.0, 1); //PB2 471
		lut->SetTableValue(38, 255.0/255.0, 0.0/255.0, 0.0/255.0, 1); //PA 488
		lut->SetTableValue(39, 0.0/255.0, 255.0/255.0, 0.0/255.0, 1); //PB1 474
		lut->SetTableValue(40, 255.0/255.0, 185.0/255.0, 15.0/255.0, 1); //PB2 471
		lut->SetTableValue(41, 89.0/255.0, 45.0/255.0, 70.0/255.0, 1); //PC 469
		lut->SetTableValue(42, 144.0/255.0, 32.0/255.0, 57.0/255.0, 1); //PD 467
		lut->SetTableValue(43, 107.0/255.0, 36.0/255.0, 44.0/255.0, 1); //PE 465
		lut->SetTableValue(44, 86.0/255.0, 34.0/255.0, 46.0/255.0, 1); //PF 463
		lut->SetTableValue(45, 194.0/255.0, 28.0/255.0, 66.0/255.0, 1); //PG 460
		lut->SetTableValue(46, 224.0/255.0, 55.0/255.0, 52.0/255.0, 1); //PH 427
		lut->SetTableValue(47, 255.0/255.0, 76.0/255.0, 35.0/255.0, 1); //TA 757
		lut->SetTableValue(48, 255.0/255.0, 125.0/255.0, 37.0/255.0, 1); //TB 764
		lut->SetTableValue(49, 255.0/255.0, 98.0/255.0, 37.0/255.0, 1); //TC 768
		lut->SetTableValue(50, 255.0/255.0, 60.0/255.0, 32.0/255.0, 1); //TD 770
		lut->SetTableValue(51, 191.0/255.0, 159.0/255.0, 108.0/255.0, 1); //TE 815
		lut->SetTableValue(52, 86.0/255.0, 35.0/255.0, 44.0/255.0, 1); //TF 463
		lut->SetTableValue(53, 177.0/255.0, 100.0/255.0, 46.0/255.0, 1); //TG 819
		lut->SetTableValue(54, 103.0/255.0, 68.0/255.0, 49.0/255.0, 1); //TH 836
		lut->SetTableValue(55, 199.0/255.0, 168.0/255.0, 113.0/255.0, 1); //TJK 816
		lut->SetTableValue(56, 190.0/255.0, 20.0/255.0, 100.0/255.0, 1); //TGa
		lut->SetTableValue(57, 10.0/255.0, 60.0/255.0, 90.0/255.0, 1); //FE2
  
		lut->Build();
  
  mapper_1->SetLookupTable( lut );
  mapper_1->ScalarVisibilityOn();
  mapper_1->UseLookupTableScalarRangeOn();
  mapper_1->Update();
  
  mapper_2->SetLookupTable( lut );
  mapper_2->ScalarVisibilityOn();
  mapper_2->UseLookupTableScalarRangeOn();
  mapper_2->Update();
  
  // Create a renderer, render window, and interactor:
  
  // Create an actor
  vtkSmartPointer<vtkActor> actor_1 = vtkSmartPointer<vtkActor>::New();
  actor_1->SetMapper(mapper_1);
  vtkSmartPointer<vtkActor> actor_2 = vtkSmartPointer<vtkActor>::New();
  actor_2->SetMapper(mapper_2);
  //actor->GetProperty()->SetRepresentationToWireframe();
  
  // A renderer and render window
  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  
  // An interactor
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow) ;

  // Create appropriate interactor/interaction style:
  //vtkSmartPointer<MouseInteractorStylePP> style =
  //  vtkSmartPointer<MouseInteractorStylePP>::New();
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  renderWindowInteractor->SetInteractorStyle(style);
  style->SetCurrentRenderer(renderer);
  
  //renderer->AddActor(actor_1);
  renderer->AddActor(actor_2);
  renderer->SetBackground(1,1,1); // Background color white
  //renderer->SetBackground( .1, .3,.2 ); // Background color dark green
  
  renderWindow->Render();
  
  renderWindowInteractor->Start();
  

  /*
  std::string filename = inputMeshFilename_1;
  vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
  reader->SetFileName( filename.c_str() );
  reader->Update();
  
  // -------------------------------------------------------------
  // Loading input mesh(es):
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input mesh 1 ... " << std::flush;
  
  vtkSmartPointer<vtkPolyData> inputMesh_1;
  
  vtkSmartPointer<vtkPolyDataReader> inputMeshReader_1 = vtkSmartPointer<vtkPolyDataReader>::New();
  inputMeshReader_1->SetFileName( inputMeshFilename_1 );
  //inputMeshReader_1->Update();
  //inputMesh_1 = inputMeshReader_1->GetOutput();
  
  std::cout << "OK." << std::endl;
  std::cout << "Loading input mesh 2 ... " << std::flush;
  
  vtkSmartPointer<vtkPolyData> inputMesh_2;
  
  vtkSmartPointer<vtkPolyDataReader> inputMeshReader_2 = vtkSmartPointer<vtkPolyDataReader>::New();
  inputMeshReader_2->SetFileName( inputMeshFilename_2 );
  inputMeshReader_2->Update();
  inputMesh_2 = inputMeshReader_2->GetOutput();
  
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Transforming mesh 2:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Computing transformation matrix ... " << std::flush;
  
  vtkSmartPointer<vtkPoints> sourcePoints =
  vtkSmartPointer<vtkPoints>::New();
  double sourcePoint1[3] = {16.6715, 34.326, 6.67693};
  sourcePoints->InsertNextPoint(sourcePoint1);
  double sourcePoint2[3] = {19.4095, 23.2767, -23.3011};
  sourcePoints->InsertNextPoint(sourcePoint2);
//double sourcePoint3[3] = {15.8685, 52.1864, 6.20179};
  double sourcePoint3[3] = {5.42256, -0.0283235, -19.6914};
  sourcePoints->InsertNextPoint(sourcePoint3);

  
  vtkSmartPointer<vtkPoints> targetPoints =
  vtkSmartPointer<vtkPoints>::New();
  double targetPoint1[3] = {-0.455298, 8.87992, -31.7223};
  targetPoints->InsertNextPoint(targetPoint1);
  double targetPoint2[3] = {27.0038, -9.17666, -24.3916};
  targetPoints->InsertNextPoint(targetPoint2);
  //double targetPoint3[3] = {5.75604, 18.2413, -44.1045};
  double targetPoint3[3] = {33.1301, -23.9376, -19.5321};
  targetPoints->InsertNextPoint(targetPoint3);
  
  vtkSmartPointer<vtkLandmarkTransform> landmarkTransform =
  vtkSmartPointer<vtkLandmarkTransform>::New();
  landmarkTransform->SetSourceLandmarks(sourcePoints);
  landmarkTransform->SetTargetLandmarks(targetPoints);
  landmarkTransform->SetModeToRigidBody();
  landmarkTransform->Update();
  
  // Display the transformation matrix that was computed
  vtkMatrix4x4* mat = landmarkTransform->GetMatrix();
  std::cout << "Matrix: " << *mat;
  
  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
  vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  //transformFilter->SetInputData( inputMesh_1 );
  //transformFilter->SetTransform(landmarkTransform);
  //transformFilter->Update();
  
  // -------------------------------------------------------------
  // Rendering data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Rendering data ... " << std::flush;
  
  // Create a mapper and actor:
  
  vtkSmartPointer<vtkPolyDataMapper> mapper_1 = vtkSmartPointer<vtkPolyDataMapper>::New();
  // mapper->SetInputConnection(cylinderSource->GetOutputPort());
  //mapper_1->SetInputData( inputMesh_2 );
  mapper_1->SetInputConnection( reader->GetOutputPort() );
  
  vtkSmartPointer<vtkActor> actor_1 = vtkSmartPointer<vtkActor>::New();
  actor_1->SetMapper( mapper_1 );
  
  vtkSmartPointer<vtkPolyDataMapper> mapper_2 = vtkSmartPointer<vtkPolyDataMapper>::New();
  // mapper->SetInputConnection(cylinderSource->GetOutputPort());
  //mapper_2->SetInputData( inputMesh_2 );
  mapper_2->SetInputConnection( transformFilter->GetOutputPort() );
  
  vtkSmartPointer<vtkActor> actor_2 = vtkSmartPointer<vtkActor>::New();
  actor_2->SetMapper( mapper_2 );
 
  // Create a renderer, render window, and interactor:
  
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer( renderer );
 
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );
    
  // Create appropriate interactor/interaction style:
    
  vtkSmartPointer<MouseInteractorStylePP> style = vtkSmartPointer<MouseInteractorStylePP>::New();
  renderWindowInteractor->SetInteractorStyle( style );
 
  // Add the actor to the scene:
  
  renderer->AddActor( actor_1 );
  //renderer->AddActor( actor_2 );
  renderer->SetBackground( .1, .3,.2 ); // Background color dark green
 
  /*
  vtkSmartPointer<vtkSphereSource> newSphereSource = vtkSmartPointer<vtkSphereSource>::New();
  newSphereSource->SetRadius( 10 );
  double pickedPosition[3];
  pickedPosition[0] = -11.9899;
  pickedPosition[1] = 16.5337;
  pickedPosition[2] = 42.1059;
  newSphereSource->SetCenter( pickedPosition );
  
  vtkSmartPointer<vtkPolyDataMapper> newMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  newMapper->SetInputConnection( newSphereSource->GetOutputPort() );
  
  vtkSmartPointer<vtkActor> newActor = vtkSmartPointer<vtkActor>::New();
  newActor->SetMapper( newMapper );
  //newActor->GetProperty()->SetColor(1,0,0);
//  this->CurrentRenderer->AddViewProp( newActor );
  renderer->AddActor( newActor );
  */
  
  // Render and interact:
  
  /*
  renderWindow->SetWindowName( argv[0] );
  renderWindow->Render();
  renderWindowInteractor->Start();
   */
 
  return EXIT_SUCCESS;
}
