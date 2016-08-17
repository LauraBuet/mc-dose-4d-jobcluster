/** \file icnsViewVTKPolyData.cxx
 *
 *  Purpose: 
 *  Load and view (a list of) vtk polydata structure. Also works
 *  with obj files, either single or multi-texture. For the applied mutli-
 *  texture obj-support used see
 *
 *  https://github.com/Scylardor/vtkTexturingHelper
 *
 *  and
 *
 *  http://scylardor.fr/2013/05/06/making-multi-texturing-work-with-vtk/
 *
 *  Optional: Use structure-specific transformations.
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
#include <vtkCamera.h>
#include <vtkMatrix4x4.h>
#include <vtkPolyDataCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <vtkInteractorStyleSwitch.h>
#include <vtkMatrixToLinearTransform.h>

// Project includes:
#include "vtkTexturingHelper.h"

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsViewVTKPolyData -I <list of vtk mesh> [...] \n"  << std::endl;
  std::cout << std::endl;
  std::cout << "-I <list of vtk meshes>             Filename list of meshes to be displayed." << std::endl;
  std::cout << "                                    Possible endings: VTK or OBJ." << std::endl;
  std::cout << "-T <list of transformations>        Filename list of transformations (4x4 matrix)" << std::endl;
  std::cout << "                                    to be applied to meshes." << std::endl;
  std::cout << "-t <list of texture file prefixes>  List of texture prefixes (e.g. ..._eco06_hell)." << std::endl;
  std::cout << "-e <list of texture file prefixes>  List of texture postfixes (e.g. .png)." << std::endl;
  std::cout << "-n <number of texture files>        Number of texture files." << std::endl;
  std::cout << std::endl;
  std::cout << "-c <list of colors>                 Colors to be assigned to actors (not for texture mode!)." << std::endl;
  std::cout << "-h                                  Print this help." << std::endl;
  std::cout << "\n";
}

// ---------------------------------------------------------------
// Additional routine declarations:
// ---------------------------------------------------------------

vtkMatrix4x4* Load4x4Matrix( std::string filename );
std::vector<float> ConvertColorStringToRGB( std::string colorString );

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
  std::cout << "icnsViewVTKPolyData" << std::endl;
  std::cout << "==========================================" << std::endl;
  
  // Initializing parameters with default values:
  
  std::vector<std::string>   inputFilenames;
  std::vector<std::string>   inputTextureFilePrefixes;
  std::vector<std::string>   inputTextureFilePostfixes;
  std::vector< unsigned int> inputNumberOfTextureFiles;
  std::vector<std::string>   inputTransforms;
  
  std::vector<std::string>   actorColors;
  
  bool readingOBJdata = false;
  
  // Initializing renderer:
  
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  
  // Reading parameters.
  // Note that the filnames of the lists are considered to *not* contain
  // a '-', because this is assumed to be reserved and to start the next
  // parameter readout.
  
  int c;
  int cnt = 0;
  
  std::cout << "Reading parameters ..." << std::endl;
  while( (c = getopt( argc, argv, "I:T:t:n:e:c:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputFilenames.push_back( argv[optind] );
          std::cout << "  Mesh [" << cnt << "]:                  " << argv[optind] << std::endl;
        }
        break;
      case 'T':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTransforms.push_back( argv[optind] );
          std::cout << "  Transform [" << cnt << "]:             " << argv[optind] << std::endl;
        }
        break;
      case 't':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTextureFilePrefixes.push_back( argv[optind] );
          std::cout << "  Texture file prefix [" << cnt << "]:   " << argv[optind] << std::endl;
        }
        break;
      case 'e':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTextureFilePostfixes.push_back( argv[optind] );
          std::cout << "  Texture file postfix [" << cnt << "]:  " << argv[optind] << std::endl;
        }
        break;
      case 'n':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputNumberOfTextureFiles.push_back( atoi(argv[optind]) );
          std::cout << "  Number of texture files [" << cnt << "]: " << argv[optind] << std::endl;
        }
        break;
      case 'c':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          actorColors.push_back( argv[optind] );
          std::cout << "  Actor color [" << cnt << "]:             " << argv[optind] << std::endl;
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
  
  // -------------------------------------------------------------
  // Check input:
  
  if( inputFilenames.empty() )
  {
    std::cerr << "ERROR: No input filename!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading mesh data:
  
  vtkSmartPointer<vtkPolyDataCollection> polyDataCollection = vtkSmartPointer<vtkPolyDataCollection>::New();
  vtkSmartPointer<vtkPolyDataReader> polyDataReader;
  
  for( unsigned iFilenames = 0; iFilenames < inputFilenames.size(); iFilenames++ )
  {
    std::string currentFilename = inputFilenames[iFilenames];
    vtkSmartPointer<vtkPolyData> currentPolyData = vtkSmartPointer<vtkPolyData>::New();
    
    // Extracting file ending and choosing approproate reader:
    
    if( currentFilename.substr( currentFilename.find_last_of(".") + 1) == "vtk" )
    {
      std::cout << "Reading VTK data: " << currentFilename << std::endl;
      polyDataReader = vtkSmartPointer<vtkPolyDataReader>::New();
      polyDataReader->SetFileName( currentFilename.c_str() );
      polyDataReader->Update();
    
      currentPolyData->DeepCopy( polyDataReader->GetOutput() );
      polyDataCollection->AddItem( currentPolyData );
      
      if( !inputTransforms.empty() && std::strcmp( inputTransforms[iFilenames].c_str(), "NONE" ) )
      {
        vtkSmartPointer<vtkMatrix4x4> currentTransformMatrix = Load4x4Matrix( inputTransforms[iFilenames] );
        
        vtkSmartPointer<vtkMatrixToLinearTransform> meshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
        meshTransform->SetInput( currentTransformMatrix );
        
        vtkSmartPointer<vtkTransformPolyDataFilter> meshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        meshTransformFilter->SetInputData( static_cast<vtkPolyData*>( polyDataCollection->GetItemAsObject( iFilenames ) ) );
        meshTransformFilter->SetTransform( meshTransform );
        meshTransformFilter->Update();
        
        polyDataCollection->ReplaceItem( iFilenames, meshTransformFilter->GetOutput() );
      }
    }
    else
    {
      std::cout << "Reading OBJ data ..." << currentFilename << std::endl;
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
      
      if( !inputTransforms.empty() )
      {
        vtkSmartPointer<vtkMatrix4x4> currentTransformMatrix = Load4x4Matrix( inputTransforms[iFilenames] );
        vtkSmartPointer<vtkMatrixToLinearTransform> meshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
        meshTransform->SetInput( currentTransformMatrix );
        
        objReader.GetActor()->SetUserTransform( meshTransform );
      }
      
      // Directly add actor to renderer:
      //currentPolyData->DeepCopy( objReader.GetPolyData() );
      renderer->AddActor( objReader.GetActor() );
    }
  }
  
  // -------------------------------------------------------------
  // Setting up data pipeline:
 
  // Generating mapper and actor for each input file and
  // adding each actor to renderer:
  
  if( !readingOBJdata )
  {
    for( unsigned iFilenames = 0; iFilenames < inputFilenames.size(); iFilenames++ )
    {
      vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      mapper->SetInputData( static_cast< vtkPolyData* >( polyDataCollection->GetItemAsObject( iFilenames ) ) );
      //mapper->SetInputConnection( inputMeshReader->GetOutputPort() );
      mapper->ScalarVisibilityOff();
    
      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper( mapper );
      
      if( !actorColors.empty() )
      {
        std::vector<float> colorRGB = ConvertColorStringToRGB( actorColors[iFilenames] );
        actor->GetProperty()->SetColor( colorRGB[0], colorRGB[1], colorRGB[2] );
      }
    
      renderer->AddActor( actor );
    }
  }
  
/*  vtkSmartPointer<vtkLookupTable> lut = vtkLookupTable::New();
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
 */
  
  // Generating renderWindow and adding renderer to window:
  
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  
  // An interactor:
  
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow( renderWindow );
  
  // Create appropriate interactor/interaction style:
  //vtkSmartPointer<MouseInteractorStylePP> style =
  //  vtkSmartPointer<MouseInteractorStylePP>::New();
  //vtkSmartPointer<vtkInteractorStyleTrackballCamera> style = vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
  //renderWindowInteractor->SetInteractorStyle(style);
  //style->SetCurrentRenderer(renderer);
  
  //renderer->AddActor( actor );
  renderer->SetBackground( 1.0, 1.0, 1.0 ); // Background color white
  renderer->GradientBackgroundOn();
  renderer->SetBackground(1,1,1);
  renderer->SetBackground2(0,0,0);
  renderer->ResetCamera();
  
  renderWindow->SetSize( 800, 600 );
  
  //vtkSmartPointer<vtkOBJExporter> exporter = vtkSmartPointer<vtkOBJExporter>::New();
  //exporter->SetInput( renderWindow );
  //exporter->SetFilePrefix( "/Users/rwerner/Documents/2016_VonEconomo3D/TEST" );
  //exporter->Write();
  
  renderWindow->Render();
  
  // -------------------------------------------------------------
  // Visualizing data (rendering):
  
  // This starts the event loop and as a side effect causes an initial render:
  std::cout << "Starting visualization!" << std::endl;
  renderWindowInteractor->Start();
  
  // Cleaning up:
  //renderer->Delete();
  //renderWindow->Delete();
  //renderWindowInteractor->Delete();
  
  return EXIT_SUCCESS;
}

// ---------------------------------------------------------------
// Definition of helper functions:
// ---------------------------------------------------------------

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

std::vector<float> ConvertColorStringToRGB( std::string colorString )
{
  std::vector<float> RGBValue(3);
  
  if( std::strcmp( colorString.c_str(), "r" ) == 0 )
  {
    RGBValue[0] = 1.0;
    RGBValue[1] = 0.0;
    RGBValue[2] = 0.0;
  }
  else if( std::strcmp( colorString.c_str(), "g" ) == 0 )
  {
    RGBValue[0] = 0.0;
    RGBValue[1] = 1.0;
    RGBValue[2] = 0.0;
  }
  if( std::strcmp( colorString.c_str(), "b" ) == 0 )
  {
    RGBValue[0] = 0.0;
    RGBValue[1] = 0.0;
    RGBValue[2] = 1.0;
  }
  
  return RGBValue;
}
