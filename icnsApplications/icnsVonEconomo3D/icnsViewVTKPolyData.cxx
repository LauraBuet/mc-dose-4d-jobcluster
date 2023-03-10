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
#include <sstream>
#include <string>
extern "C"
{
#include "getopt.h"
}

// ITK includes: None so far.

// VTK includes:
#include <vtkActor.h>
#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkMatrix4x4.h>
#include <vtkPointData.h>
#include <vtkPolyDataCollection.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRendererCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>

#include <vtkLookupTable.h>

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
  std::cout << "icnsViewVTKPolyData -I <list of vtk meshes> [...] \n"  << std::endl;
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
  std::cout << "-l <list of LUTs>                   LUTs to be assigned to actor mappers (not for texture mode!)." << std::endl;
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

vtkSmartPointer<vtkLookupTable> ReadLUT( std::string lutFilename, std::vector< std::string >& scalarNames );
vtkMatrix4x4* Read4x4Matrix( std::string filename );
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
  
  std::cout << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsViewVTKPolyData"                        << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  
  // Initializing parameters with default values:
  
  std::vector<std::string>   inputFilenames;
  std::vector<std::string>   inputTextureFilePrefixes;
  std::vector<std::string>   inputTextureFilePostfixes;
  std::vector<unsigned int>  inputNumberOfTextureFiles;
  std::vector<std::string>   inputTransforms;
  
  std::vector<std::string>   labelNames;
  
  std::vector<std::string>   actorColors;
  std::vector<std::string>   mapperLUTs;
  
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
  while( (c = getopt( argc, argv, "I:T:t:n:e:c:l:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputFilenames.push_back( argv[optind] );
          std::cout << "  Mesh [" << cnt << "]:                    " << argv[optind] << std::endl;
        }
        break;
      case 'T':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTransforms.push_back( argv[optind] );
          std::cout << "  Transform [" << cnt << "]:               " << argv[optind] << std::endl;
        }
        break;
      case 't':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTextureFilePrefixes.push_back( argv[optind] );
          std::cout << "  Texture file prefix [" << cnt << "]:     " << argv[optind] << std::endl;
        }
        break;
      case 'e':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTextureFilePostfixes.push_back( argv[optind] );
          std::cout << "  Texture file postfix [" << cnt << "]:    " << argv[optind] << std::endl;
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
      case 'l':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          mapperLUTs.push_back( argv[optind] );
          std::cout << "  Mapper LUT [" << cnt << "]:              " << argv[optind] << std::endl;
        }
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_SUCCESS;
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
      
      if( ( !inputTransforms.empty() ) &&
          ( std::strcmp( inputTransforms[iFilenames].c_str(), "NONE" ) != 0 ) )
      {
        vtkSmartPointer<vtkMatrix4x4> currentTransformMatrix = Read4x4Matrix( inputTransforms[iFilenames] );
        
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
        vtkSmartPointer<vtkMatrix4x4> currentTransformMatrix = Read4x4Matrix( inputTransforms[iFilenames] );
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
      
      vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
      actor->SetMapper( mapper );
      
      if( ( !mapperLUTs.empty() ) &&
          ( std::strcmp( mapperLUTs[iFilenames].c_str(), "NONE" ) != 0 ) )
      {
        std::cout << "Reading LUT from file: " << mapperLUTs[iFilenames] << " ... " << std::endl;
        vtkSmartPointer<vtkLookupTable> lut = ReadLUT( mapperLUTs[iFilenames], labelNames );
        
        mapper->SetLookupTable( lut );
        mapper->ScalarVisibilityOn();
        mapper->UseLookupTableScalarRangeOn();
        mapper->Update();
        
        // Checking if LUT range matches scalar range of mesh:
        
        std::map<unsigned int, unsigned int> labelStats;
        vtkSmartPointer<vtkPolyData> currentMesh = static_cast< vtkPolyData* >( polyDataCollection->GetItemAsObject( iFilenames ) );
        for( vtkIdType iPointID = 0; iPointID < currentMesh->GetNumberOfPoints(); iPointID++ )
        {
          double pointCoord[3];
          currentMesh->GetPoint( iPointID, pointCoord );
          //std::cout << "Point " << iPointID << " : (" << pointCoord[0] << " " << pointCoord[1] << " " << pointCoord[2] << ")" << std::endl;
          
          unsigned int pointScalarValue = currentMesh->GetPointData()->GetArray(0)->GetTuple1( iPointID );
          auto mapIterator = labelStats.find( pointScalarValue );
          
          // Scalar existing in map? If so, increment entry:
          if( (mapIterator != labelStats.end()) && (!labelStats.empty()) )
          {
            labelStats[ pointScalarValue ] = mapIterator->second +1;
          }
          else
          {
            labelStats[ pointScalarValue ] = 1;
          }
          mapIterator = labelStats.find( pointScalarValue );
          //std::cout << "  LABEL " << mapIterator->first << ": " << mapIterator->second << std::endl;
        }
        
        //std::cout << "In total " << labelStats.size() << " labels!" << std::endl;
        //auto mapIterator = labelStats.begin();
        //for( unsigned int iLabels = 0; iLabels < labelStats.size(); iLabels++, mapIterator++ );
        //{
        //  std::cout << "  LABEL " << mapIterator->first << ": " << mapIterator->second << std::endl;
        //}
        /*std::cout << "Distribution of scalars:" << std::endl;
        
        for( unsigned int i=0; i < pointsPerID.size(); i++ )
        {
          std::cout << "  LABEL " << i << ": " << pointsPerID[i] << " points" << std::endl;
        }
        
        unsigned int nScalarArrays = outputMesh->GetPointData()->GetNumberOfArrays();
        if( nScalarArrays > 0 )
        {
          std::cout << "Point data associated with " << nScalarArrays << " arrays." << std::endl;
          std::cout << "Scalar range: [" << outputMesh->GetScalarRange()[0] << "; " <<
          outputMesh->GetScalarRange()[1] << "] " << std::endl;
          std::cout << "Number of entries: " << outputMesh->GetPointData()->GetArray(0)->GetNumberOfTuples() << std::endl;
        }
        else
        {
          std::cout << "No scalar data available." << std::endl;
        }*/
      }
      else if ( !actorColors.empty() )
      {
        std::vector<float> colorRGB = ConvertColorStringToRGB( actorColors[iFilenames] );
        actor->GetProperty()->SetColor( colorRGB[0], colorRGB[1], colorRGB[2] );
      }
      else
      {
        mapper->ScalarVisibilityOff();
      }
    
      renderer->AddActor( actor );
    }
  }
  
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
  
  // Generating callback to reset actor transformations on return press:
  
  vtkSmartPointer<vtkCallbackCommand> clickCallback = vtkSmartPointer<vtkCallbackCommand>::New();
  clickCallback->SetCallback( ClickCallbackFunction );
  renderWindowInteractor->AddObserver( vtkCommand::KeyPressEvent, clickCallback );
  
  renderer->SetBackground( 1.0, 1.0, 1.0 ); // Background color white
  //renderer->GradientBackgroundOn();
  //renderer->SetBackground( 1, 1, 1);
  //renderer->SetBackground2( 0, 0, 0);
  renderer->ResetCamera();
  
  renderWindow->SetSize( 800, 600 );
  renderWindow->Render();
  
  // -------------------------------------------------------------
  // Visualizing data (rendering):
  
  // This starts the event loop and as a side effect causes an initial render:
  std::cout << "Starting visualization!" << std::endl;
  renderWindowInteractor->Start();
  
  return EXIT_SUCCESS;
}

// ---------------------------------------------------------------
// Definition of helper functions:
// ---------------------------------------------------------------

// ReadLUT: reads LUT in "filename".
// Requires the LUT to be of freesurfer format.

vtkSmartPointer<vtkLookupTable> ReadLUT( std::string filename, std::vector< std::string >& scalarNames )
{
  std::string line;
  ifstream txtFile( filename );
  
  // First run: Fill vector of float vectors with valid entries
  
  std::vector< std::vector< std::string > > LUTEntries;
  
  while( txtFile )
  {
    if( !getline( txtFile, line ) ) break;
    if( line.length() == 0 ) continue;
    
    // Analyze current line:
    std::istringstream currentLineStringStream( line );
    std::vector< std::string > currentLUTLine;
    
    if( std::strcmp( line.substr(0, 1).c_str(), "#" ) != 0 ) // comment lines
    {
      char delimiter = ' ';
      
      // Analyzing string:
      while( currentLineStringStream )
      {
        std::string lutLineEntryAsString;
        
        // Analyze exceptional cases:
        if( !getline( currentLineStringStream, lutLineEntryAsString, delimiter ) ) break;
        if( std::strcmp( lutLineEntryAsString.substr(0, 1).c_str(), " " ) == 0 ) continue;
        if( lutLineEntryAsString.empty() ) continue;
        if( std::strcmp( lutLineEntryAsString.substr(0, 1).c_str(), "#" ) == 0 ) break;
        
        // Otherwise:
        std::cout << lutLineEntryAsString << " x ";
        currentLUTLine.push_back( lutLineEntryAsString );
      }
      LUTEntries.push_back( currentLUTLine );
      std::cout << std::endl;
    }
  }
  txtFile.close();
  
  // Fill label names vector:
  
  for( unsigned int i = 0; i < LUTEntries.size(); i++ )
  {
    scalarNames.push_back( LUTEntries[i][1] );
  }
  
  // Set up LUT and fill it:
  
  vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
  lut->SetNumberOfTableValues( LUTEntries.size() );
  lut->SetTableRange( 0, LUTEntries.size()-1 );
	lut->Build();
  
  for( unsigned int i = 0; i < LUTEntries.size(); i++ )
  {
    std::cout << atoi( LUTEntries[i][0].c_str() ) << " -- " << std::flush;
    std::cout << scalarNames[i] << ": " << std::flush;
    std::cout << atof( LUTEntries[i][2].c_str() )/255.0 << " " <<
                 atof( LUTEntries[i][3].c_str() )/255.0 << " " <<
                 atof( LUTEntries[i][4].c_str() )/255.0 << " " <<
                 atof( LUTEntries[i][5].c_str() ) << std::endl;
  }
  
  for( unsigned int i = 0; i < LUTEntries.size(); i++ )
  {
    lut->SetTableValue( atoi( LUTEntries[i][0].c_str() ),
                        atof( LUTEntries[i][2].c_str() )/255.0,
                        atof( LUTEntries[i][3].c_str() )/255.0,
                        atof( LUTEntries[i][4].c_str() )/255.0,
                        atof( LUTEntries[i][5].c_str() ) );
  }
  
  lut->Build();
  return lut;
}


vtkMatrix4x4* Read4x4Matrix( std::string filename )
{
  std::cout << "Reading matrix from file: " << filename << "..." << std::endl;
  vtkMatrix4x4* mat = vtkMatrix4x4::New();
  
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
    //std::cout << "----------------------------------------" << std::endl;
    std::cout << "Reset actor transform" << std::endl;
    
    // Get current renderer from caller data to extract desired actors:
    
    vtkSmartPointer<vtkRenderer> currentRenderer = iren->GetRenderWindow()->GetRenderers()->GetFirstRenderer();
    std::cout << "No Renderers: " << iren->GetRenderWindow()->GetRenderers()->GetNumberOfItems() << std::endl;
    
    unsigned int nActors = currentRenderer->GetActors()->GetNumberOfItems();
    
    // Iterate over actors:
    
    for( unsigned int iActors = 0; iActors < currentRenderer->GetActors()->GetNumberOfItems(); iActors++ )
    {
      vtkSmartPointer<vtkActor> currentActor       = static_cast<vtkActor*>( currentRenderer->GetActors()->GetItemAsObject( iActors ) );
      vtkMatrix4x4* currentTransformMatrix         = currentActor->GetMatrix();
      
      // Printing matrix:
      for( unsigned int i = 0; i < 4; i++ )
      {
        for( unsigned int j = 0; j < 4; j++ )
        {
          if( i==j ) currentTransformMatrix->SetElement( i, j, 1.0f );
          if( i!=j ) currentTransformMatrix->SetElement( i, j, 0.0f );
        }
      }
      currentTransformMatrix->Modified();
    }
    
    currentRenderer->Render();
    iren->GetRenderWindow()->Render();
  } // end of callback "RETURN"
} // end of ClickCallbackFunction()
