/** \file icnsAppendVTKPolyData.cxx
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

// ITK includes: NONE SO FAR.

// VTK includes:
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkMatrix4x4.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkPolyData.h>
#include <vtkPolyDataCollection.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkTransformPolyDataFilter.h>

// Project includes: NONE SO FAR.

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsAppendPolyData -I <list of vtk meshes> -O <output mesh> \n"  << std::endl;
  std::cout << std::endl;
  std::cout << "-I <list of vtk meshes>             Filename list of meshes to be fused." << std::endl;
  std::cout << "-T <list of transformations>        Filename list of mesh transformations." << std::endl;
  std::cout << "-O <output vtk mesh>                Filename of output mesh." << std::endl;
  std::cout << std::endl;
  std::cout << "-h                                  Print this help." << std::endl;
  std::cout << std::endl;
}

// ---------------------------------------------------------------
// Additional routine declarations:
// ---------------------------------------------------------------

vtkMatrix4x4* Read4x4Matrix( std::string filename );

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
  std::cout << "icnsAppendVTKPolyData"                    << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  
  // -------------------------------------------------------------
  // Initializing parameters with default values:
  
  std::vector<std::string>   inputMeshFilenames;
  std::vector<std::string>   inputTransformationFilenames;
  std::string                outputMeshFilename;
  
  // Reading parameters:

  int c;
  unsigned int cnt = 0;
  
  std::cout << "Reading parameters ..." << std::endl;
  while( (c = getopt( argc, argv, "I:T:O:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputMeshFilenames.push_back( argv[optind] );
          std::cout << "Mesh [" << cnt << "]:                    " << argv[optind] << std::endl;
        }
        break;
      case 'T':
        optind--;
        cnt = 0;
        for( ; optind < argc && *argv[optind] != '-'; optind++, cnt++ )
        {
          inputTransformationFilenames.push_back( argv[optind] );
          std::cout << "Transform [" << cnt << "]:               " << argv[optind] << std::endl;
        }
        break;
      case 'O':
        outputMeshFilename = optarg;
        std::cout << "Output mesh filename:        " << outputMeshFilename << std::endl;
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
  
  // Check validity of arguments:
  
  if( inputMeshFilenames.empty() )
  {
    std::cerr << "ERROR: No input mesh specified!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputMeshFilename.empty() )
  {
    std::cerr << "ERROR: No output filename specified!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input meshes:
  
  vtkSmartPointer<vtkPolyDataCollection> inputMeshes = vtkSmartPointer<vtkPolyDataCollection>::New();
  vtkSmartPointer<vtkPolyDataReader> polyDataReader;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Loading input meshes ... "                  << std::endl;
  
  for( unsigned iFilenames = 0; iFilenames < inputMeshFilenames.size(); iFilenames++ )
  {
    std::string currentFilename = inputMeshFilenames[iFilenames];
    vtkSmartPointer<vtkPolyData> currentPolyData = vtkSmartPointer<vtkPolyData>::New();
    
    // Reading vtk data:
    std::cout << "VTK data:  " << currentFilename << std::endl;
    
    polyDataReader = vtkSmartPointer<vtkPolyDataReader>::New();
    polyDataReader->SetFileName( currentFilename.c_str() );
    polyDataReader->Update();
      
    currentPolyData->DeepCopy( polyDataReader->GetOutput() );
    inputMeshes->AddItem( currentPolyData );
    
    // If input transforms specified, apply transformation to meshed:
    if( ( !inputTransformationFilenames.empty() ) &&
        ( std::strcmp( inputTransformationFilenames[iFilenames].c_str(), "NONE" ) != 0 ) )
    {
      std::cout << "Transform: " << currentFilename << std::endl;

      vtkSmartPointer<vtkMatrix4x4> currentTransformMatrix = Read4x4Matrix( inputTransformationFilenames[iFilenames] );
        
      vtkSmartPointer<vtkMatrixToLinearTransform> meshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
      meshTransform->SetInput( currentTransformMatrix );
        
      vtkSmartPointer<vtkTransformPolyDataFilter> meshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
      meshTransformFilter->SetInputData( static_cast<vtkPolyData*>( inputMeshes->GetItemAsObject( iFilenames ) ) );
      meshTransformFilter->SetTransform( meshTransform );
      meshTransformFilter->Update();
        
      inputMeshes->ReplaceItem( iFilenames, meshTransformFilter->GetOutput() );
    }
  }
  
  // -------------------------------------------------------------
  // Appending loaded vtk meshed:
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Running append filter ... "                 << std::flush;
  
  vtkSmartPointer<vtkPolyData> outputMesh = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkAppendPolyData> appendFilter = vtkSmartPointer<vtkAppendPolyData>::New();
  
  for( unsigned int iInputMeshes = 0; iInputMeshes < inputMeshes->GetNumberOfItems(); iInputMeshes++ )
  {
    appendFilter->AddInputData( static_cast< vtkPolyData* >( inputMeshes->GetItemAsObject( iInputMeshes ) ) );
  }
  appendFilter->Update();
  
  std::cout << "Removing duplicate points ... "             << std::flush;
  
  vtkSmartPointer<vtkCleanPolyData> cleanFilter = vtkSmartPointer<vtkCleanPolyData>::New();
  cleanFilter->SetInputConnection( appendFilter->GetOutputPort() );
  cleanFilter->Update();
  
  outputMesh->DeepCopy( cleanFilter->GetOutput() );
  
  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Writing output mesh ... "                   << std::flush;

  vtkSmartPointer<vtkPolyDataWriter> meshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  meshWriter->SetFileName( outputMeshFilename.c_str() );
  meshWriter->SetInputData( outputMesh );
  meshWriter->Write();
    
  std::cout << "OK." << std::endl;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "icnsAppendVTKPolyData FINISHED"             << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


// ---------------------------------------------------------------
// Additional routines definitions:
// ---------------------------------------------------------------

vtkMatrix4x4* Read4x4Matrix( std::string filename )
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
