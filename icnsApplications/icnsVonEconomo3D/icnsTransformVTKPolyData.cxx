/** \file icnsTransformVTKPolyData.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2017 Department of Computational Neuroscience,
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
#include <vtkCleanPolyData.h>
#include <vtkMatrix4x4.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkPolyData.h>
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
  std::cout << "icnsTransformPolyData -I <vtk mesh> -T <transform> -O <output mesh> \n"  << std::endl;
  std::cout << std::endl;
  std::cout << "-I <input vtk meshes>               Filename of mesh to transform." << std::endl;
  std::cout << "-T <input transformation>           Transformation to be applied." << std::endl;
  std::cout << "-O <output vtk mesh>                Filename of transformed mesh." << std::endl;
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
  std::cout << "icnsTransformVTKPolyData"                    << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  
  // -------------------------------------------------------------
  // Initializing parameters with default values:
  
  std::string inputMeshFilename;
  std::string inputTransformationFilename;
  std::string outputMeshFilename;
  
  // Reading parameters:

  int c;
  unsigned int cnt = 0;
  
  std::cout << "Reading parameters ..." << std::endl;
  while( (c = getopt( argc, argv, "I:T:O:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputMeshFilename = optarg;
        std::cout << "Output mesh filename:        " << inputMeshFilename << std::endl;
        break;
      case 'T':
        inputTransformationFilename = optarg;
        std::cout << "Output mesh filename:        " << inputTransformationFilename << std::endl;
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
  
  if( inputMeshFilename.empty() )
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
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Loading input mesh ... "                    << std::flush;
  
  vtkSmartPointer<vtkPolyDataReader> polyDataReader = vtkSmartPointer<vtkPolyDataReader>::New();
  polyDataReader->SetFileName( inputMeshFilename.c_str() );
  polyDataReader->Update();
  
  vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
  polyData = polyDataReader->GetOutput();
  
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Loading input transformation:
  
  vtkSmartPointer<vtkMatrix4x4> transformMatrix = Read4x4Matrix( inputTransformationFilename );
  
  // -------------------------------------------------------------
  // Applying transformation:
  
  vtkSmartPointer<vtkMatrixToLinearTransform> meshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
  meshTransform->SetInput( transformMatrix );
  
  vtkSmartPointer<vtkTransformPolyDataFilter> meshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  meshTransformFilter->SetInputData( polyData );
  meshTransformFilter->SetTransform( meshTransform );
  meshTransformFilter->Update();
  
  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Writing output mesh ... "                   << std::flush;

  vtkSmartPointer<vtkPolyDataWriter> meshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  meshWriter->SetFileName( outputMeshFilename.c_str() );
  meshWriter->SetInputData( meshTransformFilter->GetOutput() );
  meshWriter->Write();
    
  std::cout << "OK." << std::endl;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "icnsTransformVTKPolyData FINISHED"          << std::endl;
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
