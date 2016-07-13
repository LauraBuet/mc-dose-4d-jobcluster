/** \file icnsFlipVTKPolyData.cxx
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
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkReverseSense.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

// Project includes:
#include "vtkTexturingHelper.h"

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage :\n";
  std::cout << "icnsFlipVTKPolyData -I <input vtk mesh>  \n\n";
    
  std::cout << "-I <input vtk mesh>            Filename of the mesh to be flipped.\n";
  std::cout << "-O <output vtk mesh>           Filename of the mesh to be written.\n";
    
  std::cout << "-x                             Flip along x-axis.\n";
  std::cout << "-y                             Flip along y-axis.\n";
  std::cout << "-z                             Flip along z-axis.\n";
    
  std::cout << "-h                             Print this help.\n";
  std::cout << "\n";
}


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
    
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsFlipVTKPolyData" << std::endl;
  std::cout << "==========================================" << std::endl;
    
  // Initializing parameters with default values:
    
  int c;
  
  std::string inputMeshFilename  = NULL;
  std::string outputMeshFilename = NULL;
  
  unsigned int flippingAxis[3];
  flippingAxis[0] = 0;
  flippingAxis[1] = 0;
  flippingAxis[2] = 0;
  
  bool isVTKData = true;
  bool isOBJData = false;
    
  // Reading parameters:
    
  while( (c = getopt( argc, argv, "I:O:xyzh?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputMeshFilename = optarg;
        std::cout << "Input mesh filename:   " << inputMeshFilename << std::endl;
        break;
      case 'O':
        outputMeshFilename = optarg;
        std::cout << "Output mesh filename:  " << outputMeshFilename << std::endl;
        break;
      case 'x':
        flippingAxis[0] = 1;
        std::cout << "Flipping along x axis: " << "ON" << std::endl;
        break;
      case 'y':
        flippingAxis[1] = 1;
        std::cout << "Flipping along y axis: " << "ON" << std::endl;
        break;
      case 'z':
        flippingAxis[2] = 1;
        std::cout << "Flipping along z axis: " << "ON" << std::endl;
        break;
      case 'h':
      case '?':
        PrintHelp();
        return EXIT_FAILURE;
      default:
        std::cout << "  Argument " << (char)c << " not processed!\n" << std::endl;
        return EXIT_FAILURE;
    }
  }
    
  // Check validity of arguments:
    
  if( inputMeshFilename.empty() )
  {
    std::cerr << "ERROR: No input mesh filename!" << std::endl;
    return EXIT_FAILURE;
  }
    
  if( outputMeshFilename.empty() )
  {
    std::cerr << "ERROR: No output mesh filename!" << std::endl;
    return EXIT_FAILURE;
  }
    
  if( (flippingAxis[0] == 0) && (flippingAxis[1] == 0) && (flippingAxis[2] == 0) )
  {
    std::cerr << "ERROR: No flipping axis defined!" << std::endl;
    return EXIT_FAILURE;
  }
    
  // Reading input file:
    
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input mesh ... " << std::flush;
    
  vtkSmartPointer<vtkPolyData> inputPolyData;
  vtkSmartPointer<vtkPolyData> outputPolyData;
    
  vtkSmartPointer<vtkPolyDataReader> inputMeshReader;
  vtkTexturingHelper inputMeshHelper;
    
  if( inputMeshFilename.substr( inputMeshFilename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "Reading VTK data ..." << std::endl;
    inputMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    inputMeshReader->SetFileName( inputMeshFilename.c_str() );
    inputMeshReader->Update();
    inputPolyData = inputMeshReader->GetOutput();
  }
  else
  {
    std::cout << "Reading OBJ data ..." << std::endl;
    isVTKData = false;
    isOBJData = true;
    inputMeshHelper.ReadGeometryFile( inputMeshFilename );
    inputPolyData = inputMeshHelper.GetPolyData();
  }
    
  std::cout << "OK." << std::endl;
    
  // Flip input polydata using vtkTransformPolyData.
  // First, define the appropriate transformation:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Flipping data ... " << std::flush;
  
  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->Scale( flippingAxis[0], flippingAxis[1], flippingAxis[2] );
  
  vtkSmartPointer<vtkTransformPolyDataFilter> transformPolyDataFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformPolyDataFilter->SetInputData( inputPolyData );
  transformPolyDataFilter->SetTransform( transform );
  transformPolyDataFilter->Update();
  
  // Then, reverse normals etc.:
  
  vtkSmartPointer<vtkReverseSense> reverseSense = vtkSmartPointer<vtkReverseSense>::New();
  reverseSense->SetInputConnection( transformPolyDataFilter->GetOutputPort() );
  reverseSense->ReverseCellsOn();
  reverseSense->ReverseNormalsOn();
  reverseSense->Update();
  
  outputPolyData = reverseSense->GetOutput();
  
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Finally, write data:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Writing output mesh ... " << std::flush;
  
  vtkSmartPointer<vtkPolyDataWriter> meshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  meshWriter->SetInputData( outputPolyData );
  meshWriter->SetFileName( outputMeshFilename.c_str() );
  meshWriter->Write();
  
  return EXIT_SUCCESS;
}
