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
#include <vtkMatrixToLinearTransform.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkReverseSense.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

// Project includes:
#include "vtkTexturingHelper.h"
#include "vtkOBJWriter.h"

// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << "\n";
  std::cout << "Usage :\n";
  std::cout << "icnsFlipVTKPolyData -I <input vtk mesh> -O <output vtk mesh> -T <transformation>\n" << std::endl;
    
  std::cout << "-I <input vtk mesh>            Filename of the mesh to be flipped." << std::endl;
  std::cout << "-O <output vtk mesh>           Filename of the mesh to be written." << std::endl;
  std::cout << "-W <output transformation>     Filename of transformation to be written." << std::endl;
    
  std::cout << "-x                             Flip along x-axis." << std::endl;
  std::cout << "-y                             Flip along y-axis." << std::endl;
  std::cout << "-z                             Flip along z-axis." << std::endl;
    
  std::cout << "-h                             Print this help." << std::endl;
  std::cout << std::endl;
}

// ---------------------------------------------------------------
// Additional routines:
// ---------------------------------------------------------------

vtkMatrix4x4* Load4x4Matrix( std::string filename );

// ---------------------------------------------------------------
// Main routine:
// ---------------------------------------------------------------

int main( int argc, char *argv[] )
{
  if( argc < 4 )
  {
    PrintHelp();
    return EXIT_FAILURE;
  }
  
  std::cout << " "                                        << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << "icnsFlipVTKPolyData"                      << std::endl;
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Reading parameters ..."                   << std::endl;
    
  // Initializing parameters with default values:
  
  std::string inputMeshFilename;
  std::string inputTransformationFilename;
  std::string outputMeshFilename;
  std::string outputTransformationFilename;
  
  int flippingAxis[3];
  flippingAxis[0] = 1;
  flippingAxis[1] = 1;
  flippingAxis[2] = 1;
  
  bool isVTKData = true;
  bool isOBJData = false;
    
  // Reading parameters:
  
  int c;
  while( (c = getopt( argc, argv, "I:T:O:W:xyzh?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputMeshFilename = optarg;
        std::cout << "Input mesh filename:   " << inputMeshFilename << std::endl;
        break;
      case 'T':
        inputTransformationFilename = optarg;
        std::cout << "Input transformation:  " << inputTransformationFilename << std::endl;
        break;
      case 'O':
        outputMeshFilename = optarg;
        std::cout << "Output mesh filename:  " << outputMeshFilename << std::endl;
        break;
      case 'W':
        outputTransformationFilename = optarg;
        std::cout << "Output transformation: " << outputTransformationFilename << std::endl;
        break;
      case 'x':
        flippingAxis[0] = -1;
        std::cout << "Flipping along x axis: " << "ON" << std::endl;
        break;
      case 'y':
        flippingAxis[1] = -1;
        std::cout << "Flipping along y axis: " << "ON" << std::endl;
        break;
      case 'z':
        flippingAxis[2] = -1;
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
    std::cerr << "ERROR: No input mesh filename given!" << std::endl;
    return EXIT_FAILURE;
  }
    
  if( outputMeshFilename.empty() && outputTransformationFilename.empty() )
  {
    std::cerr << "ERROR: No output filename given!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Reading input file:
    
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Loading input mesh ... " << std::flush;
    
  vtkSmartPointer<vtkPolyData> inputPolyData = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkMatrix4x4> inputTransformationMatrix;
  vtkSmartPointer<vtkPolyData> outputPolyData = vtkSmartPointer<vtkPolyData>::New();
    
  vtkSmartPointer<vtkPolyDataReader> inputMeshReader;
  vtkTexturingHelper inputMeshHelper;
    
  if( inputMeshFilename.substr( inputMeshFilename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "[VTK data] ... " << std::flush;
    inputMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    inputMeshReader->SetFileName( inputMeshFilename.c_str() );
    inputMeshReader->Update();
    
    inputPolyData->ShallowCopy( inputMeshReader->GetOutput() );
  }
  else
  {
    std::cout << "[OBJ data] ... " << std::flush;
    isVTKData = false;
    isOBJData = true;
    inputMeshHelper.ReadGeometryFile( inputMeshFilename );
    
    inputPolyData->ShallowCopy( inputMeshHelper.GetPolyData() );
  }
    
  std::cout << "OK." << std::endl;
  
  // Loading specified transforms:
  
  if( !inputTransformationFilename.empty() )
  {
    inputTransformationMatrix = Load4x4Matrix( inputTransformationFilename );
    vtkSmartPointer<vtkMatrixToLinearTransform> inputMeshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
    
    inputMeshTransform->SetInput( inputTransformationMatrix );
      
    vtkSmartPointer<vtkTransformPolyDataFilter> meshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    meshTransformFilter->SetInputData( inputPolyData);
    meshTransformFilter->SetTransform( inputMeshTransform );
    meshTransformFilter->Update();
      
    inputPolyData = meshTransformFilter->GetOutput();
  }
  
  // -------------------------------------------------------------
  // Flip input polydata using vtkTransformPolyData.
  // First, define the appropriate transformation:
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "Flipping data ... "                       << std::flush;
  
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
  
  outputPolyData->ShallowCopy( reverseSense->GetOutput() );
  
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Finally, write data:
  
  if( !outputMeshFilename.empty() )
  {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Writing output mesh ... " << std::flush;
    
    if( outputMeshFilename.substr( outputMeshFilename.find_last_of(".") + 1) == "vtk" )
    {
      std::cout << "[VTK data] ..." << std::flush;
    
      vtkSmartPointer<vtkPolyDataWriter> vtkWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
      vtkWriter->SetInputData( outputPolyData );
      vtkWriter->SetFileName( outputMeshFilename.c_str() );
      vtkWriter->Write();
      
      std::cout << "OK." << std::endl;
    }
    else
    {
      std::cout << "[OBJ data] ..." << std::flush;
    
      vtkSmartPointer<vtkOBJWriter> objWriter = vtkSmartPointer<vtkOBJWriter>::New();
      objWriter->SetInputData( outputPolyData );
      objWriter->SetFileName( outputMeshFilename.c_str() );
      objWriter->Update();
      
      std::cout << "OK." << std::endl;
    }
  }
  
  if( !outputTransformationFilename.empty() )
  {
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "Writing output transformation ... " << std::flush;
    
    vtkSmartPointer<vtkMatrix4x4> transformMatrix = transform->GetMatrix();
    if( !inputTransformationFilename.empty())
    {
      vtkMatrix4x4::Multiply4x4( transformMatrix, inputTransformationMatrix, transformMatrix );
    }
    
    ofstream txtFile;
    txtFile.open( outputTransformationFilename );
    
    for( unsigned int i = 0; i < 4; i++ )
    {
      for( unsigned int j = 0; j < 4; j++ )
      {
        txtFile << transformMatrix->GetElement( i, j );
        if( j < 3 )  txtFile << " ";
        if( j == 3 ) txtFile << "\n";
      }
    }
    
    txtFile.close();
    std::cout << "OK." << std::endl;
  }
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "icnsFlipVTKPolyData::Finished!"           << std::endl;
  std::cout << "========================================" << std::endl;
  
  return EXIT_SUCCESS;
}


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
