/** \file icnsTransferVTKPolyDataScalars.cxx
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
#include <vtkMatrix4x4.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkTransformPolyDataFilter.h>

#include <vtkPointData.h>
#include <vtkPointLocator.h>

// Project includes:
#include "vtkTexturingHelper.h" // required to load multitexture OBJs


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl; 
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsTransferVTKPolyDataScalars -S <source mesh> -T <target mesh> [...]" << std::endl;
  std::cout << std::endl;
  std::cout << "INPUT DATA" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "-S <source mesh>                    IN: Filename of source mesh (scalars to be transfered)." << std::endl;
  std::cout << "-T <target mesh>                    IN: Filename of target mesh (scalars to be defined)." << std::endl;;
  std::cout << "-W <source mesh transformation>     IN: Filename of transformation of source mesh." << std::endl;
  std::cout << std::endl;
  std::cout << "OUTPUT DATA" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "-O <recolored target mesh>          OUT: Filename of recolored vtk mesh." << std::endl;
  std::cout << std::endl;
  std::cout << "-h                                  Print this help." << std::endl;
  std::cout << std::endl;
}

// ---------------------------------------------------------------
// Additional routine declarations:
// ---------------------------------------------------------------

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
  std::cout << "==========================================" << std::endl;
  std::cout << "icnsTransferVTKPolyDataScalars            " << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  
  // -------------------------------------------------------------
  // Initializing parameters with default values:
  
  char* inputSourceMeshFilename            = NULL;
  char* inputSourceTransformationFilename  = NULL;
  char* inputTargetMeshFilename            = NULL;
  
  char* outputRecoloredTargetMeshFilename  = NULL;
  
  // Reading parameters:
  
  int c;
  while( (c = getopt( argc, argv, "S:T:W:O:h?" )) != -1 )
  {
    switch( c )
    {
      case 'S':
        inputSourceMeshFilename = optarg;
        std::cout << "Source mesh (IN):                " << inputSourceMeshFilename << std::endl;
        break;
      case 'T':
        inputTargetMeshFilename = optarg;
        std::cout << "Moving mesh (IN):                " << inputTargetMeshFilename << std::endl;
        break;
      case 'W':
        inputSourceTransformationFilename = optarg;
        std::cout << "Source mesh transformation (IN): " << inputSourceTransformationFilename << std::endl;
        break;
      case 'O':
        outputRecoloredTargetMeshFilename = optarg;
        std::cout << "Recolored target mesh (OUT):     " << outputRecoloredTargetMeshFilename << std::endl;
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
  
  if( inputSourceMeshFilename == NULL )
  {
    std::cerr << "ERROR: No source mesh specified!" << std::endl;
    return EXIT_FAILURE;
  }
  if( inputTargetMeshFilename == NULL )
  {
    std::cerr << "ERROR: No target mesh specified!" << std::endl;
    return EXIT_FAILURE;
  }
  if( outputRecoloredTargetMeshFilename == NULL )
  {
    std::cerr << "ERROR. No output file specified!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input data:
  // First: source mesh.
  
  vtkSmartPointer<vtkPolyData> sourceMesh;
  vtkSmartPointer<vtkPolyData> targetMesh;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Loading source mesh ...                   " << std::flush;
  
  vtkSmartPointer<vtkPolyDataReader> sourceMeshReader;
  vtkTexturingHelper sourceObjReader;
  std::string filename = inputSourceMeshFilename;
  
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "VTK data ... " << std::flush;
    sourceMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    sourceMeshReader->SetFileName( filename.c_str() );
    sourceMeshReader->Update();
    sourceMesh = sourceMeshReader->GetOutput();
  }
  else
  {
    std::cout << "OBJ data ... " << std::flush;
    sourceObjReader.ReadGeometryFile( filename.c_str() );
    sourceMesh = sourceObjReader.GetPolyData();
  }
  std::cout << "OK." << std::endl;
  
  // Second: target mesh.
  
  std::cout << "Loading target mesh ... " << std::flush;
  
  vtkSmartPointer<vtkPolyDataReader> targetMeshReader;
  vtkTexturingHelper targetObjReader;
  
  filename = inputTargetMeshFilename;
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "VTK data ... " << std::flush;
    targetMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    targetMeshReader->SetFileName( filename.c_str() );
    targetMeshReader->Update();
    targetMesh = targetMeshReader->GetOutput();
  }
  else
  {
    std::cout << "OBJ data ... " << std::flush;
    targetObjReader.ReadGeometryFile( filename.c_str() );
    targetMesh = targetObjReader.GetPolyData();
  }
  std::cout << "OK." << std::endl;
  
  // Third: If specified, also load initial transformation.
  
  vtkSmartPointer<vtkMatrix4x4> initialTransformationMatrix;
  if( inputSourceTransformationFilename != NULL )
  {
    std::cout << "Loading source mesh transformation ... " << std::flush;
    initialTransformationMatrix = Load4x4Matrix( inputSourceTransformationFilename );
    std::cout << "OK." << std::endl;
    
    vtkSmartPointer<vtkMatrixToLinearTransform> initialSourceMeshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
    initialSourceMeshTransform->SetInput( initialTransformationMatrix );
    
    vtkSmartPointer<vtkTransformPolyDataFilter> meshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    meshTransformFilter->SetInputData( sourceMesh );
    meshTransformFilter->SetTransform( initialSourceMeshTransform );
    meshTransformFilter->Update();
    
    sourceMesh = meshTransformFilter->GetOutput();
  }
  
  // -------------------------------------------------------------
  // Actual work: Transferring scalars
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Transferring scalar data ...              " << std::endl;
  
  vtkSmartPointer<vtkPolyData> outputMesh = vtkSmartPointer<vtkPolyData>::New();
  outputMesh->DeepCopy( targetMesh );
  
  // Build a point locator that is subsequently used:
  vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
  pointLocator->SetDataSet( sourceMesh );
  pointLocator->BuildLocator();
  
  // Iterate over points of target and output mesh, respectively:
  unsigned int nOutputPoints = outputMesh->GetNumberOfPoints();
  for( unsigned int iOutputMeshPointIDs = 0; iOutputMeshPointIDs < nOutputPoints; iOutputMeshPointIDs++ )
  {
    // Get point location of current point:
    double currentOutputMeshPointCoord[3];
    outputMesh->GetPoint( iOutputMeshPointIDs, currentOutputMeshPointCoord );
    
    // Get closest point in source mesh:
    vtkIdType currentPointID = pointLocator->FindClosestPoint( currentOutputMeshPointCoord );
    double currentPointScalarValue = sourceMesh->GetPointData()->GetArray(0)->GetTuple1( currentPointID );
    
    // Replace scalar data of current point in output mesh by scalar value of
    // closest point in source mesh:
    
    outputMesh->GetPointData()->GetArray(0)->SetTuple1( iOutputMeshPointIDs, currentPointScalarValue );
  }
  
  /*
  
  std::vector<unsigned int> pointsPerID(58,0);
  
  // Write all of the coordinates of the points in the vtkPolyData to the console.
  for(vtkIdType iPointID = 0; iPointID < outputMesh->GetNumberOfPoints(); iPointID++)
  {
    double pointCoord[3];
    outputMesh->GetPoint( iPointID, pointCoord );
    //std::cout << "Point " << iPointID << " : (" << pointCoord[0] << " " << pointCoord[1] << " " << pointCoord[2] << ")" << std::endl;
    
    unsigned int pointScalarValue = outputMesh->GetPointData()->GetArray(0)->GetTuple1( iPointID );
    pointsPerID[pointScalarValue]++;
  }
  
  std::cout << "In total " << outputMesh->GetNumberOfPoints() << " points!" << std::endl;
  std::cout << "Distribution of scalars:" << std::endl;
  
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
  }
  
  
  outputMesh->DeepCopy( targetMesh );
   */

  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Writing recolored target mesh ... "         << std::flush;
  
  vtkSmartPointer<vtkPolyDataWriter> outputMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
  outputMeshWriter->SetFileName( outputRecoloredTargetMeshFilename );
  outputMeshWriter->SetInputData( outputMesh );
  outputMeshWriter->Write();
    
  std::cout << "OK." << std::endl;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "icnsTransferVTKPolyDataScalars FINISHED   " << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


// ---------------------------------------------------------------
// Additional routines definitions:
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
  
  std::cout << *mat << std::endl;

  return mat;
}
