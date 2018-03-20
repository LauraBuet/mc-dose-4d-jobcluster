/** \file icnsEvaluateVTKPolyDataLabelAreas.cxx
 *
 *  \b Initial \b Author: Rene Werner \n\n
 *  \b Copyright (C) 2018 Department of Computational Neuroscience,
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
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkSmartPointer.h>

#include <vtkCell.h>
#include <vtkTriangle.h>

#include <vtkPointData.h>
#include <vtkPointLocator.h>

// Project includes: NONE SO FAR


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl; 
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsEvaluateVTKPolyDataLabelAreas -I <mesh> [...]" << std::endl;
  std::cout << std::endl;
  std::cout << "INPUT DATA" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "-I <source mesh>                    IN: Filename of mesh to analyze." << std::endl;
  std::cout << std::endl;
  std::cout << "OUTPUT DATA" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "-L <log file>                       OUT: NYI." << std::endl;
  std::cout << std::endl;
  std::cout << "-h                                  Print this help." << std::endl;
  std::cout << std::endl;
}

// ---------------------------------------------------------------
// Additional routine declarations:
// ---------------------------------------------------------------

//vtkMatrix4x4* Load4x4Matrix( std::string filename );

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
  std::cout << "icnsEvaluateVTKPolyDataLabelAreas         " << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  
  // -------------------------------------------------------------
  // Initializing parameters with default values:
  
  char* inputMeshFilename = NULL;
  char* outputLogFilename = NULL;
  
  // Reading parameters:
  
  int c;
  while( (c = getopt( argc, argv, "I:L:h?" )) != -1 )
  {
    switch( c )
    {
      case 'I':
        inputMeshFilename = optarg;
        std::cout << "Mesh (IN):                " << inputMeshFilename << std::endl;
        break;
      case 'L':
        outputLogFilename = optarg;
        std::cout << "Logfile (OUT):            " << outputLogFilename << std::endl;
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
  
  if( inputMeshFilename == NULL )
  {
    std::cerr << "ERROR: No source mesh specified!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input mesh data:
  
  vtkSmartPointer<vtkPolyData> inputMesh;
  vtkSmartPointer<vtkPolyDataReader> meshReader;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Loading input mesh ...                   " << std::flush;
  
  std::string filename = inputMeshFilename;
  
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "VTK data ... " << std::flush;
    meshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    meshReader->SetFileName( filename.c_str() );
    meshReader->Update();
    inputMesh = meshReader->GetOutput();
  }
  else
  {
    std::cout << "Mesh file ending not supported." << std::flush;
  }
  std::cout << "OK." << std::endl;
  
  // -------------------------------------------------------------
  // Actual work: Analyzing label areas:
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Analyzing label areas ...              " << std::endl;
  
  // Build a point locator that is subsequently used:
  vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
  pointLocator->SetDataSet( inputMesh );
  pointLocator->BuildLocator();
  
  std::cout << "In total " << inputMesh->GetNumberOfPoints() << " points." << std::endl;
  std::cout << "In total " << inputMesh->GetNumberOfCells() << " triangles." << std::endl;
  std::cout << "In total " << inputMesh->GetNumberOfPolys() << " polys." << std::endl;
  
  std::vector<unsigned int> pointsPerID(58,0);
  std::vector<float> areaPerID(58,0);
  
  // Write all of the coordinates of the points in the vtkPolyData to the console.
  for(vtkIdType iPointID = 0; iPointID < inputMesh->GetNumberOfPoints(); iPointID++)
  {
    double pointCoord[3];
    inputMesh->GetPoint( iPointID, pointCoord );
    //std::cout << "Point " << iPointID << " : (" << pointCoord[0] << " " << pointCoord[1] << " " << pointCoord[2] << ")" << std::endl;
    
    unsigned int pointScalarValue = inputMesh->GetPointData()->GetArray(0)->GetTuple1( iPointID );
    pointsPerID[pointScalarValue]++;
  }
  
  // Evalutate all triangles:
  for(vtkIdType iCellID = 0; iCellID < inputMesh->GetNumberOfCells(); iCellID++)
  {
    vtkCell* iCell = inputMesh->GetCell( iCellID );
    vtkTriangle* iTriangle = dynamic_cast< vtkTriangle* >(iCell);
    
    // coordinates of vertices of triangle:
    double p0[3];
    double p1[3];
    double p2[3];
    
    iTriangle->GetPoints()->GetPoint( 0, p0 );
    iTriangle->GetPoints()->GetPoint( 1, p1 );
    iTriangle->GetPoints()->GetPoint( 2, p2 );
    
    double iArea = vtkTriangle::TriangleArea(p0, p1, p2);
    
    vtkIdType iPointID = pointLocator->FindClosestPoint( p0 );
    double pointScalarValue = inputMesh->GetPointData()->GetArray(0)->GetTuple1( iPointID );
    areaPerID[pointScalarValue] = areaPerID[pointScalarValue] + iArea/3.0;
    
    iPointID = pointLocator->FindClosestPoint( p1 );
    pointScalarValue = inputMesh->GetPointData()->GetArray(0)->GetTuple1( iPointID );
    areaPerID[pointScalarValue] = areaPerID[pointScalarValue] + iArea/3.0;
    
    iPointID = pointLocator->FindClosestPoint( p2 );
    pointScalarValue = inputMesh->GetPointData()->GetArray(0)->GetTuple1( iPointID );
    areaPerID[pointScalarValue] = areaPerID[pointScalarValue] + iArea/3.0;
  }
  
  std::cout << "Distribution of scalars:" << std::endl;
  
  for( unsigned int i=0; i < pointsPerID.size(); i++ )
  {
    std::cout << "  LABEL " << i << ": " << pointsPerID[i] << " points; " << areaPerID[i] << " normalized area." << std::endl;
  }
  
  unsigned int nScalarArrays = inputMesh->GetPointData()->GetNumberOfArrays();
  if( nScalarArrays > 0 )
  {
    std::cout << "Point data associated with " << nScalarArrays << " arrays." << std::endl;
    std::cout << "Scalar range: [" << inputMesh->GetScalarRange()[0] << "; " <<
                                      inputMesh->GetScalarRange()[1] << "] " << std::endl;
    std::cout << "Number of entries: " << inputMesh->GetPointData()->GetArray(0)->GetNumberOfTuples() << std::endl;
  }
  else
  {
    std::cout << "No scalar data available." << std::endl;
  }
  
  // -------------------------------------------------------------
  // Writing output data:
  
  // NYI
  //std::cout << "------------------------------------------" << std::endl;
  //std::cout << "Writing log file ... "         << std::flush;
  //std::cout << "OK." << std::endl;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "icnsEvaluateVTKPolyDataLabelAreas FINISHED   " << std::endl;
  std::cout << "==========================================" << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


// ---------------------------------------------------------------
// Additional routines definitions: None so far
// ---------------------------------------------------------------

