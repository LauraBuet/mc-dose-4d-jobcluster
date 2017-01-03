/** \file icnsICPBasedSurfaceRegistration.cxx
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
#include <vtkIterativeClosestPointTransform.h>
#include <vtkLandmarkTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkMatrixToLinearTransform.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkTransformPolyDataFilter.h>

// Project includes:
#include "vtkTexturingHelper.h" // required to load multitexture OBJs


// ---------------------------------------------------------------
// Print help routine:
// ---------------------------------------------------------------

void PrintHelp()
{
  std::cout << std::endl;
  std::cout << "Usage:" << std::endl;
  std::cout << "icnsICPBasedSurfaceRegistration -F <fixed vtk mesh> -M <moving vtk mesh> [...]" << std::endl;
  std::cout << std::endl;
  std::cout << "INPUT DATA" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "-F <fixed vtk mesh>                 IN: Filename of fixed mesh." << std::endl;
  std::cout << "-M <moving vtk mesh>                IN: Filename of moving mesh." << std::endl;;
  std::cout << "-T <initial transformation>         IN: Filename of initial transformation." << std::endl;
  std::cout << std::endl;
  std::cout << "OUTPUT DATA" << std::endl;
  std::cout << "-----------------------------------------------------" << std::endl;
  std::cout << "-O <transformed moving vtk mesh>    OUT: Filename of transformed vtk mesh." << std::endl;
  std::cout << "-W <output transformation>          OUT: Filename of transformations."  << std::endl;
  std::cout << "                                         NB: If initial transformation is given, the output" << std::endl;
  std::cout << "                                         will be a concat of the initial and the actual tranform." << std::endl ;
  std::cout << std::endl;
  std::cout << "-n                                  Number of iterations (default: 100)" << std::endl;
  std::cout << "-h                                  Print this help." << std::endl;
  std::cout << "\n";
}

// ---------------------------------------------------------------
// Additional routine declarations:
// ---------------------------------------------------------------

void Write4x4Matrix( vtkMatrix4x4* matrix, std::string filename );
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
  std::cout << "icnsICPBasedSurfaceRegistration"            << std::endl;
  std::cout << "------------------------------------------" << std::endl;
  
  // -------------------------------------------------------------
  // Initializing parameters with default values:
  
  char* inputFixedMeshFilename             = NULL;
  char* inputMovingMeshFilename            = NULL;
  char* inputInitialTransformationFilename = NULL;
  
  char* outputWarpedMeshFilename           = NULL;
  char* outputTransformationFilename       = NULL;
  
  unsigned int nIterations                 = 100;
  
  // Reading parameters:
  
  int c;
  while( (c = getopt( argc, argv, "F:M:T:O:W:n:h?" )) != -1 )
  {
    switch( c )
    {
      case 'F':
        inputFixedMeshFilename = optarg;
        std::cout << "Fixed mesh (IN):               " << inputFixedMeshFilename << std::endl;
        break;
      case 'M':
        inputMovingMeshFilename = optarg;
        std::cout << "Moving mesh (IN):              " << inputMovingMeshFilename << std::endl;
        break;
      case 'T':
        inputInitialTransformationFilename = optarg;
        std::cout << "Initial transformation (IN):   " << inputInitialTransformationFilename << std::endl;
        break;
      case 'O':
        outputWarpedMeshFilename = optarg;
        std::cout << "Warped moving mesh (OUT):      " << outputWarpedMeshFilename << std::endl;
        break;
      case 'W':
        outputTransformationFilename = optarg;
        std::cout << "Computed transformation (OUT): " << outputTransformationFilename << std::endl;
        break;
      case 'n':
        nIterations = atoi( optarg );
        std::cout << "Number of ICP iterations:      " << nIterations << std::endl;
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
  
  if( inputFixedMeshFilename == NULL )
  {
    std::cerr << "ERROR: No fixed mesh specified!" << std::endl;
    return EXIT_FAILURE;
  }
  if( inputMovingMeshFilename == NULL )
  {
    std::cerr << "ERROR: No moving mesh specified!" << std::endl;
    return EXIT_FAILURE;
  }
  if( ( outputWarpedMeshFilename == NULL ) &&
      ( outputTransformationFilename = NULL ) )
  {
    std::cerr << "ERROR. No output file specified!" << std::endl;
    return EXIT_FAILURE;
  }
  
  // -------------------------------------------------------------
  // Loading input data:
  // First: fixed mesh.
  
  vtkSmartPointer<vtkPolyData> fixedMesh;
  vtkSmartPointer<vtkPolyData> movingMesh;
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Loading fixed mesh ... " << std::flush;
  
  vtkSmartPointer<vtkPolyDataReader> fixedMeshReader;
  vtkTexturingHelper fixedObjReader;
  std::string filename = inputFixedMeshFilename;
  
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "VTK data ... " << std::flush;
    fixedMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    fixedMeshReader->SetFileName( filename.c_str() );
    fixedMeshReader->Update();
    fixedMesh = fixedMeshReader->GetOutput();
  }
  else
  {
    std::cout << "OBJ data ... " << std::flush;
    fixedObjReader.ReadGeometryFile( filename.c_str() );
    fixedMesh = fixedObjReader.GetPolyData();
  }
  std::cout << "OK." << std::endl;
  
  // Second: moving mesh.
  
  std::cout << "Loading moving mesh ... " << std::flush;
  
  vtkSmartPointer<vtkPolyDataReader> movingMeshReader;
  vtkTexturingHelper movingObjReader;
  
  filename = inputMovingMeshFilename;
  if( filename.substr( filename.find_last_of(".") + 1) == "vtk" )
  {
    std::cout << "VTK data ... " << std::flush;
    movingMeshReader = vtkSmartPointer<vtkPolyDataReader>::New();
    movingMeshReader->SetFileName( filename.c_str() );
    movingMeshReader->Update();
    movingMesh = movingMeshReader->GetOutput();
  }
  else
  {
    std::cout << "OBJ data ... " << std::flush;
    movingObjReader.ReadGeometryFile( filename.c_str() );
    movingMesh = movingObjReader.GetPolyData();
  }
  std::cout << "OK." << std::endl;
  
  // Third: If specified, also load initial transformation.
  
  vtkSmartPointer<vtkMatrix4x4> initialTransformationMatrix;
  if( inputInitialTransformationFilename != NULL )
  {
    std::cout << "Loading moving mesh ... " << std::flush;
    initialTransformationMatrix = Load4x4Matrix( inputInitialTransformationFilename );
    std::cout << "OK." << std::endl;
    
    vtkSmartPointer<vtkMatrixToLinearTransform> initialMovingMeshTransform = vtkSmartPointer<vtkMatrixToLinearTransform>::New();
    initialMovingMeshTransform->SetInput( initialTransformationMatrix );
    
    vtkSmartPointer<vtkTransformPolyDataFilter> meshTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    meshTransformFilter->SetInputData( movingMesh );
    meshTransformFilter->SetTransform( initialMovingMeshTransform );
    meshTransformFilter->Update();
    
    movingMesh = meshTransformFilter->GetOutput();
  }
  
  // -------------------------------------------------------------
  // Performing ICP transformation:
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Preparing ICP filter ... "                  << std::endl;
  
  vtkSmartPointer<vtkIterativeClosestPointTransform> icp = vtkSmartPointer<vtkIterativeClosestPointTransform>::New();
  
  //std::cout << "-> FIXED MESH:" << std::endl;
  //fixedMesh->Print( std::cout );
  icp->SetSource( fixedMesh );
  
  //std::cout << "-> TARGET MESH:" << std::endl;
  //movingMesh->Print( std::cout );
  icp->SetTarget( movingMesh );
  
  //icp->DebugOn();
  icp->StartByMatchingCentroidsOn(); // WAS COMMENTED
  
  icp->SetMaximumNumberOfIterations( nIterations );
  //icp->SetMaximumNumberOfLandmarks(target->GetNumberOfPoints());
  icp->SetCheckMeanDistance(1);
  icp->SetMaximumMeanDistance(0.0000001);
  icp->GetLandmarkTransform()->SetModeToRigidBody();
  //icp->GetLandmarkTransform()->SetModeToAffine(); // WAS RIGID
  icp->Modified();
  icp->Update();
  icp->Inverse();
  
  // Get the resulting transformation matrix
  // (this matrix takes the target points to the source = fixed points):
  vtkSmartPointer<vtkMatrix4x4> transformMatrix = icp->GetMatrix();
  std::cout << "The resulting transform matrix is: " << *transformMatrix << std::endl;

  // -------------------------------------------------------------
  // Writing output data:
  
  std::cout << "------------------------------------------" << std::endl;
  std::cout << "Writing output image:"                      << std::flush;

  // Writing transformation.
  // Case 1: no input transform given.
  
  if( (outputTransformationFilename != NULL) && (inputInitialTransformationFilename == NULL) )
  {
    std::cout << "-> Saving transformation ... " << std::flush;
    
    Write4x4Matrix( transformMatrix, outputTransformationFilename );
    
    std::cout << "OK." << std::endl;
    std::cout << "[NB: transform is defined to warp MOVING MESH points.]" << std::endl;
  }
  
  // Case 2: input transform given.
  // Thus: Concatenation of input transform and ICP result required.
  
  if( (outputTransformationFilename != NULL) && (inputInitialTransformationFilename != NULL) )
  {
    std::cout << "-> Saving transformation ... " << std::endl;
    std::cout << "   Preloading input transform: " << inputInitialTransformationFilename << std::endl;
    
    vtkMatrix4x4::Multiply4x4( transformMatrix, initialTransformationMatrix, transformMatrix );
    Write4x4Matrix( transformMatrix, outputTransformationFilename );
    
    std::cout << "OK." << std::endl;
    std::cout << "[NB: transform is defined to warp MOVING MESH points.]" << std::endl;
  }
 
  // Writing transformed moving mesh:
  
  if( outputWarpedMeshFilename != NULL )
  {
    std::cout << "-> Saving warped moving mesh ... " << std::flush;
    
    // Transforming moving polydata:
    
    vtkSmartPointer<vtkTransformPolyDataFilter> icpTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    #if VTK_MAJOR_VERSION <= 5
    icpTransformFilter->SetInput( movingMesh );
    #else
    icpTransformFilter->SetInputData( movingMesh );
    #endif
    icpTransformFilter->SetTransform( icp );
    icpTransformFilter->Update();
    
    // Writing polydata:
    
    vtkSmartPointer<vtkPolyDataWriter> warpedMeshWriter = vtkSmartPointer<vtkPolyDataWriter>::New();
    warpedMeshWriter->SetFileName( outputWarpedMeshFilename );
    warpedMeshWriter->SetInputData( icpTransformFilter->GetOutput() );
    warpedMeshWriter->Write();
    
    std::cout << "OK." << std::endl;
  }
  
  std::cout << "----------------------------------------" << std::endl;
  std::cout << "icnsICPBasedSurfaceRegistration FINISHED" << std::endl;
  std::cout << "========================================" << std::endl;
  std::cout << std::endl;

  return EXIT_SUCCESS;
}


// ---------------------------------------------------------------
// Additional routines definitions:
// ---------------------------------------------------------------

void Write4x4Matrix( vtkMatrix4x4* matrix, std::string filename )
{
  ofstream txtFile;
  txtFile.open( filename );
  
  for( unsigned int i = 0; i < 4; i++ )
  {
    for( unsigned int j = 0; j < 4; j++ )
    {
      txtFile << matrix->GetElement( i, j );
      if( j < 3 )  txtFile << " ";
      if( j == 3 ) txtFile << "\n";
    }
  }
  
  txtFile.close();
} // end of Write4x4Matrix()

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
