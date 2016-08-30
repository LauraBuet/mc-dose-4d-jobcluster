/** \file imiSurfaceConstrainedDiffusionRegistrationFilter.cxx
 *
 *  \b Initial \b Author: Jan Ehrhardt\n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#include <vtkPointLocator.h>
#include <vtkCellLocator.h>
#include <vtkMath.h>
#include <vtkTriangle.h>
#include <vtkTriangleFilter.h>
#include <vtkCellArray.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkSmartPointer.h>

#include "icnsSurfaceConstrainedDiffusionRegistrationFilter.h"

using namespace imi;

///////////////////////////////////////////////////////////////////////////////////////////
//
// imiSurfaceConstrainedDiffusionRegistrationFilter::New()
//
///////////////////////////////////////////////////////////////////////////////////////////
imiSurfaceConstrainedDiffusionRegistrationFilter* imiSurfaceConstrainedDiffusionRegistrationFilter
::New()
{
  return new imiSurfaceConstrainedDiffusionRegistrationFilter();
}

imiSurfaceConstrainedDiffusionRegistrationFilter
::imiSurfaceConstrainedDiffusionRegistrationFilter()
{
  m_bPreComputeForceField = true;
  m_SmoothingSigma = 2.0;

  m_spTargetPointLocator = vtkSmartPointer<vtkPointLocator>::New();
  m_spReferencePointLocator = vtkSmartPointer<vtkPointLocator>::New();

  m_spDeformedTargetSurface = NULL;
  m_spInitialTransform = NULL;

  m_XForce = NULL;
  m_YForce = NULL;
  m_ZForce = NULL;
  m_XForceSmooth = NULL;
  m_YForceSmooth = NULL;
  m_ZForceSmooth = NULL;
  m_LastXForce = NULL;
  m_LastYForce = NULL;
  m_LastZForce = NULL;
  m_NumberOfForceVectors = 0;

  m_bUsePrecomputedNeighborhoodLists = false;
}

imiSurfaceConstrainedDiffusionRegistrationFilter
::~imiSurfaceConstrainedDiffusionRegistrationFilter()
{
  if (!ClearForces())
  {
    std::cerr << "Can not free forces memory!" << std::endl;
  }
}

bool imiSurfaceConstrainedDiffusionRegistrationFilter
::ClearForces()
{
  if( m_XForce != NULL )
    delete m_XForce;
  if( m_YForce != NULL )
    delete m_YForce;
  if( m_ZForce != NULL )
    delete m_ZForce;
  if( m_XForceSmooth != NULL )
    delete m_XForceSmooth;
  if( m_YForceSmooth != NULL )
    delete m_YForceSmooth;
  if( m_ZForceSmooth != NULL )
    delete m_ZForceSmooth;
  if( m_LastXForce != NULL )
    delete m_LastXForce;
  if( m_LastYForce != NULL )
    delete m_LastYForce;
  if( m_LastZForce != NULL )
    delete m_LastZForce;
  m_XForce = NULL;
  m_YForce = NULL;
  m_ZForce = NULL;
  m_XForceSmooth = NULL;
  m_YForceSmooth = NULL;
  m_ZForceSmooth = NULL;
  m_LastXForce = NULL;
  m_LastYForce = NULL;
  m_LastZForce = NULL;

  m_NumberOfForceVectors = 0;

  return true;
}

vtkSmartPointer<vtkPolyData> imiSurfaceConstrainedDiffusionRegistrationFilter
::GetTransformedTargetSurface()
{
  if( GetTargetSurface() == NULL )
    return NULL;

  if( m_spDeformedTargetSurface == NULL )
  {
    m_spDeformedTargetSurface = vtkSmartPointer<vtkPolyData>::New();
    m_spDeformedTargetSurface->DeepCopy( GetTargetSurface() );
  }
  return m_spDeformedTargetSurface;
}

///////////////////////////////////////////////////
//
//  InitializeRegistration()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::InitializeRegistration()
{
  if( m_bIsInitialized )
  {
    return true;
  }

  if( !imiSurfaceRegistrationFilter::InitializeRegistration() )
  {
    return false;
  }

  std::cout << " imiSurfaceConstrainedDiffusionRegistrationFilter::InitializeRegistration()..." << std::endl;

  if( m_spDeformedTargetSurface == NULL )
  {
    vtkSmartPointer<vtkPolyData> dummy = GetTransformedTargetSurface();
    if( dummy == NULL )
    {
      std::cerr <<  "Can't create deformed Source!" << std::endl;
      return false;
    }
  }

  std::cout << "  Parameters: Sigma=" << GetSmoothingSigma() << "  Iterations=" << GetNumberOfIterations() << std::endl;

  std::cout << "  Building point locators..." << std::endl;

  // build the locator
  m_spTargetPointLocator->SetDataSet( GetTargetSurface() );
  m_spTargetPointLocator->Initialize();
  m_spTargetPointLocator->BuildLocator();

  m_spReferencePointLocator->SetDataSet( GetReferenceSurface() );
  m_spReferencePointLocator->Initialize();
  m_spReferencePointLocator->BuildLocator();

  // delete memory for the forces
  if( !ClearForces() )
  {
    std::cerr <<  "Cannot free forces memory" << std::endl;
    return false;
  }

  // delete neighbourhood-lists
  m_TargetNeighborIdListVector.clear();
  m_TargetNeighborWeightListVector.clear();

  std::cout << "  Allocating force array..." << std::endl;
  // allocate space for forces
  m_NumberOfForceVectors = GetTargetSurface()->GetNumberOfPoints();
  m_XForce = new double[m_NumberOfForceVectors];
  m_YForce = new double[m_NumberOfForceVectors];
  m_ZForce = new double[m_NumberOfForceVectors];
  m_XForceSmooth = new double[m_NumberOfForceVectors];
  m_YForceSmooth = new double[m_NumberOfForceVectors];
  m_ZForceSmooth = new double[m_NumberOfForceVectors];

  // initialize force field with zero
  std::cout << "   Setting initial force field to zero..." << std::endl;

  for( int i = 0; i < m_NumberOfForceVectors; i++ )
  {
    m_XForce[i] = 0;
    m_YForce[i] = 0;
    m_ZForce[i] = 0;
    m_XForceSmooth[i] = 0;
    m_YForceSmooth[i] = 0;
    m_ZForceSmooth[i] = 0;
  }

  if( m_spInitialTransform.GetPointer() != NULL )
  {
    if( !InitializeDeformedTargetSurfaceByInitialTransform() )
    {
      std::cerr <<  "InitializeDeformedTargetSurfaceByInitialTransform() failed!" << std::endl;
      return false;
    }
  }

  m_bIsInitialized = true;

  return true;
}
///////////////////////////////////////////////////
//
//  UpdateRegistration()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::UpdateRegistration()
{
  std::cout << " Update SurfaceConstrainedDiffusion ..." << std::endl;

  /*
   if( CheckGlobalDebugLevel( 9 ) )
  {
    char filename[200];
    sprintf( filename, "%sIteration_%.3d_DeformedTargetSurface.vtk",
        m_sDebugPathName.c_str(), m_iCurrentIteration + 1 );
    imiDEBUGINFO( 9, "    Writing deformed target surface to file " << filename );

    WriteSurface( GetTransformedTargetSurface(), filename );
  }
  */

  // Use initial force field in first iteration
  if( !ComputeForceField() )
  {
    std::cerr <<  "Computing force field failed in iteration " << m_iCurrentIteration << "!" << std::endl;
    return false;
  }

  if( !SmoothForceField() )
  {
    std::cerr <<  "Smoothing force field failed in iteration " << m_iCurrentIteration << "!" << std::endl;
    return false;
  }

  if( !ApplyForceField() )
  {
    std::cerr <<  "Applying force field failed in iteration " << m_iCurrentIteration << "!" << std::endl;
    return false;
  }

  // Compute the mean change in the force vectors.
  m_dLastRMSValue = ComputeRMSValue();

  return true;
}

///////////////////////////////////////////////////
//
//  ApplyForceField()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::ApplyForceField()
{
  return ApplyForceFieldSurface(
      GetTargetSurface(),
      GetTransformedTargetSurface() );
}

///////////////////////////////////////////////////
//
//  ApplyForceFieldSurface()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::ApplyForceFieldSurface(
    vtkSmartPointer<vtkPolyData> inputSurface,
    vtkSmartPointer<vtkPolyData> deformedInputSurface)
{
  std::cout << " ApplyForceFieldSurface..." << std::endl;

  if( inputSurface->GetNumberOfPoints() != deformedInputSurface->GetNumberOfPoints() )
  {
    std::cerr << 
        "imiSurfaceConstrainedDiffusionRegistrationFilter::ApplyForceFieldSurface() "
        "Input and deformed have different number of points!" << std::endl;
    return false;
  }

  if( inputSurface->GetNumberOfPoints() != m_NumberOfForceVectors )
  {
    std::cerr << 
        "imiSurfaceConstrainedDiffusionRegistrationFilter::ApplyForceFieldSurface() "
        "Input and force array have different number of points!" << std::endl;
    return false;
  }

  double inputPoint[3];

  int numPoints = inputSurface->GetNumberOfPoints();
  double normSum = 0.0;

  for (int i = 0; i < numPoints; i++)
  {
    inputSurface->GetPoint( i, inputPoint );

    inputPoint[0] += m_XForceSmooth[i];
    inputPoint[1] += m_YForceSmooth[i];
    inputPoint[2] += m_ZForceSmooth[i];

    normSum += vcl_sqrt(
        vnl_math_sqr( m_XForceSmooth[i] ) +
        vnl_math_sqr( m_YForceSmooth[i] ) +
        vnl_math_sqr( m_ZForceSmooth[i] ) );

    deformedInputSurface->GetPoints()->SetPoint( i, inputPoint );
  }
  deformedInputSurface->Modified();

  std::cout << "  Norm: " << normSum / numPoints << std::endl;

  return true;
}

///////////////////////////////////////////////////
//
//  SmoothForceField()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::SmoothForceField()
{
  return SmoothForceFieldSurface();
}

///////////////////////////////////////////////////
//
//  SmoothForceFieldSurface()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::SmoothForceFieldSurface()
{
  std::cout << " SmoothForceFieldSurface (s=" << GetSmoothingSigma() << ") ... " << std::endl;

  // Parameters for debug output
  double meanGaussSum = 0.0;
  double meanFoundPts = 0.0;
  double normSum = 0.0;

  // Get number of surface points.
  vtkSmartPointer<vtkPolyData> targetSurface = GetTargetSurface();
  unsigned int numPoints = targetSurface->GetNumberOfPoints();

  //
  // Pre-compute neighborhood lists if requested.
  //
  if( m_bUsePrecomputedNeighborhoodLists )
  {
    //
    // Check if BuildNeighbourhoodIDList() has to be called.
    //
    if( GetCurrentIteration() == 0 || m_TargetNeighborIdListVector.size() != numPoints )
    {
      if( !BuildNeighbourhoodIDList(
              targetSurface,
              GetSmoothingSigma(),
              m_TargetNeighborIdListVector,
              m_TargetNeighborWeightListVector ) )
      {
        std::cerr <<  "SmoothForceFieldSurface(): Cannot build Neighborhood ID Lists!" << std::endl;
        return false;
      }
      // make a small check
      if( m_TargetNeighborIdListVector.size() != numPoints ||
          m_TargetNeighborWeightListVector.size() != numPoints )
      {
        std::cerr <<  "SmoothForceFieldSurface(): Building Neighborhood ID Lists failed!" << std::endl;
        return false;
      }
    }
  }

  // Define thread parameters.
  icnsImageThreads* imageThreads = new icnsImageThreads();
  unsigned int numberOfThreads = imageThreads->GetNumberOfThreads();


  // ==========================================
  // Calculate smoothed points.
  imiSmoothPointsParameters* forceParameters;
  unsigned int threadRange = numPoints / numberOfThreads;
  for( unsigned int threadNr = 0; threadNr < numberOfThreads; ++threadNr )
  {
    // Initialize thread parameters.
    forceParameters = new imiSmoothPointsParameters();

    // Calculate thread range.
    forceParameters->from = threadNr * threadRange;
    forceParameters->to = (threadNr == numberOfThreads - 1) ?
        numPoints : (threadNr + 1) * threadRange;

    // Initializing thread parameters.
    forceParameters->idListVector = m_TargetNeighborIdListVector;
    forceParameters->weightListVector = m_TargetNeighborWeightListVector;
    forceParameters->surface = targetSurface;
    forceParameters->pointLocator = m_spTargetPointLocator;
    forceParameters->sigma = this->GetSmoothingSigma();
    forceParameters->xForce = m_XForce;
    forceParameters->yForce = m_YForce;
    forceParameters->zForce = m_ZForce;
    forceParameters->xForceSmooth = m_XForceSmooth;
    forceParameters->yForceSmooth = m_YForceSmooth;
    forceParameters->zForceSmooth = m_ZForceSmooth;
    forceParameters->numPoints = numPoints;

    forceParameters->meanGaussSum = 0.0;
    forceParameters->meanFoundPts = 0.0;
    forceParameters->normSum = 0.0;

//    forceParameters->surface->Update();
//    forceParameters->pointLocator->BuildLocator();

    // Set thread method.
    if( m_bUsePrecomputedNeighborhoodLists )
    {
      imageThreads->SetThreadMethod( threadNr, &T_SmoothPointsPrecomputed, forceParameters );
    }
    else
    {
      imageThreads->SetThreadMethod( threadNr, &T_SmoothPoints, forceParameters );
    }
  }

  // Execute the threads.
  imageThreads->ExecuteThreads();

  // Clean up.
  for( unsigned int threadNr = 0; threadNr < numberOfThreads; ++threadNr )
  {
    // Get thread parameters.
    forceParameters = (imiSmoothPointsParameters*)
        imageThreads->GetThreadParameters( threadNr );

    // Get output of thread for debug information.
    meanGaussSum += forceParameters->meanGaussSum;
    meanFoundPts += forceParameters->meanFoundPts;
    normSum += forceParameters->normSum;

    // Delete thread parameters.
    delete imageThreads->GetThreadParameters( threadNr );
  }
  delete imageThreads;

  /*
  if( CheckGlobalDebugLevel( 7 ) )
  {
    meanGaussSum /= numPoints;
    meanFoundPts /= numPoints;

    imiDEBUGINFO( 7, "          Norm: " << normSum / numPoints );
    imiDEBUGINFO( 7, "  meanGaussSum: " << meanGaussSum << " meanPts found: " << meanFoundPts );
  }
  */

  return true;
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// T_SmoothPointsPrecomputed()
//
///////////////////////////////////////////////////////////////////////////////////////////
IMI_THREAD_RETURN_TYPE imiSurfaceConstrainedDiffusionRegistrationFilter
::T_SmoothPointsPrecomputed( void* vargs )
{
  imiSmoothPointsParameters* args = (imiSmoothPointsParameters*) vargs;

  std::vector<NeighborhoodIDListType>     idListVector     = args->idListVector;
  std::vector<NeighborhoodWeightListType> weightListVector = args->weightListVector;

  double *xForce = args->xForce;
  double *yForce = args->yForce;
  double *zForce = args->zForce;
  double *xForceSmooth = args->xForceSmooth;
  double *yForceSmooth = args->yForceSmooth;
  double *zForceSmooth = args->zForceSmooth;

  NeighborhoodIDListType     *idList     = NULL; // IDs of the neighboring points
  NeighborhoodWeightListType *weightList = NULL; // Gaussian weights of the neighboring points

  double newForce[3];

  for( unsigned int i = args->from; i < args->to; ++i )
  {
    idList = &idListVector[i];
    weightList = &weightListVector[i];

    // small check
    if( idList->size() != weightList->size() )
    {
      std::cerr <<  "SmoothForceFieldSurface(): Neighborhood ID List and weight list don't have the same size!" << std::endl;
      return false;
    }

    int foundPoints = idList->size();

    newForce[0] = 0;
    newForce[1] = 0;
    newForce[2] = 0;

    double gaussFac;
    double gaussSum = 0.0;

    // iterate through the list
    NeighborhoodIDListType::iterator idListIterator;
    NeighborhoodWeightListType::iterator weightListIterator;

    for( idListIterator = idList->begin(), weightListIterator = weightList->begin();
         idListIterator != idList->end();
         ++idListIterator, ++weightListIterator )
    {
      // get neighborhood surface point
      NeighborhoodIDType id = *idListIterator;

      // get gaussian weight for this neighborhood point
      gaussFac = static_cast<double>(*weightListIterator);
      gaussSum += gaussFac;

      // add all weighted forces
      newForce[0] += xForce[id] * gaussFac;
      newForce[1] += yForce[id] * gaussFac;
      newForce[2] += zForce[id] * gaussFac;
    }

    // normalize with the gaussSum to preserve length
    xForceSmooth[i] = newForce[0] / gaussSum;
    yForceSmooth[i] = newForce[1] / gaussSum;
    zForceSmooth[i] = newForce[2] / gaussSum;

    // output debug information
    /*
    if( CheckGlobalDebugLevel( 7 ) )
    {
      args->normSum += vcl_sqrt(
          xForceSmooth[i] * xForceSmooth[i] +
          yForceSmooth[i] * yForceSmooth[i] +
          zForceSmooth[i] * zForceSmooth[i] );

      args->meanGaussSum += gaussSum;
      args->meanFoundPts += foundPoints;
    }
    */
  }

  return IMI_THREAD_RETURN_VALUE;
}


///////////////////////////////////////////////////////////////////////////////////////////
//
// T_SmoothPoints()
//
///////////////////////////////////////////////////////////////////////////////////////////
IMI_THREAD_RETURN_TYPE imiSurfaceConstrainedDiffusionRegistrationFilter
::T_SmoothPoints( void* vargs )
{
  imiSmoothPointsParameters* args = (imiSmoothPointsParameters*) vargs;

  vtkPolyData*     surface      = args->surface;
  vtkPointLocator* pointLocator = args->pointLocator;

  double *xForce = args->xForce;
  double *yForce = args->yForce;
  double *zForce = args->zForce;
  double *xForceSmooth = args->xForceSmooth;
  double *yForceSmooth = args->yForceSmooth;
  double *zForceSmooth = args->zForceSmooth;

  const double sigma = args->sigma;
  const double sigma2 = vnl_math_sqr( sigma );
  const double radius = 2*sigma;
  const double radius2 = vnl_math_sqr( radius );

  double newForce[3];
  double surfacePoint1[3];
  double surfacePoint2[3];

  double dist2 = 0.0;
  double gaussFac = 0.0;
  double gaussSum = 0.0;

  // Define neighborhood ID list
  vtkIdList *idList = vtkIdList::New();

  for( unsigned int i = args->from; i < args->to; ++i )
  {
    surface->GetPoint( i, surfacePoint1 );

    idList->Reset();
    pointLocator->FindPointsWithinRadius(
        radius, surfacePoint1, idList );

    int foundPoints = idList->GetNumberOfIds();
    if( foundPoints < 1 )
    {
      std::cerr <<  "No Points in Radius " << radius << "found !" << std::endl;
      std::cerr <<  "Increase sigma!" << std::endl;
      return false;
    }

    newForce[0] = 0.0;
    newForce[1] = 0.0;
    newForce[2] = 0.0;

    gaussSum = 0.0;

    for( int j = 0; j < foundPoints; j++ )
    {

      int id = idList->GetId( j );

      if( id >= static_cast<int>(args->numPoints) || id < 0 )
      {
        std::cerr <<  "Invalid Id: " << id << std::endl;
        continue;
      }

      surface->GetPoint( id, surfacePoint2 );

      dist2 = vtkMath::Distance2BetweenPoints( surfacePoint1, surfacePoint2 );
      if( dist2 > radius2 )
      {
        std::cerr <<  "Distance > Radius: " << dist2 << std::endl;
        continue;
      }

      gaussFac = vcl_exp( -(dist2 / (2.0 * sigma2)) );
      gaussSum += gaussFac;

      newForce[0] += xForce[id] * gaussFac;
      newForce[1] += yForce[id] * gaussFac;
      newForce[2] += zForce[id] * gaussFac;
    }

    // normalize with the gaussSum to preserve length
    xForceSmooth[i] = newForce[0] / gaussSum;
    yForceSmooth[i] = newForce[1] / gaussSum;
    zForceSmooth[i] = newForce[2] / gaussSum;

    // output debug information
    /*
    if( CheckGlobalDebugLevel( 7 ) )
    {
      args->normSum += vcl_sqrt(
          vnl_math_sqr( xForceSmooth[i] ) +
          vnl_math_sqr( yForceSmooth[i] ) +
          vnl_math_sqr( zForceSmooth[i] ) );

      args->meanGaussSum += gaussSum;
      args->meanFoundPts += foundPoints;
    }
    */
  }

  // Delete the list
  idList->Delete();

  return IMI_THREAD_RETURN_VALUE;
}


///////////////////////////////////////////////////
//
//  ComputeForceField()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::ComputeForceField()
{
  std::cout << "  Computing force field..." << std::endl;

  // Define thread parameters.
  icnsImageThreads* imageThreads = new icnsImageThreads();
  unsigned int numberOfThreads = imageThreads->GetNumberOfThreads();

  // ==========================================
  // Calculate force field.
  int numPoints = this->GetTargetSurface()->GetNumberOfPoints();
  imiComputeForceFieldParameters* parameters;
  unsigned int threadRange = numPoints / numberOfThreads;

  for( unsigned int threadNr = 0; threadNr < numberOfThreads; ++threadNr )
  {
    // Initialize thread parameters.
    parameters = new imiComputeForceFieldParameters();

    // Calculate thread range.
    parameters->from = threadNr * threadRange;
    parameters->to = (threadNr == numberOfThreads - 1) ?
        numPoints : (threadNr + 1) * threadRange;

    // Initializing thread parameters.
    parameters->referenceSurface = this->GetReferenceSurface();
    parameters->targetSurface = this->GetTargetSurface();
    parameters->deformedTargetSurface = this->GetTransformedTargetSurface();
    parameters->referencePointLocator = m_spReferencePointLocator;
    parameters->xForce = m_XForce;
    parameters->yForce = m_YForce;
    parameters->zForce = m_ZForce;
    parameters->normSum = 0.0;

    // Set thread method.
    imageThreads->SetThreadMethod( threadNr, &T_ComputeForceField, parameters );
  }

  // Execute the threads.
  imageThreads->ExecuteThreads();

  // Clean up.
  double normSum = 0.0;
  for( unsigned int threadNr = 0; threadNr < numberOfThreads; ++threadNr )
  {
    // Get thread parameters.
    parameters = (imiComputeForceFieldParameters*)
        imageThreads->GetThreadParameters( threadNr );

    // Get output of thread for debug information.
    normSum += parameters->normSum;

    // Delete thread parameters.
    delete imageThreads->GetThreadParameters( threadNr );
  }
  delete imageThreads;

  std::cout << "  Norm: " << normSum / numPoints << std::endl;

  return true;
}


///////////////////////////////////////////////////
//
//  T_ComputeForceField()
//
///////////////////////////////////////////////////
IMI_THREAD_RETURN_TYPE imiSurfaceConstrainedDiffusionRegistrationFilter
::T_ComputeForceField( void *vargs )
{
  imiComputeForceFieldParameters* args = (imiComputeForceFieldParameters*) vargs;

  vtkPolyData* referenceSurface = args->referenceSurface;
  vtkPolyData* targetSurface = args->targetSurface;
  vtkPolyData* deformedTargetSurface = args->deformedTargetSurface;

  vtkPointLocator* referencePointLocator = args->referencePointLocator;

  double *xForce = args->xForce;
  double *yForce = args->yForce;
  double *zForce = args->zForce;

  int referencePointId;
  double referencePoint[3];
  double targetPoint[3];
  double deformedTargetPoint[3];

  for( unsigned int i = args->from; i < args->to; ++i )
  {
    targetSurface->GetPoint( i, targetPoint );
    deformedTargetSurface->GetPoint( i, deformedTargetPoint );

    // Find closest reference point to actual deformed target point.
    referencePointId = referencePointLocator->FindClosestPoint(
        deformedTargetPoint );
    if( referencePointId < 0 )
    {
      std::cerr <<  "No nearest neighbor found for target point ["
          << deformedTargetPoint[0] << ","
          << deformedTargetPoint[1] << ","
          << deformedTargetPoint[2] << "]!" << std::endl;
      continue;
    }

    // Compute deformation vector to move target point to reference point.
    referenceSurface->GetPoint(referencePointId, referencePoint);

    xForce[i] = referencePoint[0] - targetPoint[0];
    yForce[i] = referencePoint[1] - targetPoint[1];
    zForce[i] = referencePoint[2] - targetPoint[2];

    /*
    if( CheckGlobalDebugLevel( 6 ) )
    {
    args->normSum += vcl_sqrt(
        vnl_math_sqr( xForce[i] ) +
        vnl_math_sqr( yForce[i] ) +
        vnl_math_sqr( zForce[i] ) );
    }
    */
  }

  return IMI_THREAD_RETURN_VALUE;
}


///////////////////////////////////////////////////
//
//  ComputeRMSValue()
//
///////////////////////////////////////////////////
double imiSurfaceConstrainedDiffusionRegistrationFilter
::ComputeRMSValue()
{
  std::cout << " ComputeRMSValue()... " << std::endl;

  bool bFirstCheck = false;
  double meanDist = 999999;
  double sumDist = 0;

  if( m_LastXForce == NULL )
  {
    bFirstCheck = true;
    m_LastXForce = new double[m_NumberOfForceVectors];
  }
  if( m_LastYForce == NULL )
  {
    bFirstCheck = true;
    m_LastYForce = new double[m_NumberOfForceVectors];
  }
  if( m_LastZForce == NULL )
  {
    bFirstCheck = true;
    m_LastZForce = new double[m_NumberOfForceVectors];
  }

  if( !bFirstCheck && m_iCurrentIteration > 0 )
  {
    // compute difference to last update
    double dist2;
    for( int i = 0; i < m_NumberOfForceVectors; i++ )
    {
      dist2 = vnl_math_sqr( m_LastXForce[i] - m_XForceSmooth[i] ) +
              vnl_math_sqr( m_LastYForce[i] - m_YForceSmooth[i] ) +
              vnl_math_sqr( m_LastZForce[i] - m_ZForceSmooth[i] );
      sumDist += vcl_sqrt( dist2 );

      //copy forces
      m_LastXForce[i] = m_XForceSmooth[i];
      m_LastYForce[i] = m_YForceSmooth[i];
      m_LastZForce[i] = m_ZForceSmooth[i];

    }
    meanDist = sumDist / m_NumberOfForceVectors;

  }
  else
  {
    // for first time take length of force vectors
    double dist2;
    for (int k = 0; k < m_NumberOfForceVectors; k++)
    {
      dist2 = vnl_math_sqr( m_XForceSmooth[k] ) +
              vnl_math_sqr( m_YForceSmooth[k] ) +
              vnl_math_sqr( m_ZForceSmooth[k] );

      sumDist += vcl_sqrt( dist2 );

      // copy forces
      m_LastXForce[k] = m_XForceSmooth[k];
      m_LastYForce[k] = m_YForceSmooth[k];
      m_LastZForce[k] = m_ZForceSmooth[k];
    }
    meanDist = sumDist / m_NumberOfForceVectors;
  }

  std::cout << " Distance: " << meanDist << "(mean) " << sumDist << "(sum)" << std::endl;

  return meanDist;
}


///////////////////////////////////////////////////
//
//  InitializeDeformedTargetSurfaceByInitialTransform()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::InitializeDeformedTargetSurfaceByInitialTransform()
{
  std::cout << " InitializeDeformedTargetSurfaceByInitialTransform ... " << std::endl;

  if( m_spInitialTransform.GetPointer() == NULL )
  {
    std::cerr << "No initial transform set!" << std::endl;
    return false;
  }

  vtkSmartPointer<vtkPolyData> targetSurface = GetTargetSurface();

  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter;

  transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetInputData( targetSurface );
  transformFilter->SetTransform( m_spInitialTransform );
  transformFilter->Update();

  vtkSmartPointer<vtkPolyData> transformedTargetSurface =
      transformFilter->GetOutput();

  //small check
  if( transformedTargetSurface->GetNumberOfPoints() != m_NumberOfForceVectors
      || targetSurface->GetNumberOfPoints() != m_NumberOfForceVectors )
  {
    std::cerr <<  "Something wrong by initial transform!" << std::endl;
    return false;
  }

  // we initialize now the deformed target surface with this transformed surface
  m_spDeformedTargetSurface = vtkSmartPointer<vtkPolyData>::New();
  m_spDeformedTargetSurface->DeepCopy( transformedTargetSurface );

  return true;
}


///////////////////////////////////////////////////
//
//  SetInitialTransform()
//
///////////////////////////////////////////////////
void imiSurfaceConstrainedDiffusionRegistrationFilter
::SetInitialTransform(
    AffineTransformPointerType affineITKTransform)
{
  std::cout << "  SetInitialTransform from ITK affine transform ..." << std::endl;
  typedef AffineTransformType::ParametersType AffineTransformParametersType;

  AffineTransformParametersType parameters =
      affineITKTransform->GetParameters();

  vtkMatrix4x4* matrix = vtkMatrix4x4::New();
  matrix->Identity();

  int k = 0;
  for( int i = 0; i < 3; i++)
  {
    for( int j = 0; j < 3; j++)
    {
      matrix->SetElement( i, j, parameters[k] );
      k++;
    }
    matrix->SetElement( i, 3, parameters[9 + i] );
  }

  m_spInitialTransform = vtkSmartPointer<vtkTransform>::New();
  m_spInitialTransform->Identity();
  m_spInitialTransform->SetMatrix( matrix );
}


///////////////////////////////////////////////////
//
//  SetInitialTransform()
//
///////////////////////////////////////////////////
void imiSurfaceConstrainedDiffusionRegistrationFilter
::SetInitialTransform(
    vtkSmartPointer<vtkTransform> initialTransform)
{
  m_spInitialTransform = initialTransform;
}


///////////////////////////////////////////////////
//
//  GetInitialTransform()
//
///////////////////////////////////////////////////
vtkSmartPointer<vtkTransform> imiSurfaceConstrainedDiffusionRegistrationFilter
::GetInitialTransform()
{
  return m_spInitialTransform;
}


///////////////////////////////////////////////////
//
//  BuildNeighbourhoodIDList()
//
///////////////////////////////////////////////////
bool imiSurfaceConstrainedDiffusionRegistrationFilter
::BuildNeighbourhoodIDList(
    vtkSmartPointer<vtkPolyData> surface,
    double sigma,
    std::vector<NeighborhoodIDListType> &idListVector,
    std::vector<NeighborhoodWeightListType> &weightVector )
{

  //
  // define a limit for the size of the neighboorhood lists
  // to avoid swapping (set it to 1GB)
  //
#define NEIGHBORHOOD_LIST_MEMORY_LIMIT 1073741824

  if( surface.GetPointer() == NULL )
  {
    return false;
  }

  if( surface->GetNumberOfPoints() < 10 )
  {
    std::cout << "BuildNeighbourhoodIDList: Not enough surface points!" << std::endl;
    return false;
  }

  double sigma2 = vnl_math_sqr( sigma );
  double radius = 2 * sigma;

  vtkSmartPointer<vtkPointLocator> surfacePointLocator;

  if( surface == GetReferenceSurface() )
  {
    std::cout << " Build neighborhood ID list for reference surface... " << std::endl;
    surfacePointLocator = m_spReferencePointLocator;
  }
  else if( surface == GetTargetSurface() )
  {
    std::cout << " Build neighborhood ID list for target surface... " << std::endl;
    surfacePointLocator = m_spTargetPointLocator;
  }
  else
  {
    std::cerr <<  "Build neighborhood ID list only implemented for reference and target surface!" << std::endl;
    return false;
  }

  const unsigned int numberOfPoints = surface->GetNumberOfPoints();

  idListVector.clear();
  weightVector.clear();

  vtkIdList *neighbourIdsList = vtkIdList::New();

  double surfacePoint[3];
  unsigned long totalPointsCounter = 0;
  for( unsigned int i = 0; i < numberOfPoints; i++ )
  {
    // get the current surface point
    surface->GetPoint( i, surfacePoint );

    // get a list with Ids in the neighbourhood
    neighbourIdsList->Reset();
    surfacePointLocator->FindPointsWithinRadius( radius, surfacePoint, neighbourIdsList );
    int foundPoints = neighbourIdsList->GetNumberOfIds();
    if( foundPoints < 1 )
    {
      std::cerr <<  "BuildNeighbourhoodIDList(): No neighborhood points in radius " << radius << "found!" << std::endl;
      std::cerr <<  "Increase sigma!" << std::endl;
      return false;
    }

    // make a list with the gaussian weights
    NeighborhoodIDListType idList;
    idList.clear();
    NeighborhoodWeightListType weightList;
    weightList.clear();

    double dist2;
    double gaussFac;
    double gaussSum = 0;
    double neighbourhoodPoint[3];
    for( int j = 0; j < foundPoints; j++ )
    {
      int id = neighbourIdsList->GetId( j );
      if( id >= static_cast<int>(numberOfPoints) || id < 0 )
      {
        std::cerr <<  "Invalid Id: " << id << std::endl;
        continue;
      }

      surface->GetPoint( id, neighbourhoodPoint );

      dist2 = vtkMath::Distance2BetweenPoints( surfacePoint, neighbourhoodPoint );
      if( dist2 > vnl_math_sqr( radius ) )
      {
        std::cerr <<  "Distance > Radius: " << dist2 << std::endl;
        continue;
      }

      gaussFac = vcl_exp( -(dist2 / (2.0 * sigma2)) );
      gaussSum += gaussFac;

      idList.push_back( static_cast<NeighborhoodIDType>(id) );
      weightList.push_back( static_cast<NeighborhoodWeightType>(gaussFac) );
      totalPointsCounter++;
    }

    if( idList.size() != weightList.size() )
    {
      std::cerr << "BuildNeighbourhoodIDList(): Number of IDs and weights do not fit for point id "<<i<<"!" << std::endl;
      return false;
    }
    idListVector.push_back( idList );
    weightVector.push_back( weightList );

    double currentMemoryAllocated = totalPointsCounter * sizeof(float) + totalPointsCounter * sizeof(int);
    if(currentMemoryAllocated > NEIGHBORHOOD_LIST_MEMORY_LIMIT )
    {
      std::cerr << "BuildNeighbourhoodIDList(): the current memory allocated is larger than "<<NEIGHBORHOOD_LIST_MEMORY_LIMIT<<" !" << std::endl;
      std::cerr << "BuildNeighbourhoodIDList(): We can not pre-compute the NeighbourhoodIDList for large surfaces !" << std::endl;
      std::cerr << "BuildNeighbourhoodIDList(): set UsePrecomputedNeighborhoodLists to false!" << std::endl;

      idListVector.clear();
      weightVector.clear();

      return false;
    }

  }

  if( idListVector.size() != numberOfPoints
      || weightVector.size() != numberOfPoints )
  {
    std::cerr << "BuildNeighbourhoodIDList(): Number of ID-lists or weight-lists do not match number of points!" << std::endl;
    return false;
  }

  neighbourIdsList->Delete();

  // check allocated memory
  double totaleMemoryAllocated = totalPointsCounter * sizeof(float) + totalPointsCounter * sizeof(int);
  std::cout << "  Number of surface points: "<<numberOfPoints<<" ( "<<numberOfPoints/1000000.0<<" Mio)." << std::endl;
  std::cout << "  Total number of points in the neighborhood: "<<totalPointsCounter<<" ( "<<totalPointsCounter/1000000<<" Mio)." << std::endl;
  std::cout << "  Total memory allocated: "<<totaleMemoryAllocated<<" ( "<<totaleMemoryAllocated/(1024*1024)<<" MByte)." << std::endl;

  return true;
}
