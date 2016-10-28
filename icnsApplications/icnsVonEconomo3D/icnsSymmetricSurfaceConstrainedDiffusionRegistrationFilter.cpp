/** \file imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter.cxx
 *
 *  \b Initial \b Author: Jan Ehrhardt\n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#include "vtkPointLocator.h"
#include "vtkCellLocator.h"
#include "vtkMath.h"
#include "vtkTriangle.h"
#include "vtkTriangleFilter.h"
#include "vtkCellArray.h"
#include "vtkTransformPolyDataFilter.h"

#include "icnsSymmetricSurfaceConstrainedDiffusionRegistrationFilter.h"

using namespace imi;

///////////////////////////////////////////////////////////////////////////////////////////
//
// imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::New()
//
///////////////////////////////////////////////////////////////////////////////////////////
imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter* imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::New()
{
  return new imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter();
}

imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter()
{
  m_bCurrentDirection = true;

  m_spDeformedReferenceSurface = NULL;

  m_XReferenceForce = NULL;
  m_YReferenceForce = NULL;
  m_ZReferenceForce = NULL;
  m_XReferenceForceSmooth = NULL;
  m_YReferenceForceSmooth = NULL;
  m_ZReferenceForceSmooth = NULL;
  m_NumberOfReferenceForceVectors = 0;

  m_ReferenceForceMagnitude = NULL;
  m_TargetForceMagnitude = NULL;
}

imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::~imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter()
{
  if (!ClearForces())
  {
    std::cerr << "Can not free forces memory!" << std::endl;
  }
}

bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::ClearForces()
{
  if (!imiSurfaceConstrainedDiffusionRegistrationFilter::ClearForces())
  {
    std::cerr << "imiSurfaceConstrainedDiffusionRegistrationFilter::ClearForces() failed!" << std::endl;
    return false;
  }

  if (m_XReferenceForce != NULL)
    delete m_XReferenceForce;
  if (m_YReferenceForce != NULL)
    delete m_YReferenceForce;
  if (m_ZReferenceForce != NULL)
    delete m_ZReferenceForce;

  if (m_XReferenceForceSmooth != NULL)
    delete m_XReferenceForceSmooth;
  if (m_YReferenceForceSmooth != NULL)
    delete m_YReferenceForceSmooth;
  if (m_ZReferenceForceSmooth != NULL)
    delete m_ZReferenceForceSmooth;

  if (m_ReferenceForceMagnitude != NULL)
    delete m_ReferenceForceMagnitude;
  if (m_TargetForceMagnitude != NULL)
    delete m_TargetForceMagnitude;

  m_XReferenceForce = NULL;
  m_YReferenceForce = NULL;
  m_ZReferenceForce = NULL;
  m_XReferenceForceSmooth = NULL;
  m_YReferenceForceSmooth = NULL;
  m_ZReferenceForceSmooth = NULL;
  m_ReferenceForceMagnitude = NULL;
  m_TargetForceMagnitude = NULL;

  m_NumberOfReferenceForceVectors = 0;

  return true;
}

vtkSmartPointer<vtkPolyData> imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::GetTransformedReferenceSurface()
{
  if (GetReferenceSurface() == NULL)
    return NULL;

  if (m_spDeformedReferenceSurface == NULL)
  {
    m_spDeformedReferenceSurface = vtkSmartPointer<vtkPolyData>::New();
    m_spDeformedReferenceSurface->DeepCopy(GetReferenceSurface());
  }
  return m_spDeformedReferenceSurface;
}

///////////////////////////////////////////////////
//
//  InitializeRegistration()
//
///////////////////////////////////////////////////
bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::InitializeRegistration()
{
  if( m_bIsInitialized )
  {
    return true;
  }

  if( !imiSurfaceConstrainedDiffusionRegistrationFilter::InitializeRegistration() )
  {
    return false;
  }

  std::cout << " imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::InitializeRegistration() ..." << std::endl;

  if( m_spDeformedReferenceSurface == NULL )
  {
    vtkSmartPointer<vtkPolyData> dummy = GetTransformedReferenceSurface();
    if( dummy == NULL )
    {
      std::cerr << "Can't create deformed Reference !" << std::endl;
      return false;
    }
  }

  // memory for the forces will be deleted by call of the virtual function ClearForces from
  // imiSurfaceConstrainedDiffusionRegistrationFilter::InitializeRegistration()

  // delete neighbourhood-lists
  m_ReferenceNeighbourIdListVector.clear();
  m_ReferenceNeighbourWeightListVector.clear();

  std::cout << "  Allocate reference force array ..." << std::endl;
  // allocate space for reference forces
  m_NumberOfReferenceForceVectors = GetReferenceSurface()->GetNumberOfPoints();
  m_XReferenceForce = new double[m_NumberOfReferenceForceVectors];
  m_YReferenceForce = new double[m_NumberOfReferenceForceVectors];
  m_ZReferenceForce = new double[m_NumberOfReferenceForceVectors];
  m_XReferenceForceSmooth = new double[m_NumberOfReferenceForceVectors];
  m_YReferenceForceSmooth = new double[m_NumberOfReferenceForceVectors];
  m_ZReferenceForceSmooth = new double[m_NumberOfReferenceForceVectors];

  // allocate space for force magnitude arrays
  m_ReferenceForceMagnitude = new double[m_NumberOfReferenceForceVectors];
  m_TargetForceMagnitude = new double[m_NumberOfForceVectors];

  std::cout << "   Setting initial reference force field to zero ..." << std::endl;
  // initialize reference forces
  for( int i = 0; i < m_NumberOfReferenceForceVectors; i++ )
  {
    m_XReferenceForce[i] = 0;
    m_YReferenceForce[i] = 0;
    m_ZReferenceForce[i] = 0;
    m_XReferenceForceSmooth[i] = 0;
    m_YReferenceForceSmooth[i] = 0;
    m_ZReferenceForceSmooth[i] = 0;
  }

  if( m_spInitialTransform.GetPointer() != NULL )
  {
    if( !InitializeDeformedReferenceSurfaceByInitialTransform() )
    {
      std::cerr << "InitializeDeformedTargetSurfaceByInitialTransform() failed !" << std::endl;
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
bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::UpdateRegistration()
{
  std::cout << " Update SymmetricSurfaceConstrainedDiffusionRegistrationFilter ..." << std::endl;

  /*
  if( CheckGlobalDebugLevel( 9 ) )
  {
    char filename[200];
    sprintf( filename, "%sIteration_%.3d_DeformedTargetSurface.vtk",
        m_sDebugPathName.c_str(), m_iCurrentIteration + 1 );
    imiDEBUGINFO( 9, "    Writing deformed target surface to file " << filename );

    WriteSurface( GetTransformedTargetSurface(), filename );

    sprintf( filename, "%sIteration_%.3d_DeformedReferenceSurface.vtk",
        m_sDebugPathName.c_str(), m_iCurrentIteration + 1 );
    imiDEBUGINFO( 9, "    Writing deformed reference surface to file " << filename );

    WriteSurface( GetTransformedReferenceSurface(), filename );
  }
  */

  // Use initial force field in first iteration.
  if( !ComputeForceField() )
  {
    std::cerr << "Computing force field failed in iteration " << m_iCurrentIteration << "!" << std::endl;
    return false;
  }

  if( !SmoothForceField() )
  {
    std::cerr << "Smoothing force field failed in iteration " << m_iCurrentIteration << "!" << std::endl;
    return false;
  }

  if( !ApplyForceField() )
  {
    std::cerr << "Applying force field failed in iteration " << m_iCurrentIteration << "!" << std::endl;
    return false;
  }

  // Update only for iterations with the same direction.
  if( m_bCurrentDirection )
  {
    // Compute the mean change in the force vectors.
    m_dLastRMSValue = ComputeRMSValue();
  }

  return true;
}

///////////////////////////////////////////////////
//
//  ApplyForceField()
//
///////////////////////////////////////////////////
bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
::ApplyForceField()
{
  if( !ApplyForceFieldSurface( GetTargetSurface(), GetTransformedTargetSurface() ) )
  {
    std::cerr << "Applying target surface forces failed !" << std::endl;
    return false;
  }
  if( !ApplyReferenceForceFieldToSurface( GetReferenceSurface(), GetTransformedReferenceSurface() ) )
  {
    std::cerr << "Applying target surface forces failed !" << std::endl;
    return false;
  }
  return true;
}

///////////////////////////////////////////////////
//
//  ApplyForceFieldSurface()
//
///////////////////////////////////////////////////
bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
::ApplyReferenceForceFieldToSurface(
    vtkSmartPointer<vtkPolyData> inputSurface,
    vtkSmartPointer<vtkPolyData> deformedInputSurface)
{
  std::cout << " ApplyForceFieldSurface ...  " << std::endl;

  if( inputSurface->GetNumberOfPoints() != deformedInputSurface->GetNumberOfPoints() )
  {
    std::cerr << 
        "imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::ApplyReferenceForceFieldToSurface() "
        "Input and deformed have different number of points!" << std::endl;
    return false;
  }

  if( inputSurface->GetNumberOfPoints() != m_NumberOfReferenceForceVectors )
  {
    std::cerr << 
        "imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter::ApplyReferenceForceFieldToSurface() "
        "Input and force array have different number of points!" << std::endl;
    return false;
  }

  double inputPoint[3];

  int numPoints = inputSurface->GetNumberOfPoints();
  double normSum = 0.0;

  for (int i = 0; i < numPoints; i++)
  {
    inputSurface->GetPoint(i, inputPoint);

    inputPoint[0] += m_XReferenceForceSmooth[i];
    inputPoint[1] += m_YReferenceForceSmooth[i];
    inputPoint[2] += m_ZReferenceForceSmooth[i];

    normSum += vcl_sqrt(
        vnl_math_sqr( m_XReferenceForceSmooth[i] ) +
        vnl_math_sqr( m_YReferenceForceSmooth[i] ) +
        vnl_math_sqr( m_ZReferenceForceSmooth[i] ) );

    deformedInputSurface->GetPoints()->SetPoint(i, inputPoint);
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
bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
::SmoothForceField()
{
  if( !SmoothForceFieldSurface() )
  {
    std::cerr << "Smoothing target surface forces failed!" << std::endl;
    return false;
  }

  if( !SmoothReferenceForceFieldSurface() )
  {
    std::cerr << "Smoothing reference surface forces failed!" << std::endl;
    return false;
  }

  return true;
}

///////////////////////////////////////////////////
//
//  SmoothReferenceForceFieldSurface()
//
///////////////////////////////////////////////////
bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
::SmoothReferenceForceFieldSurface()
{
  std::cout << " SmoothReferenceForceFieldSurface (s=" << GetSmoothingSigma() << ") ... " << std::endl;

  double meanGaussSum = 0;
  double meanFoundPts = 0;
  double normSum = 0;

  vtkSmartPointer<vtkPolyData> referenceSurface = GetReferenceSurface();
  unsigned int numPoints = referenceSurface->GetNumberOfPoints();

  //
  // Pre-compute neighborhood lists if requested.
  //
  if( m_bUsePrecomputedNeighborhoodLists )
  {
    //
    // Check if BuildNeighbourhoodIDList() has to be called.
    //
    if( GetCurrentIteration() == 0 || m_ReferenceNeighbourIdListVector.size() != numPoints )
    {
      if( !BuildNeighbourhoodIDList(
               referenceSurface,
               GetSmoothingSigma(),
               m_ReferenceNeighbourIdListVector,
               m_ReferenceNeighbourWeightListVector ) )
      {
        std::cerr << "SmoothReferenceForceFieldSurface(): Can not build neighborhood ID Lists !" << std::endl;
        return false;
      }
      // make a small check
      if( m_ReferenceNeighbourIdListVector.size() != numPoints ||
          m_ReferenceNeighbourWeightListVector.size() != numPoints )
      {
        std::cerr << "SmoothReferenceForceFieldSurface(): Building neighborhood ID Lists failed !" << std::endl;
        return false;
      }
    }
  }

  // Define thread parameters.
  icnsImageThreads* imageThreads = new icnsImageThreads();
  unsigned int numberOfThreads = imageThreads->GetNumberOfThreads();

  // ==========================================
  // Calculate image to be smoothed from vector components.
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
    forceParameters->idListVector = m_ReferenceNeighbourIdListVector;
    forceParameters->weightListVector = m_ReferenceNeighbourWeightListVector;
    forceParameters->surface = referenceSurface;
    forceParameters->pointLocator = m_spReferencePointLocator;
    forceParameters->sigma = this->GetSmoothingSigma();
    forceParameters->xForce = m_XReferenceForce;
    forceParameters->yForce = m_YReferenceForce;
    forceParameters->zForce = m_ZReferenceForce;
    forceParameters->xForceSmooth = m_XReferenceForceSmooth;
    forceParameters->yForceSmooth = m_YReferenceForceSmooth;
    forceParameters->zForceSmooth = m_ZReferenceForceSmooth;
    forceParameters->numPoints = numPoints;

    forceParameters->meanGaussSum = 0.0;
    forceParameters->meanFoundPts = 0.0;
    forceParameters->normSum = 0.0;

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
    forceParameters = ((imiSmoothPointsParameters*)
        imageThreads->GetThreadParameters( threadNr ));

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


///////////////////////////////////////////////////
//
//  ComputeForceField()
//
///////////////////////////////////////////////////
bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
::ComputeForceField()
{
  std::cout << "  Computing symmetric force field..." << std::endl;

  // Initialize reference force magnitude. TODO right reference force
  int numReferencePoints = this->GetReferenceSurface()->GetNumberOfPoints();
  int numTargetPoints = this->GetTargetSurface()->GetNumberOfPoints();
  if( m_bCurrentDirection )
  {
    for( int i = 0; i < numReferencePoints; i++ )
    {
      m_ReferenceForceMagnitude[i] = -1.0;
    }
  }
  else
  {
    for( int i = 0; i < numTargetPoints; i++ )
    {
      m_TargetForceMagnitude[i] = -1.0;
    }
  }

  // Define thread parameters.
  icnsImageThreads* imageThreads = new icnsImageThreads();
  unsigned int numberOfThreads = imageThreads->GetNumberOfThreads();

  //
  // CALCULATE FORCES FOR POINTS ON TARGET SURFACE
  //
  int numPoints = (m_bCurrentDirection) ? numTargetPoints : numReferencePoints;
  unsigned int threadRange = numPoints / numberOfThreads;

  imiComputeSymmetricForceFieldParameters* parameters;
  for( unsigned int threadNr = 0; threadNr < numberOfThreads; ++threadNr )
  {
    // Initialize thread parameters.
    parameters = new imiComputeSymmetricForceFieldParameters();

    // Calculate thread range.
    parameters->from = threadNr * threadRange;
    parameters->to = (threadNr == numberOfThreads - 1) ?
        numPoints : (threadNr + 1) * threadRange;

    // Initializing thread parameters.
    if( m_bCurrentDirection )
    {
      parameters->referenceSurface = this->GetReferenceSurface();
      parameters->targetSurface = this->GetTargetSurface();
      parameters->deformedReferenceSurface = this->GetTransformedReferenceSurface();
      parameters->deformedTargetSurface = this->GetTransformedTargetSurface();
      parameters->referencePointLocator = m_spReferencePointLocator;
      parameters->targetPointLocator = m_spTargetPointLocator;
      parameters->xForce = m_XForce;
      parameters->yForce = m_YForce;
      parameters->zForce = m_ZForce;
      parameters->xReferenceForce = m_XReferenceForce;
      parameters->yReferenceForce = m_YReferenceForce;
      parameters->zReferenceForce = m_ZReferenceForce;

      parameters->referenceForceMagnitude = m_ReferenceForceMagnitude;
    }
    else
    {
      parameters->referenceSurface = this->GetTargetSurface();
      parameters->targetSurface = this->GetReferenceSurface();
      parameters->deformedReferenceSurface = this->GetTransformedTargetSurface();
      parameters->deformedTargetSurface = this->GetTransformedReferenceSurface();
      parameters->referencePointLocator = m_spTargetPointLocator;
      parameters->targetPointLocator = m_spReferencePointLocator;
      parameters->xForce = m_XReferenceForce;
      parameters->yForce = m_YReferenceForce;
      parameters->zForce = m_ZReferenceForce;
      parameters->xReferenceForce = m_XForce;
      parameters->yReferenceForce = m_YForce;
      parameters->zReferenceForce = m_ZForce;

      parameters->referenceForceMagnitude = m_TargetForceMagnitude;
    }

    parameters->normSum = 0.0;

    // Set thread method.
    imageThreads->SetThreadMethod( threadNr, &T_ComputeSymmetricForceField, parameters );
  }

  // Execute the threads.
  imageThreads->ExecuteThreads();

  // Clean up.
  double normSum = 0.0;
  for( unsigned int threadNr = 0; threadNr < numberOfThreads; ++threadNr )
  {
    // Get thread parameters.
    parameters = (imiComputeSymmetricForceFieldParameters*)
        imageThreads->GetThreadParameters( threadNr );

    // Get output of thread for debug information.
    normSum += parameters->normSum;
  }

  std::cout << "  Norm (Target)   : " << normSum / numPoints << std::endl;

  //
  // CALCULATE FORCES FOR REMAINING POINTS ON REFERENCE SURFACE
  //
  numPoints = (m_bCurrentDirection) ? numReferencePoints : numTargetPoints;
  threadRange = numPoints / numberOfThreads;

  for( unsigned int threadNr = 0; threadNr < numberOfThreads; ++threadNr )
  {
    // Initialize thread parameters.
    parameters = (imiComputeSymmetricForceFieldParameters*)
        imageThreads->GetThreadParameters( threadNr );

    // Calculate new thread range.
    parameters->from = threadNr * threadRange;
    parameters->to = (threadNr == numberOfThreads - 1) ?
        numPoints : (threadNr + 1) * threadRange;

    // Set norm sum back.
    parameters->normSum = 0.0;

    // Set thread method.
    imageThreads->SetThreadMethod( threadNr, &T_ComputeSymmetricForceFieldRemainingPoints, parameters );
  }

  // Execute the threads.
  imageThreads->ExecuteThreads();

  // Clean up.
  normSum = 0.0;
  for( unsigned int threadNr = 0; threadNr < numberOfThreads; ++threadNr )
  {
    // Get thread parameters.
    parameters = (imiComputeSymmetricForceFieldParameters*)
        imageThreads->GetThreadParameters( threadNr );

    // Get output of thread for debug information.
    normSum += parameters->normSum;

    // Delete thread parameters.
    delete imageThreads->GetThreadParameters( threadNr );
  }
  delete imageThreads;

  std::cout << "  Norm (Reference): " << normSum / numPoints << std::endl;

  // Swap direction.
  m_bCurrentDirection = !m_bCurrentDirection;

  return true;
}


///////////////////////////////////////////////////
//
//  T_ComputeSymmetricForceField()
//
///////////////////////////////////////////////////
IMI_THREAD_RETURN_TYPE imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
::T_ComputeSymmetricForceField( void *vargs )
{
  imiComputeSymmetricForceFieldParameters* args = (imiComputeSymmetricForceFieldParameters*) vargs;

  // Read parameters for thread. Target and reference might be swapped,
  // depending on current direction.
  vtkPolyData* targetSurface = args->targetSurface;
  vtkPolyData* referenceSurface = args->referenceSurface;
  vtkPolyData* deformedTargetSurface = args->deformedTargetSurface;
  //vtkPolyData* deformedReferenceSurface = args->deformedReferenceSurface;

  //vtkPointLocator* targetPointLocator = args->targetPointLocator;
  vtkPointLocator* referencePointLocator = args->referencePointLocator;

  double *xTargetForce = args->xForce;
  double *yTargetForce = args->yForce;
  double *zTargetForce = args->zForce;

  double *xReferenceForce = args->xReferenceForce;
  double *yReferenceForce = args->yReferenceForce;
  double *zReferenceForce = args->zReferenceForce;

  double *referenceForceMagnitude = args->referenceForceMagnitude;

  // Other variable declaration.
  int referencePointId;
  double referencePoint[3];
  double targetPoint[3];
  double deformedTargetPoint[3];

  // Go through all points of the target surface, compute nearest neighbor of
  // current deformed target. Compute deformation vector from original target
  // to nearest neighbor and set this vector as current force. Set the inverse
  // vector as force to the point on the reference surface.
  //std::cout << " From " << args->from  << " To " << args->to << " Norm " << args->normSum << std::endl;
  for( unsigned int i = args->from; i < args->to; ++i )
  {
    targetSurface->GetPoint( i, targetPoint );
    deformedTargetSurface->GetPoint( i, deformedTargetPoint );

    // Find closest reference point to actual deformed target point.
    referencePointId = referencePointLocator->FindClosestPoint( deformedTargetPoint );

    if( referencePointId >= 0 )
    {
      referenceSurface->GetPoint( referencePointId, referencePoint );

      // Compute deformation vector to move target point to reference point.
      xTargetForce[i] = referencePoint[0] - targetPoint[0];
      yTargetForce[i] = referencePoint[1] - targetPoint[1];
      zTargetForce[i] = referencePoint[2] - targetPoint[2];

      // Compute magnitude of deformation vector and update sum.
      double forceMagnitude = vcl_sqrt(
          vnl_math_sqr( xTargetForce[i] ) +
          vnl_math_sqr( yTargetForce[i] ) +
          vnl_math_sqr( zTargetForce[i] ) );

      // Add the largest force to the corresponding reference point.
      if( forceMagnitude > referenceForceMagnitude[referencePointId] )
      {
        // This force is inverse to the target force.
        xReferenceForce[referencePointId] = targetPoint[0] - referencePoint[0];
        yReferenceForce[referencePointId] = targetPoint[1] - referencePoint[1];
        zReferenceForce[referencePointId] = targetPoint[2] - referencePoint[2];

        referenceForceMagnitude[referencePointId] = forceMagnitude;
      }

      // Return sum of magnitudes.
      args->normSum += forceMagnitude;
    }
    else
    {
      std::cerr << "No nearest neighbor found for target point ["
          << deformedTargetPoint[0] << ","
          << deformedTargetPoint[1] << ","
          << deformedTargetPoint[2] << "]!" << std::endl;

      // set forces to 0 for this point
      xTargetForce[i] = 0;
      yTargetForce[i] = 0;
      zTargetForce[i] = 0;
    }
  }

  return IMI_THREAD_RETURN_VALUE;
}


///////////////////////////////////////////////////
//
//  T_ComputeForceFieldNearestNeighbor()
//
///////////////////////////////////////////////////
IMI_THREAD_RETURN_TYPE imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
::T_ComputeSymmetricForceFieldRemainingPoints( void *vargs )
{
  imiComputeSymmetricForceFieldParameters* args = (imiComputeSymmetricForceFieldParameters*) vargs;

  // Read parameters for thread. Target and reference might be swapped,
  // depending on current direction.
  vtkPolyData* targetSurface = args->targetSurface;
  vtkPolyData* referenceSurface = args->referenceSurface;
  //vtkPolyData* deformedTargetSurface = args->deformedTargetSurface;
  vtkPolyData* deformedReferenceSurface = args->deformedReferenceSurface;

  vtkPointLocator* targetPointLocator = args->targetPointLocator;
  //vtkPointLocator* referencePointLocator = args->referencePointLocator;

  double *xReferenceForce = args->xReferenceForce;
  double *yReferenceForce = args->yReferenceForce;
  double *zReferenceForce = args->zReferenceForce;

  double *referenceForceMagnitude = args->referenceForceMagnitude;

  // Other variable declaration.
  int targetPointId;
  double referencePoint[3];
  double targetPoint[3];
  double deformedReferencePoint[3];

  // Go through all points of the reference surface, for all points where no
  // force was computed in the last loop. Compute nearest neighbor of current
  // deformed reference to target surface. Compute deformation vector from
  // original reference point to nearest neighbor and set this vector as current force.
  //std::cout << " From " << args->from  << " To " << args->to << " Norm " << args->normSum << std::endl;
  for( unsigned int i = args->from; i < args->to; ++i )
  {
    // Check if already a force is set for this point.
    if( referenceForceMagnitude[i] < 0 )
    {
      referenceSurface->GetPoint( i, referencePoint );
      deformedReferenceSurface->GetPoint( i, deformedReferencePoint );

      // Find closest reference point to actual deformed target point.
      targetPointId = targetPointLocator->FindClosestPoint( deformedReferencePoint );
      if( targetPointId >= 0 )
      {
        targetSurface->GetPoint( targetPointId, targetPoint );

        // Compute force for reference surface
        xReferenceForce[i] = targetPoint[0] - referencePoint[0];
        yReferenceForce[i] = targetPoint[1] - referencePoint[1];
        zReferenceForce[i] = targetPoint[2] - referencePoint[2];

        // Compute magnitude of deformation vector and update sum.
        double forceMagnitude = vcl_sqrt(
            vnl_math_sqr( xReferenceForce[i] ) +
            vnl_math_sqr( yReferenceForce[i] ) +
            vnl_math_sqr( zReferenceForce[i] ) );

        referenceForceMagnitude[i] = forceMagnitude;

        // Return sum of magnitudes.
        args->normSum += forceMagnitude;
      }
      else
      {
        std::cerr << "No nearest neighbor found for deformed reference point ["
            << deformedReferencePoint[0] << ","
            << deformedReferencePoint[1] << ","
            << deformedReferencePoint[2] << "]!" << std::endl;

        // set forces to 0 for this point
        xReferenceForce[i] = 0;
        yReferenceForce[i] = 0;
        zReferenceForce[i] = 0;

        referenceForceMagnitude[i] = 0;
      }
    }
  }

  return IMI_THREAD_RETURN_VALUE;
}

///////////////////////////////////////////////////
//
//  InitializeDeformedReferenceSurfaceByInitialTransform()
//
///////////////////////////////////////////////////
bool imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
::InitializeDeformedReferenceSurfaceByInitialTransform()
{
  std::cout << " InitializeDeformedTargetSurfaceByInitialTransform ... " << std::endl;

  if (m_spInitialTransform.GetPointer() == NULL)
  {
    std::cerr << "No initial transform set!" << std::endl;
    return false;
  }

  vtkSmartPointer<vtkPolyData> referenceSurface = GetReferenceSurface();

  vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter;

  transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformFilter->SetInputData(referenceSurface);
  transformFilter->SetTransform(m_spInitialTransform->GetInverse());
  transformFilter->Update();

  vtkSmartPointer<vtkPolyData> transformedReferenceSurface =
      transformFilter->GetOutput();

  //small check
  if( transformedReferenceSurface->GetNumberOfPoints() != m_NumberOfReferenceForceVectors ||
      referenceSurface->GetNumberOfPoints() != m_NumberOfReferenceForceVectors )
  {
    std::cerr << "Something wrong by initial transform!" << std::endl;
    return false;
  }

  // we initialize now the deformed target surface with this transformed surface
  m_spDeformedReferenceSurface = vtkSmartPointer<vtkPolyData>::New();
  m_spDeformedReferenceSurface->DeepCopy(transformedReferenceSurface);

  return true;
}
