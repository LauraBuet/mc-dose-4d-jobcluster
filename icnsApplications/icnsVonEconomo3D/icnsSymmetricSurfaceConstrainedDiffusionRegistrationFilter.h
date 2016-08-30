/** \file imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter.h
 *
 *  \b Initial \b Author: Jan Ehrhardt \n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter_h
#define __imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter_h

// System includes:
#include <string>
#include <fstream>

// ITK includes:

// VTK includes:
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPointLocator.h"

// Project includes:
#include "icnsSurfaceConstrainedDiffusionRegistrationFilter.h"

namespace imi
{
/** \class imi::imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter
 *  \brief A base class for a affine ICP registration.
 *
 *  <small> <!--Copyright Information: -->
 *  \b Initial \b Author: Alexander Schmidt-Richberg, Jan Ehrhardt \n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck\n
 *  </small>
 ***************************************************************************/
class imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter :
  public imiSurfaceConstrainedDiffusionRegistrationFilter
{
public:
  
  /** \brief Create an instance of imiRegistrationFilter,
   *  use Delete() to destroy. */
  static imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter* New();

  /** TODO
   *
   * \return */
  virtual vtkSmartPointer<vtkPolyData> GetTransformedReferenceSurface();

protected:
  /** Constructor. */
  imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter();

  /** Destructor. */
  virtual ~imiSymmetricSurfaceConstrainedDiffusionRegistrationFilter();

  /** free allocated arrays for force computation */
  virtual bool ClearForces();

  /** \brief Initialize the registration handler with the image data.
   *  \return Returns false if initialization fails (missing data). */
  virtual bool InitializeRegistration();

  /** \brief Execute one registration iteration.
   *
   *  Executes one iteration and updates the values of the registration and
   *  returns the RMS change.
   *  \return Success of the update.*/
  virtual bool UpdateRegistration();

  virtual bool ApplyForceField();

  bool ApplyReferenceForceFieldToSurface(
      vtkSmartPointer<vtkPolyData> inputSurface,
      vtkSmartPointer<vtkPolyData> deformedInputSurface );

  virtual bool SmoothForceField();

  bool SmoothReferenceForceFieldSurface();

  virtual bool ComputeForceField();

  bool ComputeSymmetricForceFieldByNearestNeighbour( bool bDirectionFlag );

  bool InitializeDeformedReferenceSurfaceByInitialTransform();

protected:
  /** An inner class to store parameters for multithreaded function call. */
  class imiComputeSymmetricForceFieldParameters: public imiComputeForceFieldParameters
  {
  public:
    // Surfaces.
    vtkSmartPointer<vtkPolyData>     deformedReferenceSurface;

    // The reference point locator.
    vtkSmartPointer<vtkPointLocator> targetPointLocator;

    // Reference smoothing forces.
    double *xReferenceForce;
    double *yReferenceForce;
    double *zReferenceForce;

    double *referenceForceMagnitude;

  };
  /** Method for multithreaded calculation of the final field. */
  static IMI_THREAD_RETURN_TYPE T_ComputeSymmetricForceField( void *args );

  /** Method for multithreaded calculation of the final field. */
  static IMI_THREAD_RETURN_TYPE T_ComputeSymmetricForceFieldRemainingPoints( void *args );

  bool m_bCurrentDirection;
  //
  // Surface data
  //
  vtkSmartPointer<vtkPolyData> m_spDeformedReferenceSurface;

  double *m_XReferenceForce;
  double *m_YReferenceForce;
  double *m_ZReferenceForce;

  double *m_XReferenceForceSmooth;
  double *m_YReferenceForceSmooth;
  double *m_ZReferenceForceSmooth;

  double *m_ReferenceForceMagnitude;
  double *m_TargetForceMagnitude;
  int m_NumberOfReferenceForceVectors;

  //
  // to be fast for smoothing, pre-compute neighbourhoods and weights
  //
  std::vector<NeighborhoodIDListType> m_ReferenceNeighbourIdListVector;
  std::vector<NeighborhoodWeightListType> m_ReferenceNeighbourWeightListVector;

};

}

#endif
