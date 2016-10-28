/** \file imiSurfaceConstrainedDiffusionRegistrationFilter.h
 *
 *  \b Initial \b Author: Jan Ehrhardt \n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiSurfaceConstrainedDiffusionRegistrationFilter_h
#define __imiSurfaceConstrainedDiffusionRegistrationFilter_h

// System includes:
#include <string>
#include <fstream>
#include <vector>
#include <list>

// ITK includes:
#include "itkAffineTransform.h"

// VTK includes:
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkPointLocator.h"
#include "vtkTransform.h"

// Project includes:
#include "icnsSurfaceRegistrationFilter.h"
#include "icnsImageThreads.h"

namespace imi
{
/** \class imi::imiSurfaceConstrainedDiffusionRegistrationFilter
 *  \brief A base class for a affine ICP registration.
 *
 *  <small> <!--Copyright Information: -->
 *  \b Initial \b Author: Alexander Schmidt-Richberg, Jan Ehrhardt \n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck\n
 *  </small>
 ***************************************************************************/
class imiSurfaceConstrainedDiffusionRegistrationFilter : public imiSurfaceRegistrationFilter
{
public:
  /** type definitions. */
  typedef itk::AffineTransform<double, 3>  AffineTransformType;
  typedef AffineTransformType::Pointer             AffineTransformPointerType;

  /** \brief Create an instance of imiRegistrationFilter,
   *  use Delete() to destroy. */
  static imiSurfaceConstrainedDiffusionRegistrationFilter* New();

  void SetInitialTransform( AffineTransformPointerType initialITKTransform );
  void SetInitialTransform( vtkSmartPointer<vtkTransform> initialTransform );
  vtkSmartPointer<vtkTransform> GetInitialTransform();

  void SetSmoothingSigma( const double s ) { m_SmoothingSigma = s; }
  double GetSmoothingSigma() const { return m_SmoothingSigma; }

  /** \brief Select if pre-computed neighborhoods are used.
   *
   *  Instead of finding the relevant neighbors of each surface point in each
   *  smoothing step we pre-compute all neighbors and save the neighbor-IDs for
   *  each point in lists. This speeds up the smoothing step, however a huge
   *  amount of memory is needed and this will fail for large surfaces
   *  and/or large sigmas. */
  void SetUsePreComputedNeighborhoodLists( const bool flag ) { m_bUsePrecomputedNeighborhoodLists = flag; }
  bool GetUsePreComputedNeighborhoodLists() const { return m_bUsePrecomputedNeighborhoodLists; }

  virtual vtkSmartPointer<vtkPolyData> GetTransformedTargetSurface();

protected:
  typedef unsigned int                    NeighborhoodIDType;
  typedef std::vector<NeighborhoodIDType> NeighborhoodIDListType;

  typedef float                               NeighborhoodWeightType;
  typedef std::vector<NeighborhoodWeightType> NeighborhoodWeightListType;

  /** Constructor. */
  imiSurfaceConstrainedDiffusionRegistrationFilter();

  /** Destructor. */
  virtual ~imiSurfaceConstrainedDiffusionRegistrationFilter();

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

  bool ApplyForceFieldSurface(
      vtkSmartPointer<vtkPolyData> inputSurface,
      vtkSmartPointer<vtkPolyData> deformedInputSurface );

  virtual bool SmoothForceField();

  bool SmoothForceFieldSurface();

  virtual bool ComputeForceField();

  double ComputeRMSValue();

  bool InitializeDeformedTargetSurfaceByInitialTransform();

  bool BuildNeighbourhoodIDList(
      vtkSmartPointer<vtkPolyData> surface,
      double sigma,
      std::vector<NeighborhoodIDListType> &idListVector,
      std::vector<NeighborhoodWeightListType> &weightVector );


  /** An inner class to store parameters for multithreaded function call. */
  class imiSmoothPointsParameters: public icnsImageThreadParameters
  {
  public:
    // Lists only used for pre-computed neighborhood lists.
    std::vector<NeighborhoodIDListType>     idListVector;
    std::vector<NeighborhoodWeightListType> weightListVector;

    // Required if neighborhood lists are not pre-computed.
    vtkSmartPointer<vtkPolyData>     surface;
    vtkSmartPointer<vtkPointLocator> pointLocator;
    double sigma;

    // Smoothing forces.
    double *xForce;
    double *yForce;
    double *zForce;
    double *xForceSmooth;
    double *yForceSmooth;
    double *zForceSmooth;

    int numPoints;

    // Parameters for debug output.
    double normSum;
    double meanGaussSum;
    double meanFoundPts;
  };

  /** An inner class to store parameters for multithreaded function call. */
  class imiComputeForceFieldParameters: public icnsImageThreadParameters
  {
  public:
    // Surfaces.
    vtkSmartPointer<vtkPolyData>     referenceSurface;
    vtkSmartPointer<vtkPolyData>     targetSurface;
    vtkSmartPointer<vtkPolyData>     deformedTargetSurface;

    // The reference point locator.
    vtkSmartPointer<vtkPointLocator> referencePointLocator;

    // Smoothing forces.
    double *xForce;
    double *yForce;
    double *zForce;

    // Parameters for debug output.
    double normSum;
  };

  /** Method for multithreaded calculation of the final field. */
  static IMI_THREAD_RETURN_TYPE T_SmoothPointsPrecomputed( void *args );

  /** Method for multithreaded calculation of the final field. */
  static IMI_THREAD_RETURN_TYPE T_SmoothPoints( void *args );

  /** Method for multithreaded calculation of the final field. */
  static IMI_THREAD_RETURN_TYPE T_ComputeForceField( void *args );

  double m_SmoothingSigma;
  bool m_bPreComputeForceField;

  //
  // Surface data
  //
  vtkSmartPointer<vtkPolyData> m_spDeformedTargetSurface;

  // the initial transformation
  vtkSmartPointer<vtkTransform> m_spInitialTransform;

  //
  // Locators
  //
  vtkSmartPointer<vtkPointLocator> m_spTargetPointLocator;
  vtkSmartPointer<vtkPointLocator> m_spReferencePointLocator;

  int m_DeformFieldRepresent;
  int m_ProjectionMethod;

  double *m_XForce;
  double *m_YForce;
  double *m_ZForce;

  double *m_XForceSmooth;
  double *m_YForceSmooth;
  double *m_ZForceSmooth;

  double *m_LastXForce;
  double *m_LastYForce;
  double *m_LastZForce;
  int m_NumberOfForceVectors;

  //
  // to be fast for smoothing, pre-compute neighborhoods and weights
  //
  bool m_bUsePrecomputedNeighborhoodLists;
  std::vector<NeighborhoodIDListType>     m_TargetNeighborIdListVector;
  std::vector<NeighborhoodWeightListType> m_TargetNeighborWeightListVector;
};

}

#endif
