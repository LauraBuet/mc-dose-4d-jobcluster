/** \file imiSurfaceRegistrationFilter.h
 *
 *  \b Initial \b Author: Jan Ehrhardt \n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiSurfaceRegistrationFilter_h
#define __imiSurfaceRegistrationFilter_h

// System includes:
#include <string>
#include <fstream>

// VTK includes:
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>

// Project includes:

namespace imi
{
/** \class imi::imiSurfaceRegistrationFilter
 *  \brief A base class for surface registration.
 *
 *  <small> <!--Copyright Information: -->
 *  \b Initial \b Author: Jan Ehrhardt \n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck\n
 *  </small>
 ***************************************************************************/
class imiSurfaceRegistrationFilter
{
public:
  //
  // Methods for set/get target reference etc
  //
  /** \brief Set the target surface.
   *  \param  targetSurface A pointer to the target surface. */
  void SetTargetSurface( vtkSmartPointer<vtkPolyData> targetSurface );
  vtkSmartPointer<vtkPolyData> GetTargetSurface() const {return m_spTargetSurface;}

  /** \brief Set the reference surface.
   *  \param  referenceSurface A pointer to the reference surface. */
  void SetReferenceSurface( vtkSmartPointer<vtkPolyData> referenceSurface );
  vtkSmartPointer<vtkPolyData> GetReferenceSurface() const {return m_spReferenceSurface;}

  /** \brief Set number of iterations.
   *  \param  iterations The number of iterations. */
  void SetNumberOfIterations( const int iterations );
  int GetNumberOfIterations() const {return m_iNumberOfIterations;}

  //
  // Methods for execution of registration
  /** \brief Execute registration
   *
   *  Executes registration. First calls initialize and then UpdateRegistration()
   *  for  NumberOfIterations iterations.
   *  \return true or false. */
  virtual bool ExecuteRegistration();

  /** \brief Get the number of the current iteration. */
  int GetCurrentIteration() const { return m_iCurrentIteration; }

  /** \brief Get the RMS value of the last iteration. */
  double GetLastRMSValue() const { return m_dLastRMSValue; }

  virtual vtkSmartPointer<vtkPolyData> GetTransformedTargetSurface() = 0;

  /** \brief Abort registration at end of current iteration. */
  void AbortRegistration();

  /** \brief Returns true, if registration was aborted. */
  bool GetAbortRegistration() const { return m_bAbortRegistration; }

  //
  // Log methods.
  /** Set a filename for the log output. If no filename is set, no log-file is written. */
  void SetLogFileName(const char* logFileName);
  /** Write iteration summary to the log-file (if open). */
  virtual void WriteLogInformation(int iteration);
  /** Write this string to the log-file (if open). */
  virtual void WriteToLogFile(std::string s);
  /** Write header information to the log-file (if open). */
  virtual void WriteLogHeader();

  void SetDebugPathName( std::string pathName ) { m_sDebugPathName = pathName; }

  void WriteSurface(vtkSmartPointer<vtkPolyData> polys, const char* filename);

protected:
  /** Constructor. */
  imiSurfaceRegistrationFilter();

  /** Destructor. */
  virtual ~imiSurfaceRegistrationFilter();

  /** \brief Initialize the registration handler with the image data.
   *  \return Returns false if initialization fails (missing data). */
  virtual bool InitializeRegistration();

  /** \brief Execute one registration iteration.
   *
   *  Executes one iteration and updates the values of the registration and
   *  returns the RMS change.
   *  \return Success of the update. */
  virtual bool UpdateRegistration() = 0;

  /** \brief Check if the stop-criterion is fulfilled.
   *  Uses m_pStopCriterion to test weather the stop criterion is fulfilled. */
  virtual bool CheckStopCriterion(int currentIteration, double metricValue);

  /** \brief Finish the execution of the registration.
   *  \return Success of finishing procedure. */
  virtual bool FinishRegistration();

  //
  // Surface data
  //
  vtkSmartPointer<vtkPolyData> m_spReferenceSurface;
  vtkSmartPointer<vtkPolyData> m_spTargetSurface;
  //
  // Other member variables
  bool m_bIsInitialized;
  double m_dLastRMSValue;
  int m_iNumberOfIterations;
  int m_iCurrentIteration;
  bool m_bAbortRegistration;
  double m_EpsilonDist;

  std::string m_sDebugPathName;
  std::string m_sLogFileName;
  std::ofstream m_LogFile;
};

}

#endif
