/** \file imiSurfaceRegistrationFilter.cxx
 *
 *  \b Initial \b Author: Jan Ehrhardt\n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// System includes:
#include <fstream>

// VTK includes:
#include "vtkMarchingCubes.h"
#include "vtkImageThreshold.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkProperty.h"
#include "vtkPolyDataWriter.h"
#include "vtkDecimatePro.h"
#include "vtkLandmarkTransform.h"
#include "vtkMatrix4x4.h"

// ITK includes:
//#include "itkRoundedLinearInterpolateImageFunction.h"

// Project includes:
#include "icnsSurfaceRegistrationFilter.h"



using namespace imi;

///////////////////////////////////////////////////////////////////////////////////////////
//
// Constructor
//
///////////////////////////////////////////////////////////////////////////////////////////
imiSurfaceRegistrationFilter::imiSurfaceRegistrationFilter()
{

  m_iNumberOfIterations = 50;
  m_iCurrentIteration = 0;
  m_dLastRMSValue = 0;
  m_bAbortRegistration = false;
  m_bIsInitialized = false;

  
  m_EpsilonDist = 0.0001;
  //m_EpsilonDist = 0.00001;
  //m_EpsilonDist = 0.00000001;

  m_sDebugPathName = "";
  m_sLogFileName = "";
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// Destructor

//
///////////////////////////////////////////////////////////////////////////////////////////
imiSurfaceRegistrationFilter::~imiSurfaceRegistrationFilter()
{}


///////////////////////////////////////////////////////////////////////////////////////////
//
// SetTargetSurface
//
///////////////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::SetTargetSurface( vtkSmartPointer<vtkPolyData> targetSurface )
{
  m_spTargetSurface = targetSurface;
  m_bIsInitialized = false;
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// SetReferenceImage
//
///////////////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::SetReferenceSurface( vtkSmartPointer<vtkPolyData> referenceSurface )
{
  m_spReferenceSurface = referenceSurface;
  m_bIsInitialized = false;
}


////////////////////////////////////////////////////////////////////////////////
//
// SetNumberOfIterations
//
////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::SetNumberOfIterations( const int iterations )
{
  if( iterations > 0 )
  {
    m_iNumberOfIterations = iterations;
  }
  else
  {
    m_iNumberOfIterations = 0;
  }
}


///////////////////////////////////////////////////
//
//  InitializeRegistration()
//
///////////////////////////////////////////////////
bool imiSurfaceRegistrationFilter::InitializeRegistration()
{
  if( m_bIsInitialized )
  {
    return true;
  }

  std::cout << " Initialize Surface Registration ... " << std::endl;

  if ( m_spReferenceSurface == NULL )
  {
    std::cerr << "No reference surface set!" << std::endl;
    return false;
  }
  if ( m_spTargetSurface == NULL )
  {
    std::cerr << "No target surface set!" << std::endl;
    return false;
  }

  // reset values
  m_iCurrentIteration = 0;
  m_bAbortRegistration = false;
  m_dLastRMSValue = 0;

  // Open log file if necessary.
  if( !m_sLogFileName.empty() )
  {
    std::cout << "  Open log file." << std::endl;
    m_LogFile.open( m_sLogFileName.c_str(), std::ios::app );
    if( !m_LogFile.is_open() )
    {
      std::cerr << "Can't open log file!" << std::endl;
      m_LogFile.close();
    }
  }

  // Write the header of the log file.
  WriteLogHeader();

  // Set initialized.
  m_bIsInitialized = true;

  // Return success of initialization.
  return true;
}

///////////////////////////////////////////////////
//
//  ExecuteRegistration()
//
///////////////////////////////////////////////////
bool imiSurfaceRegistrationFilter::ExecuteRegistration()
{
  std::cout << "Starting registration." << std::endl;

  bool bAbortRegistration = false;

  //////////////////////////////////////////////
  // Call initialization.
  if( !InitializeRegistration() )
  {
    std::cerr << "Initialization failed! Aborting registration." << std::endl;
    return false;
  }

  //////////////////////////////////////////////
  // Start iteration.
  for( m_iCurrentIteration = 0; m_iCurrentIteration < m_iNumberOfIterations; m_iCurrentIteration++ )
  {
    std::cout << "  Iteration " << m_iCurrentIteration + 1 << ":" << std::endl;

    if( !UpdateRegistration() )
    {
      std::cerr << "UpdateRegistration() failed!" << std::endl;
      return false;
    }

    std::cout << " SSD-Value:" << m_dLastRMSValue << std::endl;

    // Log iteration validation.
    WriteLogInformation( m_iCurrentIteration );

    // Check if registration aborted.
    bAbortRegistration = GetAbortRegistration();
    if( bAbortRegistration )
    {
      std::cerr << "Registration execution aborted! " << std::endl;
      break;
    }

    // Check if stop criterion fulfilled.
    if( CheckStopCriterion( m_iCurrentIteration, m_dLastRMSValue ) )
    {
      std::cout << "Stopped Registration after iteration " << m_iCurrentIteration << ".\n" << std::endl;
      break;
    }
    std::cout << "Iteration " << m_iCurrentIteration + 1 << " finished.\n" << std::endl;
  }

  // Finish registration.
  if( !FinishRegistration() )
  {
    std::cerr << "Finishing of registration failed!" << std::endl;
    return false;
  }

  return !bAbortRegistration;
}

////////////////////////////////////////////////////////////////////////////////
//
// FinishRegistration
//
////////////////////////////////////////////////////////////////////////////////
bool imiSurfaceRegistrationFilter::FinishRegistration()
{
  std::cout << "  Finish Registration ..." << std::endl;

  // Close log-file.
  if( m_LogFile.is_open() )
  {
    std::cout << "   Close Log-File." << std::endl;
    m_LogFile.close();
  }

  return true;
}

////////////////////////////////////////////////////////////////////////////////
//
// AbortRegistration()
//
////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::AbortRegistration()
{
  m_bAbortRegistration = true;
}

////////////////////////////////////////////////////////////////////////////////
//
// WriteSurface()
//
////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::WriteSurface(vtkSmartPointer<vtkPolyData> polys, const char* filename)

{
  std::cout << "  Write surface to " << filename << " ... " << std::endl;

  vtkPolyDataWriter* writer = vtkPolyDataWriter::New();
  writer->SetInputData( polys );
  writer->SetFileName( filename );
  writer->Write();
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// CheckStopCriterion
//
///////////////////////////////////////////////////////////////////////////////////////////
bool imiSurfaceRegistrationFilter::CheckStopCriterion( int it, double metricValue )
{
  if( metricValue <  m_EpsilonDist && it > 1 )
  {
    return true;
  }

  return false;
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// SetLogFileName
//
///////////////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::SetLogFileName( const char* logFileName )
{
  if( logFileName != NULL )
    m_sLogFileName = logFileName;
  else
    m_sLogFileName = "";
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// WriteLogHeader
//
///////////////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::WriteLogHeader()
{
  if( m_LogFile.is_open() )
  {
    m_LogFile << "\n# Start Surface Registration" << std::endl;
    m_LogFile << "#    reference surface points: " << GetReferenceSurface()->GetNumberOfPoints() << std::endl;
    m_LogFile << "#       target surface points: " << GetTargetSurface()->GetNumberOfPoints() << std::endl;
    m_LogFile << "#                  iterations: " << GetNumberOfIterations() << std::endl;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// WriteLogInformation
//
///////////////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::WriteLogInformation( int iteration )
{
  if( m_LogFile.is_open() )
  {
    m_LogFile << "\n" << iteration << "\t" << GetLastRMSValue() << std::flush;
  }
}

///////////////////////////////////////////////////////////////////////////////////////////
//
// WriteToLogFile
//
///////////////////////////////////////////////////////////////////////////////////////////
void imiSurfaceRegistrationFilter::WriteToLogFile( std::string s )
{
  if( m_LogFile.is_open() )
  {
    m_LogFile << s << std::flush;
  }
}
