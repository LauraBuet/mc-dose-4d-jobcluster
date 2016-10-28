/** \file imiImageThreads.cxx
 *
 *  \b Initial \b Author: Alexander Schmidt-Richberg, Jan Ehrhardt\n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#include "icnsImageThreads.h"

unsigned int icnsImageThreads::m_iGlobalMaxThreads = IMI_DEFAULT_THREADNUM;

// Constructor, initializations
icnsImageThreads::icnsImageThreads()
{
  // Calculating number of processors depending on system.
  int num = 1;
#ifdef _WIN32
  {
    SYSTEM_INFO sysInfo;
    GetSystemInfo(&sysInfo);
    num = sysInfo.dwNumberOfProcessors;
  }
  //std::cout << "Number of processors: " << num << " (Win)" << std::endl;
#elif defined(__APPLE__)
  int mib[2] = {CTL_HW, HW_NCPU};
  size_t dataLen = sizeof(int); // 'num' is an 'int'
  int result = sysctl(mib, 2, &num, &dataLen, NULL, 0);
  if (result == -1)
  {
    num = 1;
  }
  //std::cout << "Number of processors: " << num << " (Apple)" << std::endl;
#elif defined(_SC_NPROCESSORS_ONLN)
  num = sysconf( _SC_NPROCESSORS_ONLN );
  //std::cout << "Number of processors: " << num << " (SC1)" << std::endl;
#elif defined(_SC_NPROC_ONLN)
  num = sysconf( _SC_NPROC_ONLN );
  //std::cout << "Number of processors: " << num << " (SC2)" << std::endl;
#endif
  m_iNumberOfThreads = num;

  if(m_iNumberOfThreads > m_iGlobalMaxThreads) m_iNumberOfThreads = m_iGlobalMaxThreads;

  for(int i=0; i<IMI_MAX_NUMBER_OF_THREADS; i++)
  {
    m_ThreadParameter[i] = NULL;
    m_ThreadMethods[i] = NULL;
  }
}

// Destructor, nothing to be done
icnsImageThreads::~icnsImageThreads(void)
{
//  for(int i=0; i<IMI_MAX_NUMBER_OF_THREADS; i++)
//  {
//    if( m_ThreadParameter[i] != NULL )
//    {
//      delete m_ThreadParameter[i];
//    }
//  }
}

void icnsImageThreads::SetGlobalMaxThreads(unsigned int num)
{
  if( num > IMI_MAX_NUMBER_OF_THREADS )
  {
    std::cout<<"\nWarning: Number of threads too big, using maximum of "
      << IMI_MAX_NUMBER_OF_THREADS << " instead." << std::endl;
    m_iGlobalMaxThreads = IMI_MAX_NUMBER_OF_THREADS;
  }
  else if(num <= 0)
  {
      std::cout<<"\nWarning: Number of threads must be bigger then 0, using one thread instead. " << std::endl;
    m_iGlobalMaxThreads = 1;
  }
  else
  {
    m_iGlobalMaxThreads = num;
  }
}

unsigned int icnsImageThreads::SetNumberOfThreads(unsigned int num)
{
  if( num > m_iGlobalMaxThreads )
  {
    std::cout << "\nWarning: Number of threads bigger then global maximum, using global maximum of "
    << IMI_MAX_NUMBER_OF_THREADS << " instead." << std::endl;
    m_iNumberOfThreads = m_iGlobalMaxThreads;
  }
  else if(num <= 0)
  {
    std::cout << "\nWarning: Number of threads must be bigger then 0, using one thread instead." << std::endl;
    m_iNumberOfThreads = 1;
  }
  else
  {
    m_iNumberOfThreads = num;
  }

  return m_iNumberOfThreads;
}

bool icnsImageThreads::SetThreadMethod(
    unsigned int index,
    ImageThreadMethodType method,
    icnsImageThreadParameters *param)
{
  // You can only set the method for 0 through NumberOfThreads-1.
  if ( index >= GetNumberOfThreads() )
  {
    std::cerr << "Can't set thread-method " << index << " with a thread count of "
    << GetNumberOfThreads() << "."  << std::endl;
    return false;
  }
  else
  {
    m_ThreadMethods[index] = method;
    m_ThreadParameter[index] = param;
  }

  return true;
}

icnsImageThreadParameters* icnsImageThreads::GetThreadParameters( unsigned int index )
{
  if( index >= GetNumberOfThreads() )
  {
    std::cout << "\nWarning: Thread not existing, cannot return parameters!" << std::endl;
    return NULL;
  }
  else
  {
    return m_ThreadParameter[index];
  }
}

// -------------------------------------------------------------------
// Include the system depending thread - code
// -------------------------------------------------------------------
#ifdef IMI_USE_PTHREADS
#include "icnsImageThreads_pthread.cxx"
#else
#include "icnsImageThreads_sequ.cxx"
#endif
