/** \file imiImageThreads.h
 *
 *  \b Initial \b Author: Alexander Schmidt-Richberg, Jan Ehrhardt \n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __icnsImageThreads_h
#define __icnsImageThreads_h

// System
#ifdef _WIN32
#include <windows.h>
#include <winbase.h>
#endif

#ifdef __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#include <iostream>

/** Maximal number of threads supported. */
#define IMI_MAX_NUMBER_OF_THREADS 32

/** Default number of threads. */
#define IMI_DEFAULT_THREADNUM 4

/** \brief Base class for all parameters given to the tread calling functions. */
class icnsImageThreadParameters
{
public:
  /** The number of threads started. */
  unsigned int maxNumOfThreads;

  /** the id of this thread. */
  unsigned int threadNum;

  /** The index of the startpixel. */
  unsigned int from;

  /**The index of the endpixel */
  unsigned int to;
};

/** Error code definition. use the test (return value >= THR_FATAL_ERRORS) to
 *  check if you have to discard your data. THR_NOERROR announces that
 *  everything is okay. */
enum ThrErrorCode {
  THR_NOERROR            = 0,   //
  THR_UNDEFINED          = -1,  //
  THR_NO_THREAD_METHOD,         // No method set fot a thread.
  THR_NOT_CREATED,              // Thread has not been created.
  THR_WAIT_ERROR                // An error occurred while waiting for the
                                // end of a thread.
};

// For pthreads the picThreadFunctionType has return value void*,
// otherwise the return type is void.
#ifdef _WIN32
// TODO?
#elif defined(__APPLE__)
  #define IMI_USE_PTHREADS
#elif defined(_SC_NPROCESSORS_ONLN)
  #define IMI_USE_PTHREADS
#elif defined(_SC_NPROC_ONLN)
  #define IMI_USE_PTHREADS
#endif

#ifdef IMI_USE_PTHREADS
  typedef void *(*ImageThreadMethodType)(void *);
  #define IMI_THREAD_RETURN_VALUE  NULL
  #define IMI_THREAD_RETURN_TYPE   void *
#else
  typedef void (*ImageThreadMethodType)(void *);
  #define IMI_THREAD_RETURN_VALUE
  #define IMI_THREAD_RETURN_TYPE   void
#endif

/** \brief A class to handle multithreading.
 *
 *  This class helps to use multiple threads within a program. it is designed to
 *  support parallelization of algorithms but can be used for arbitrary work.
 *
 *  To multithread your program, you have to build an instance of this class,
 *  inherit your own parameter class from picThreadParameter and setup the
 *  thread calling functions and its parameters in SetThreadMethod() (for
 *  all indices from 0 to GetMaxNoThreads() must a function be set !!)
 *  Then call ExecuteThreads().
 *
 *  The Parameter will be casted to void* and its up to you, to do a correct
 *  recast !
 *
 *  If you have defined the PIC_CONFIG_FILE environment variable, you can put
 *  a number into this file, and this number will treated as number of threads.
 *  You can change this number during the program executes, and send a
 *  signal 16 to the application (kill -16 pid). */
class icnsImageThreads
{
public:
  /** \brief Constructor including initializations. */
  icnsImageThreads();

  /** \brief Destructor, kills running threads. */
  ~icnsImageThreads();

  /** \brief Get the global default-number of all threads.
   *
   *  This number is the number of processors or the number in the PIC_CONFIG_FILE. */
  static unsigned int GetGlobalMaxThreads(){return m_iGlobalMaxThreads;}

  /** \brief Set the number of threads globaly. */
  static void SetGlobalMaxThreads(unsigned int num);

  /** \brief Change the number of threads for this object.
   *
   *  By default m_iGlobalMaxThreads is used. Must be smaller or equal to
   *  m_iGlobalMaxThreads. */
  unsigned int SetNumberOfThreads(unsigned int num);

  /** \brief Get the number of threads for this object. */
  unsigned int GetNumberOfThreads() {return m_iNumberOfThreads;};

  /** \brief Setup the functions and parameter.
   *
   *  This function must be called for every thread, indicated by the index. */
  bool SetThreadMethod(
      unsigned int index,
      ImageThreadMethodType method,
      icnsImageThreadParameters *param);

  /** \brief get parameters of a certain thread. */
  icnsImageThreadParameters* GetThreadParameters( unsigned int index );

  /** \brief Runs the threads, after all functions are set up. */
  int ExecuteThreads();

private:
  /** Global max thread number. */
  static unsigned int m_iGlobalMaxThreads;

  /** Maximal number of threads of the instance. */
  unsigned int m_iNumberOfThreads;

  /** Array to store the thread parameters. */
  icnsImageThreadParameters* m_ThreadParameter[IMI_MAX_NUMBER_OF_THREADS];

  /** Array to store pointers to the methods. */
  ImageThreadMethodType m_ThreadMethods[IMI_MAX_NUMBER_OF_THREADS];
};

#endif // ICNSIMAGETHREADS_H_
