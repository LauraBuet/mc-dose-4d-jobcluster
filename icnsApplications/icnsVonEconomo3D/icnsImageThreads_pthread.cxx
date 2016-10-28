/** \file imiImageThreads_pthread.cxx
 *
 *  \b Initial \b Author: Alexander Schmidt-Richberg, Jan Ehrhardt\n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// System includes:
#include <pthread.h>

// Project includes:
#include "icnsImageThreads.h"

//using namespace imi;

int icnsImageThreads::ExecuteThreads()
{
  int ret;
  int thr_error = THR_NOERROR;
  unsigned int thread_loop;
  pthread_t process_id[IMI_MAX_NUMBER_OF_THREADS];
  pthread_attr_t attr;

  pthread_attr_init(&attr);

#ifndef __sgi
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
#endif

  // Test if for all threads a function is set.
  for ( thread_loop = 0; thread_loop < GetNumberOfThreads(); thread_loop++ )
  {
    if ( m_ThreadMethods[thread_loop] == (ImageThreadMethodType)NULL)
    {
      std::cerr << "No thread-method set for index: " << thread_loop << std::endl;
      return THR_NO_THREAD_METHOD;
    }
  }

  // Loop over the threads and start the methods sequentially.
  for ( thread_loop = 1; thread_loop < GetNumberOfThreads(); thread_loop++ )
  {
    m_ThreadParameter[thread_loop]->maxNumOfThreads = GetNumberOfThreads();
    m_ThreadParameter[thread_loop]->threadNum = thread_loop;

    ret = pthread_create(
        &(process_id[thread_loop]),
        &attr,
        m_ThreadMethods[thread_loop],
        (void *) m_ThreadParameter[thread_loop] );

    if(ret != 0)
    {
      std::cerr << "Execution of pthread_create failed !" << std::endl;
      return THR_NOT_CREATED;
    }
  }
  // Start the first thread-method.
  m_ThreadParameter[0]->maxNumOfThreads = GetNumberOfThreads();
  m_ThreadParameter[0]->threadNum = 0;
  m_ThreadMethods[0]((void *)m_ThreadParameter[0]);

  // The parent thread has finished its method - so now it waits for each of the
  // other processes (created with sproc) to exit.
  for ( thread_loop = 1; thread_loop < GetNumberOfThreads(); thread_loop++ )
  {
    int ret = pthread_join(  process_id[thread_loop], NULL );
    if(ret != 0)
    {
      std::cerr << "Error waiting for the termination of a thread (pthread_join)!" << std::endl;
      thr_error = THR_WAIT_ERROR;
    }
  }

  return thr_error;
}
