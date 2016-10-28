/** \file imiImageThreads_sequ.cxx
 *
 *  \b Initial \b Author: Alexander Schmidt-Richberg, Jan Ehrhardt\n\n
 *  \b Copyright (C) 2010 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// Project includes:
#include "icnsImageThreads.h"

using namespace imi;

int imiImageThreads::ExecuteThreads()
{
  unsigned int thread_loop;

  // Test if for all threads a function is set.
  for ( thread_loop = 0; thread_loop < GetNumberOfThreads(); thread_loop++ )
  {
    if ( m_ThreadMethods[thread_loop] == (ImageThreadMethodType)NULL )
    {
      std::cerr << "No thread-method set for index: " << thread_loop << std::endl;
      return THR_NO_THREAD_METHOD;
    }
  }

  // Loop over the threads and start the methods sequentially.
  for ( thread_loop = 0; thread_loop < GetNumberOfThreads(); thread_loop++ )
  {
    m_ThreadParameter[thread_loop]->maxNumOfThreads = GetNumberOfThreads();
    m_ThreadParameter[thread_loop]->threadNum = thread_loop;

    m_ThreadMethods[thread_loop]( (void *) m_ThreadParameter[thread_loop] );
  }

  return THR_NOERROR;
}
