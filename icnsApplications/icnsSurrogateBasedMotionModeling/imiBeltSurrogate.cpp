/** \file imiBeltSurrogate.cpp
 *
 *  \b Initial \b Author: Matthias Wilms, Maximilian Blendowski \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// Project includes
#include "imiBeltSurrogate.h"

using namespace imi;

// ----------------------------
//   Typedefs:
// ----------------------------


// --------------------------------------------------------
//   Constructor / Destructor / No Copy-Constructor needed
// --------------------------------------------------------

imiBeltSurrogate::imiBeltSurrogate()
{
  m_roiXLow = 0;
  m_roiXHigh = 1;
  m_roiYLow = 0;
  m_roiYHigh = 1;
}

imiBeltSurrogate::~imiBeltSurrogate()
{
}

imiBeltSurrogate* imiBeltSurrogate::New()
{
  return new imiBeltSurrogate();
}

//--------------------------------------------------------
// Methods:
//--------------------------------------------------------

//inherited methods:
bool imiBeltSurrogate::ComputeMeasurements()
{
  m_measurementMatrix = VnlMatrixType( 2, m_inputImages.size() );
  ImageType::IndexType setIndex;

  const float voxelVolume=m_inputImages[0]->GetSpacing()[0]*m_inputImages[0]->GetSpacing()[1]*m_inputImages[0]->GetSpacing()[2];

  long curr_voxelNr = 0;
  long prev_voxelNr = 0;
  for(unsigned int pic = 0; pic < m_inputImages.size(); pic++ )
  {
    //go over phases
    for(unsigned int z_ind = m_beltPosition; z_ind < m_beltPosition + m_beltSize; z_ind++ )
    {
      //go over z-direction for belt voxels
      setIndex[2] = z_ind;
      for( unsigned int x_dir = m_roiXLow; x_dir <= m_roiXHigh; x_dir++ )
      {
        //go over x-direction for belt voxels
        setIndex[0] = x_dir;
        for( unsigned int y_dir = m_roiYLow; y_dir <= m_roiYHigh; y_dir++ )
        {
          //go over y-direction for belt voxels
          setIndex[1] = y_dir;
          if( m_inputImages[pic]->GetPixel( setIndex ) != 0 )
          {
            curr_voxelNr++;
          }
        }
      }
    }
    //scale to mm^3
    curr_voxelNr=curr_voxelNr*voxelVolume;
    //iteration for one phase is complete: filling the regressor matrix
    m_measurementMatrix[0][pic] = curr_voxelNr;
    m_measurementMatrix[1][pic] = curr_voxelNr - prev_voxelNr;
    prev_voxelNr = curr_voxelNr;
    curr_voxelNr = 0;
  }
  m_measurementMatrix[1][0] = m_measurementMatrix[1][0] - m_measurementMatrix[0][m_inputImages.size() - 1];//only reasonable in the case of 4D CT data but not for 4D MRI sequences.
  return true;
}
