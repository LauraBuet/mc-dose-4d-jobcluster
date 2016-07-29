/** \file imiDiaphragmSurrogate.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// Project includes
#include "imiDiaphragmSurrogate.h"

using namespace imi;

// ----------------------------
//   Typedefs:
// ----------------------------

// --------------------------------------------------------
//   Constructor / Destructor / No Copy-Constructor needed
// --------------------------------------------------------

imiDiaphragmSurrogate::imiDiaphragmSurrogate()
{
  m_threshold = 650;
  m_headFirst = false;
}

imiDiaphragmSurrogate::~imiDiaphragmSurrogate()
{
}

imiDiaphragmSurrogate* imiDiaphragmSurrogate::New()
{
  return new imiDiaphragmSurrogate();
}

//--------------------------------------------------------
// Methods:
//--------------------------------------------------------

//inherited methods:
bool imiDiaphragmSurrogate::ComputeMeasurements()
{
  m_measurementMatrix = VnlMatrixType( 1, m_inputImages.size() );
  ImageType::IndexType setIndex;

  const float zSpacing = m_inputImages[0]->GetSpacing()[2];

  MatrixValueType posSum = 0.;

  int currVoxel = 0;
  int prevVoxel = 0;
  setIndex[1] = m_coroSlice;
  int zIncrement;

  if( m_headFirst )
  {
    zIncrement = -1;
  }
  else
  {
    zIncrement = 1;
  }

  //go over phases
  for( unsigned int pic = 0; pic < m_inputImages.size(); pic++ )
  {
    posSum = 0.;
    for( unsigned int starts = 0; starts < m_xStart.size(); starts++ )
    {
      //go over x-direction
      for( unsigned int x = m_xStart[starts]; x < m_xStart[starts] + m_xSize; x++ )
      {
        //std::cout<<"New: "<<std::endl;
        setIndex[0] = x;
        //go over z-direction (bottom -> top; assumption: z=0 is anatomically located under the dome of the diaphragm )
        bool found = false;
        for( int z = m_zStart[starts]; z < m_inputImages[0]->GetLargestPossibleRegion().GetSize( 2 ) && z >= 0 && found == false; z += zIncrement )
        {
          setIndex[2] = z;
          //std::cout<<"z: "<<z<<std::endl;
          currVoxel = m_inputImages[pic]->GetPixel( setIndex );
          if( currVoxel <= m_threshold && ((m_maskImages.size() >0 && m_maskImages[pic]->GetPixel( setIndex ) > 0) || m_maskImages.size() ==0) )
          {
            MatrixValueType pos = z - zIncrement*(static_cast<float>( m_threshold ) - currVoxel) / (prevVoxel-currVoxel);
            //std::cout<<"pos: "<<pos<<std::endl;
            posSum += pos * zSpacing;
            found = true;
          }
          else
          {
            prevVoxel = currVoxel;
          }

        }
      }
    }
    m_measurementMatrix[0][pic] = posSum / (m_xSize * m_xStart.size());
  }
  return true;
}
