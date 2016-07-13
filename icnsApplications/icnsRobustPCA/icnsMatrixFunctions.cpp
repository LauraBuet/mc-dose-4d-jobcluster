/** \file imiMatrixFunctions.cxx
 *
 *  \b Initial \b Author: Jan Ehrhardt \n\n
 *  \b Copyright (C) 2014 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// std includes
#include <string>

// Armadillo includes:

// VNL includes:
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// ITK includes:
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageDuplicator.h"

// Project includes:
// #include "imiMacro.h"
#include "icnsMatrixFunctions.h"
//#include "imiImageFunctions.h"

using namespace imi;

////////////////////////////////////////////////////////////////
//
// WriteMatrixToMatlab( ArmaMatrixType )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::WriteMatrixToMatlab( const ArmaMatrixType &M, std::string filename, std::string matrixName )
{
  VnlMatrixType vnlMat;
  if( !ConvertArmadilloToVnl( M, vnlMat ) )
  {
    std::cerr << "imiMatrixFunctions::WriteMatrixToMatlab(): can not convert armadillo matrix to VNL type!" << std::endl;
    return false;
  }

  return WriteMatrixToMatlab( vnlMat, filename, matrixName );
}

////////////////////////////////////////////////////////////////
//
// WriteMatrixToMatlab( VnlMatrixType )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::WriteMatrixToMatlab( const VnlMatrixType &M, std::string filename, std::string matrixName )
{
  std::cout << "  Write matrix to file " << filename << " ... " << std::endl;

  vnl_matlab_filewrite matrixWrite( filename.c_str() );
  matrixWrite.write( M, matrixName.c_str() );

  return true;
}
////////////////////////////////////////////////////////////////
//
// ReadMatrixFromMatlab( ArmaMatrixType )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ReadMatrixFromMatlab( ArmaMatrixType &M, std::string filename )
{
  VnlMatrixType vnlMat;

  if( !ReadMatrixFromMatlab( vnlMat, filename ) )
  {
    std::cerr << "imiMatrixFunctions::ReadMatrixFromMatlab(): can not read matlab file!" << std::endl;
    return false;
  }

  return ConvertVnlToArmadillo( vnlMat, M );
}

////////////////////////////////////////////////////////////////
//
// ReadMatrixFromMatlab( VnlMatrixType )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ReadMatrixFromMatlab( VnlMatrixType &M, std::string filename )
{
  std::cout << "  Reading file " << filename << " ... " << std::endl;

  vcl_ifstream fid;
  fid.open( filename.c_str() );
  vnl_matlab_read_or_die( fid, M );

  std::cout << "  Completed. \n" << std::endl;

  return true;
}

////////////////////////////////////////////////////////////////
//
// ConvertArmadilloToVnl( )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ConvertArmadilloToVnl( const ArmaMatrixType &inMat, VnlMatrixType &outMat )
{
  outMat.set_size( inMat.n_rows, inMat.n_cols );

  for( unsigned int i = 0; i < outMat.rows(); i++ )
    for( unsigned int j = 0; j < outMat.cols(); j++ )
      outMat( i, j ) = inMat( i, j );

  return true;
}

////////////////////////////////////////////////////////////////
//
// ConvertVnlToArmadillo( )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ConvertVnlToArmadillo( const VnlMatrixType &inMat, ArmaMatrixType &outMat )
{
  outMat.set_size( inMat.rows(), inMat.cols() );

  for( unsigned int i = 0; i < outMat.n_rows; i++ )
    for( unsigned int j = 0; j < outMat.n_cols; j++ )
      outMat( i, j ) = inMat( i, j );

  return true;
}

////////////////////////////////////////////////////////////////
//
// ComuteLaplacianMatrix( )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ComputeLaplacianMatrix( const ArmaMatrixType &inMat, ArmaMatrixType &outL )
{
  ArmaVectorType tempD;

  return ComputeLaplacianMatrix( inMat, outL, tempD );
}
////////////////////////////////////////////////////////////////
//
// ComuteLaplacianMatrix( )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ComputeLaplacianMatrix( const ArmaMatrixType &inMat, ArmaMatrixType &outL, ArmaVectorType &outD )
{
  // Check if inMat is quadratic!
  const unsigned int n = inMat.n_cols;
  if( inMat.n_rows != n )
  {
    std::cerr << "imiMatrixFunctions::ComuteLaplacianMatrix(): input matrix is not quadratic!" << std::endl;
    return false;
  }
  //
  // Compute diagonal matrix D with col-sums of A in diagonal
  //
  const double epsilon = 1.0e-8;
  outD = arma::sum( inMat, 1 );

  for( unsigned int i = 0; i < outD.n_elem; i++ )
    outD( i ) = 1.0 / std::sqrt( outD( i ) ) + epsilon;

  //
  // Compute Laplacian matrix L
  //
  outL = arma::eye( n, n ) - arma::diagmat( outD ) * inMat * arma::diagmat( outD );

  return true;
}

////////////////////////////////////////////////////////////////
//
// ConvertVectorFieldToColumnVector( )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ConvertVectorFieldToColumnVector( const DisplacementFieldPointerType vectorField, const BinarySegmentationPointerType maskImage, unsigned int vectorSize,
    ArmaVectorType &outVec, bool bUseMagnitudeOnly )
{
  typedef itk::ImageRegionConstIterator<BinarySegmentationType> ConstSegmentationIteratorType;

  //
  // Check Input
  //
  if( vectorField.IsNull() )
  {
    std::cerr << "icnsMatrixFunctions::ConvertVectorFieldToColumnVector(): vectorField is NULL!" << std::endl;
    return false;
  }

  const unsigned int numComponents = vectorField->GetNumberOfComponentsPerPixel();
  const bool bUseMask = maskImage.IsNotNull();

  if( bUseMask )
  {
    // Check if mask and vector field have the same size
    if( vectorField->GetRequestedRegion().GetSize() != maskImage->GetRequestedRegion().GetSize() )
    {
      std::cerr << "icnsMatrixFunctions::ConvertVectorFieldToColumnVector(): vectorField and mask must have same size!" << std::endl;
      return false;
    }
  }

  // if length is not given, pre-compute length of output vector
  if( vectorSize == 0 )
  {
    if( !bUseMask )
    {
      // if no mask is used, compute number of image pixels * number of components per pixel
      DisplacementFieldType::SizeType size = vectorField->GetLargestPossibleRegion().GetSize();
      vectorSize = 1;
      for( unsigned int dim = 0; dim < vectorField->GetImageDimension(); dim++ )
        vectorSize *= size[dim];

      if( !bUseMagnitudeOnly )
        vectorSize *= numComponents;
    }
    else
    {
      ConstSegmentationIteratorType maskIt;

      maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );

      unsigned int noOfMaskVoxels = 0;
      maskIt.GoToBegin();
      while( !maskIt.IsAtEnd() )
      {
        if( maskIt.Get() != 0 )
          noOfMaskVoxels++;
        ++maskIt;
      }

      if( bUseMagnitudeOnly )
        vectorSize = noOfMaskVoxels * numComponents;
      else
        vectorSize = noOfMaskVoxels * numComponents;

    }
  }

  //imiDEBUGINFO( 9, "ConvertVectorFieldToColumnVector(): generate vector of length "<<vectorSize<<((bUseMask) ? "with" : "without")<<" using a mask." );

  if( outVec.size() != vectorSize )
  {
    outVec.resize( vectorSize );
  }

  //
  // Now, write field information into vector:
  //

  typedef itk::ImageRegionConstIterator<DisplacementFieldType> ConstFieldIteratorType;

  ConstFieldIteratorType inputIt( vectorField, vectorField->GetLargestPossibleRegion() );
  inputIt.GoToBegin();
  unsigned int vector_pos = 0;
  DisplacementFieldType::PixelType vectorValue;
  if( !bUseMask )
  {
    while( !inputIt.IsAtEnd() )
    {
      // if no mask is used go through all pixels

      vectorValue = inputIt.Get();

      if( bUseMagnitudeOnly )
      {
        outVec( vector_pos ) = vectorValue.GetNorm();
        vector_pos++;
      }
      else
      {
        for( unsigned int d = 0; d < vectorField->GetImageDimension(); d++ )
        {
          outVec( vector_pos ) = vectorValue[d];
          vector_pos++;
        }
      }
      ++inputIt;
    }
  }
  else
  {

    // if mask is used go through all pixels with mask value != 0

    ConstSegmentationIteratorType maskIt;
    maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );
    maskIt.GoToBegin();

    maskIt.GoToBegin();
    while( !inputIt.IsAtEnd() )
    {
      if( maskIt.Get() != 0.0 )
      {
        vectorValue = inputIt.Get();
        for( unsigned int d = 0; d < maskImage->GetImageDimension(); d++ )
        {
          outVec[vector_pos] = vectorValue[d];
          vector_pos++;
        }
      }
      ++inputIt;
      ++maskIt;
    }
  }

  // Final check:
  if( vector_pos != outVec.size() )
  {
    std::cerr << "  Matrix not correctly filled with field values!" << std::endl;
    return false;
  }

  return true;
}
////////////////////////////////////////////////////////////////
//
// VectorizeVectorFieldPatch( )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ConvertVectorFieldToColumnVector( const ImageField2DTypePointer vectorField, unsigned int vectorSize, ArmaVectorType &outVec, bool bUseMagnitudeOnly )
{
  //
  // Check Input
  //
  if( vectorField.IsNull() )
  {
    std::cerr << "imiMatrixFunctions::ConvertVectorFieldToColumnVector(): vectorField is NULL!" << std::endl;
    return false;
  }

  const unsigned int numComponents = vectorField->GetNumberOfComponentsPerPixel();

  // if length is not given, pre-compute length of output vector
  if( vectorSize == 0 )
  {
    // if no mask is used, compute number of image pixels * number of components per pixel
    ImageField2DType::SizeType size = vectorField->GetLargestPossibleRegion().GetSize();
    vectorSize = 1;
    for( unsigned int dim = 0; dim < vectorField->GetImageDimension(); dim++ )
      vectorSize *= size[dim];

    if( !bUseMagnitudeOnly )
      vectorSize *= numComponents;
  }

  //imiDEBUGINFO( 9, "ConvertVectorFieldToColumnVector(): generate vector of length "<<vectorSize<<" without using a mask." );

  if( outVec.size() != vectorSize )
  {
    outVec.resize( vectorSize );
  }

  //
  // Now, write field information into vector:
  //

  typedef itk::ImageRegionConstIterator<ImageField2DType> ConstFieldIteratorType;

  ConstFieldIteratorType inputIt( vectorField, vectorField->GetLargestPossibleRegion() );
  inputIt.GoToBegin();
  unsigned int vector_pos = 0;
  ImageField2DType::PixelType vectorValue;

  while( !inputIt.IsAtEnd() )
  {
    // if no mask is used go through all pixels

    vectorValue = inputIt.Get();
    if( bUseMagnitudeOnly )
    {
      outVec( vector_pos ) = vectorValue.GetNorm();
      vector_pos++;
    }
    else
    {
      for( unsigned int d = 0; d < vectorField->GetImageDimension(); d++ )
      {
        outVec( vector_pos ) = vectorValue[d];
        vector_pos++;
      }
    }
    ++inputIt;
  }

  // Final check:
  if( vector_pos != outVec.size() )
  {
    std::cerr << "  Matrix not correctly filled with field values!" << std::endl;
    return false;
  }

  return true;
}

////////////////////////////////////////////////////////////////
//
// ConvertImageToColumnVector( )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ConvertImageToColumnVector( const ImagePointerType image, const BinarySegmentationPointerType mask, ArmaVectorType &outVec )
{
  typedef itk::ImageRegionConstIterator<BinarySegmentationType> ConstSegmentationIteratorType;

  //
  // Check Input
  //
  if( image.IsNull() )
  {
    std::cerr << "imiMatrixFunctions::ConvertImageToColumnVector(): image is NULL!" << std::endl;
    return false;
  }

  const bool bUseMask = mask.IsNotNull();

  if( bUseMask )
  {
    // Check if mask and vector field have the same size
    if( image->GetRequestedRegion().GetSize() != mask->GetRequestedRegion().GetSize() )
    {
      std::cerr << "icnsMatrixFunctions::ConvertImageToColumnVector(): image and mask must have same size!" << std::endl;
      return false;
    }
  }


  //
  // Now, write image information into vector:
  //

  typedef itk::ImageRegionConstIterator<ImageType> ConstImageIteratorType;

  ConstImageIteratorType inputIt( image, image->GetLargestPossibleRegion() );
  inputIt.GoToBegin();
  unsigned int vector_pos = 0;
  ImageType::PixelType value;
  if( !bUseMask )
  {
    while( !inputIt.IsAtEnd() )
    {
      // if no mask is used go through all pixels

      value = inputIt.Get();
      outVec( vector_pos ) = value;
      vector_pos++;

      ++inputIt;
    }
  }
  else
  {

    // if mask is used go through all pixels with mask value != 0

    ConstSegmentationIteratorType maskIt;
    maskIt = ConstSegmentationIteratorType( mask, mask->GetLargestPossibleRegion() );
    maskIt.GoToBegin();

    maskIt.GoToBegin();
    while( !inputIt.IsAtEnd() )
    {
      if( maskIt.Get() != 0.0 )
      {
        value = inputIt.Get();
        outVec( vector_pos ) = value;
        vector_pos++;
      }
      ++inputIt;
      ++maskIt;
    }
  }

  // Final check:
  if( vector_pos != outVec.size() )
  {
    std::cout << vector_pos << std::endl;
    std::cout << outVec.size() << std::endl;
    std::cerr << "  Vector not correctly filled with field values!" << std::endl;
    return false;
  }


  return true;
}

////////////////////////////////////////////////////////////////
//
// ConvertColumnVectorToImage( )
//
////////////////////////////////////////////////////////////////
bool imiMatrixFunctions::ConvertColumnVectorToImage( const ImagePointerType refImage, ImagePointerType &outputImage, const BinarySegmentationPointerType mask, const ArmaVectorType &outVec )
{
  //allocate output image
  typedef itk::ImageDuplicator< ImageType > DuplicatorType;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(refImage);
  duplicator->Update();
  outputImage = duplicator->GetOutput();

  const bool bUseMask = mask.IsNotNull();

  if( bUseMask )
  {
    // Check if mask and vector field have the same size
    if( outputImage->GetRequestedRegion().GetSize() != mask->GetRequestedRegion().GetSize() )
    {
      std::cerr << "icnsMatrixFunctions::ConvertImageToColumnVector(): image and mask must have same size!" << std::endl;
      return false;
    }
  }

  typedef itk::ImageRegionIterator<ImageType> ImageIteratorType;
  typedef itk::ImageRegionIterator<BinarySegmentationType> ConstSegmentationIteratorType;

  ImageIteratorType inputIt( outputImage, outputImage->GetLargestPossibleRegion() );
  inputIt.GoToBegin();
  unsigned int vector_pos = 0;
  ImageType::PixelType value;

  if( !bUseMask )
  {
    while( !inputIt.IsAtEnd() )
    {
      // if no mask is used go through all pixels

      value = outVec( vector_pos );
      inputIt.Set(value);
      vector_pos++;

      ++inputIt;
    }
  }
  else
  {

    // if mask is used go through all pixels with mask value != 0

    ConstSegmentationIteratorType maskIt;
    maskIt = ConstSegmentationIteratorType( mask, mask->GetLargestPossibleRegion() );
    maskIt.GoToBegin();

    maskIt.GoToBegin();
    while( !inputIt.IsAtEnd() )
    {
      if( maskIt.Get() != 0.0 )
      {
        value = outVec( vector_pos );
        inputIt.Set(value);
        vector_pos++;
      } else
      {
        inputIt.Set(0);
      }
      ++inputIt;
      ++maskIt;
    }
  }

  // Final check:
  if( vector_pos != outVec.size() )
  {
    std::cerr << "  Vector not correctly filled with field values!" << std::endl;
    return false;
  }


  return true;
}
