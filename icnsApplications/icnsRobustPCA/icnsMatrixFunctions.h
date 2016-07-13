/** \file imiMatrixFunctions.h
 *
 *  \b Initial \b Author: Jan Ehrhardt \n\n
 *  \b Copyright (C) 2014 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiMatrixFunctions_h
#define __imiMatrixFunctions_h

// std includes
#include <string>

// Armadillo includes:
#include <armadillo>

// VNL includes:
#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// ITK includes:
#include "itkConstNeighborhoodIterator.h"

// IMI includes:
//#include "imiMacro.h"
//#include "imiITKImageTypeDefinitions.h"

namespace imi
{
/** \class imi::imiMatrixFunctions
 *  \brief Header file for a class collecting general static functions for matrix processing.
 *
 *  <small> <!--Copyright Information: -->
 *  \b Initial \b Author: Jan Ehrhardt \n
 *  \b Copyright (C) 2014 Institute of Medical Informatics,
 *     University of Luebeck\n
 *  </small>
 ***************************************************************************/
class imiMatrixFunctions
{
public:
  typedef itk::Vector<float, 2> Vector2DPixelType;
  typedef itk::Image<Vector2DPixelType, 2> ImageField2DType;
  typedef ImageField2DType::Pointer        ImageField2DTypePointer;
  
  typedef itk::Vector<float, 3> DisplacementFieldPixelType;
  typedef itk::Image<DisplacementFieldPixelType, 3> DisplacementFieldType;
  typedef DisplacementFieldType::Pointer DisplacementFieldPointerType;
  
  typedef short ImagePixelType;
  typedef itk::Image<ImagePixelType, 3> ImageType;
  typedef ImageType::Pointer ImagePointerType;
  
  typedef unsigned char BinarySegmentationPixelType;
  typedef itk::Image<BinarySegmentationPixelType, 3> BinarySegmentationType;
  typedef BinarySegmentationType::Pointer BinarySegmentationPointerType;
  
  typedef arma::vec ArmaVectorType;
  typedef arma::uvec ArmaIndexVectorType;
  typedef arma::mat ArmaMatrixType;

  typedef double VnlMatrixValueType;
  typedef vnl_vector<VnlMatrixValueType> VnlVectorType;
  typedef vnl_matrix<VnlMatrixValueType> VnlMatrixType;

  static bool WriteMatrixToMatlab(const ArmaMatrixType &M, std::string filename,  std::string matrixName );
  static bool WriteMatrixToMatlab(const VnlMatrixType &M, std::string filename,  std::string matrixName );

  static bool ReadMatrixFromMatlab(ArmaMatrixType &M, std::string filename);
  static bool ReadMatrixFromMatlab(VnlMatrixType &M, std::string filename);

  static bool ConvertArmadilloToVnl(const ArmaMatrixType &inMat, VnlMatrixType &outMat);
  static bool ConvertVnlToArmadillo(const VnlMatrixType &inMat, ArmaMatrixType &outMat);

  // Compute L = Id - D*A*D (with D=diag(column_sum(A)))
  static bool ComputeLaplacianMatrix(const ArmaMatrixType &inMat, ArmaMatrixType &outL);
  // Compute L = Id - D*A*D (with D=diag(column_sum(A))) return L and D
  static bool ComputeLaplacianMatrix(const ArmaMatrixType &inMat, ArmaMatrixType &outL, ArmaVectorType &outD);

  // Compute L = Id - D*A*D (with D=diag(column_sum(A))) return L and D
  //static bool ComputeLaplacianMatrix(const ArmaMatrixType &inMat, ArmaMatrixType &outL, ArmaVectorType &outD);

  static bool ConvertVectorFieldToColumnVector(const DisplacementFieldPointerType vectorField, const BinarySegmentationPointerType mask, unsigned int vectorSize, ArmaVectorType &outVec, bool bUseMagnitudeOnly=false);
  static bool ConvertVectorFieldToColumnVector(const ImageField2DTypePointer vectorField, unsigned int vectorSize, ArmaVectorType &outVec, bool bUseMagnitudeOnly=false);

  static bool ConvertImageToColumnVector(const ImagePointerType image, const BinarySegmentationPointerType mask, ArmaVectorType &outVec);
  static bool ConvertColumnVectorToImage(const ImagePointerType refImage, ImagePointerType& outputImage, const BinarySegmentationPointerType mask, const ArmaVectorType &outVec);


};

} // namespace

#endif // __imiMatrixFunctions_h
