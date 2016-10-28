/** \file imiSurrogateBasedMotionPredictionTypeDefinitions.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiSurrogateBasedMotionPredictionTypeDefinitions_h
#define __imiSurrogateBasedMotionPredictionTypeDefinitions_h

#include <vnl/vnl_matlab_filewrite.h>
#include <vnl/vnl_matlab_read.h>

// ITK includes:
#include <itkImageFileReader.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkImageFileWriter.h>
#include <itkVector.h>

// Global typedefs:
const unsigned int Dimension = 3;

typedef float                                             MatrixValueType;
typedef vnl_vector<MatrixValueType>                       VnlVectorType;
typedef vnl_matrix<MatrixValueType>                       VnlMatrixType;
typedef vnl_diag_matrix<MatrixValueType>                  VnlDiagMatrixType;

typedef itk::Image<short, 3>                              ImageType;
typedef ImageType::Pointer                                ImagePointerType;
typedef itk::Vector<float, 3>                             DisplacementFieldPixelType;
typedef itk::Image<DisplacementFieldPixelType, Dimension> DisplacementFieldType;
typedef DisplacementFieldType::Pointer                    DisplacementFieldPointerType;
typedef itk::ImageFileReader<DisplacementFieldType>       FieldReaderType;


//Helpers. Should be placed somewhere else...

namespace icnsSurrogateBasedMotionPredictionHelpers
{

  inline bool ReadMatrixFromMatlabFileDouble( vnl_matrix<MatrixValueType>& vnlMatrix, std::string filename )
  {
    std::cout << "  Reading file " << filename << " ... " << std::flush;

    vcl_ifstream fid;
    fid.open( filename.c_str() );
    vnl_matrix<double> tempDoubleMatrix;
    vnl_matlab_read_or_die( fid, tempDoubleMatrix );
    
    vnlMatrix.set_size( tempDoubleMatrix.rows(), tempDoubleMatrix.cols() );
    for( unsigned int k = 0; k < tempDoubleMatrix.rows(); k++ )
    {
      for( unsigned int l = 0; l < tempDoubleMatrix.cols(); l++ )
      {
        vnlMatrix[k][l] = static_cast<MatrixValueType>( tempDoubleMatrix[k][l] );
      }
    }
    
    std::cout << "Dim: " << vnlMatrix.rows() << "x" << vnlMatrix.cols() << std::endl;
    //vcl_cerr << vnlMatrix;
    return true;
  }

  inline bool ReadMatrixFromMatlabFileRealType( vnl_matrix<MatrixValueType>& vnlMatrix, std::string filename )
  {
    std::cout << "  Reading file " << filename << " ... " << std::flush;

    vcl_ifstream fid;
    fid.open( filename.c_str() );
    vnl_matlab_read_or_die( fid, vnlMatrix );
      
    std::cout << "Dim: " << vnlMatrix.rows() << "x" << vnlMatrix.cols() << std::endl;
      
    return true;
  }

  inline std::string GetFileEnding( std::string filename )
  {
    std::string pattern = ".";
    std::string::size_type pos = filename.rfind( pattern, filename.length() );
    std::string fileTypePattern = filename.substr( pos + 1, 3 );

    return fileTypePattern;
  }

  inline bool GenerateMatrixFromVectorFields( vnl_matrix<MatrixValueType>& outputMatrix, std::vector<std::string> fieldFilenames, std::string maskFilename )
  {
    /* ******************************************************************
     * If a mask file is given (and this is highly recommended for
     * memory efficiency reasons ...), then use this file for masking
     * the relevant regions.
     * ******************************************************************/

    bool useMask = false;
    unsigned int noOfMaskVoxels = 0;
    ImagePointerType maskImage;

    typedef itk::ImageRegionConstIterator<ImageType> ConstSegmentationIteratorType;
    ConstSegmentationIteratorType maskIt;

    if( !maskFilename.empty() )
    {
      typedef itk::ImageFileReader<ImageType> MaskReaderType;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader->SetFileName( maskFilename );
      try
      {
        maskReader->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        std::cerr <<  " Load failed with exception:" << excep << std::endl;
      }
      useMask = true;
      maskImage = maskReader->GetOutput();
      maskImage->Update();

      maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );
      maskIt.GoToBegin();
      while( !maskIt.IsAtEnd() )
      {
        if( maskIt.Get() != 0 ) noOfMaskVoxels++;
        ++maskIt;
      }

      std::cout << "  Using specified mask ..." << std::endl;
      std::cout << "  Number of mask voxels: " << noOfMaskVoxels << std::endl;
    }

    /* ******************************************************************
     * Now coming up for the core functionality:
     * ******************************************************************/

    unsigned int noOfRegressandObservations  = fieldFilenames.size();
    unsigned int noOfRegressandVectorEntries = 0;

    DisplacementFieldType::SizeType fieldSize;
    ImageType::SizeType maskSize;

    for( unsigned int i = 0; i < noOfRegressandObservations; i++ )
    {
      FieldReaderType::Pointer fieldReader = FieldReaderType::New();
      fieldReader->SetFileName( fieldFilenames[i] );
      try
      {
        fieldReader->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        std::cerr << " Load failed with exception:" << excep << std::endl;
        return EXIT_FAILURE;
      }
      DisplacementFieldPointerType inputField = fieldReader->GetOutput();

      // If first field is loaded, allocate vnl regressand matrix:
      if( i == 0 )
      {
        if( useMask )
        {
          noOfRegressandVectorEntries = noOfMaskVoxels * Dimension;
          maskSize = maskImage->GetLargestPossibleRegion().GetSize();
          fieldSize = inputField->GetLargestPossibleRegion().GetSize();
        }
        else
        {
          fieldSize = inputField->GetLargestPossibleRegion().GetSize();
          unsigned int noOfVoxels = 1;
          for( unsigned int d = 0; d < Dimension; d++ )
          {
            noOfVoxels *= fieldSize[d];
          }
          noOfRegressandVectorEntries = noOfVoxels * Dimension;
        }

        // Allocate matrix:
        std::cout << "  Allocate VNL FLOAT matrix (size: " << noOfRegressandVectorEntries << " x " << noOfRegressandObservations << ") ..." << std::endl;
        if( !outputMatrix.set_size( noOfRegressandVectorEntries, noOfRegressandObservations ) )
        {
          std::cout << "  Cannot allocate matrix of requested size !" << std::endl;
          return EXIT_FAILURE;
        }
        outputMatrix.fill( 0.0 );
      }
      // Otherwise check size. Has to be the same for all fields.
      else
      {
        if( inputField->GetLargestPossibleRegion().GetSize() != fieldSize )
        {
          std::cout << "  Loaded fields have different size!" << std::endl;
          return EXIT_FAILURE;
        }
        if( useMask && (inputField->GetLargestPossibleRegion().GetSize() != maskSize) )
        {
          std::cerr << "  Loaded field and mask have different size!" << std::endl;
          return EXIT_FAILURE;
        }
      }

      // Now, write field information into vnl matrix:
      std::cout << "  Converting itk field to vnl matrix ..." << std::endl;

      typedef itk::ImageRegionConstIterator<DisplacementFieldType> ConstFieldIteratorType;
      ConstFieldIteratorType inputIt( inputField, inputField->GetLargestPossibleRegion() );
      inputIt.GoToBegin();
      unsigned int vector_pos = 0;
      DisplacementFieldType::PixelType vectorValue;

      // First version (= default version): No mask used ...
      if( !useMask )
      {
        while( !inputIt.IsAtEnd() )
        {
          vectorValue = inputIt.Get();
          for( unsigned int d = 0; d < Dimension; d++ )
          {
            outputMatrix[vector_pos][i] = vectorValue[d];
            vector_pos++;
          }
          ++inputIt;
        }
      }
        
      // Second version: Use mask ...
      else
      {
        maskIt.GoToBegin();
        while( !inputIt.IsAtEnd() )
        {
          if( maskIt.Get() != 0.0 )
          {
            vectorValue = inputIt.Get();
            for( unsigned int d = 0; d < Dimension; d++ )
            {
              outputMatrix[vector_pos][i] = vectorValue[d];
              vector_pos++;
            }
          }
          ++inputIt;
          ++maskIt;
        }
      }

      // Final check:
      if( vector_pos != outputMatrix.rows() )
      {
        std::cout << "  Matrix not correctly filled with field values!" << std::endl;
        return false;
      }
    } // End for-loop over fields.
    return true;
  }

  inline bool GenerateVectorFieldFromVnlVector( vnl_matrix<MatrixValueType>& outputAsVnlVector, DisplacementFieldPointerType& outputField, std::string refFieldname, std::string maskFilename )
  {
    //---------------------------------------------------
    // Load reference field which serves as a template.
    //---------------------------------------------------

    std::cout << "  Loading template field ..." << std::endl;

    typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
    FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( refFieldname.c_str() );
    try
    {
      fieldReader->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << " Load failed with exception:" << excep << std::endl;
      return EXIT_FAILURE;
    }
    outputField = fieldReader->GetOutput();
    outputField->DisconnectPipeline();
    DisplacementFieldType::SizeType fieldSize = outputField->GetLargestPossibleRegion().GetSize();

    //---------------------------------------------------
    // If specified: load mask for conversion:
    //---------------------------------------------------

    bool useMask = false;
    unsigned int noOfMaskVoxels = 0;
    ImagePointerType maskImage;

    typedef itk::ImageRegionConstIterator<ImageType> ConstSegmentationIteratorType;
    ConstSegmentationIteratorType maskIt;

    if( !maskFilename.empty() )
    {
      typedef itk::ImageFileReader<ImageType> MaskReaderType;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader->SetFileName( maskFilename );
      try
      {
        maskReader->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        std::cerr << "  Load failed with exception:" << excep << std::endl;
      }
      useMask = true;
      maskImage = maskReader->GetOutput();
      maskImage->Update();

      maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );
      maskIt.GoToBegin();
      while( !maskIt.IsAtEnd() )
      {
        if( maskIt.Get() != 0 )
          noOfMaskVoxels++;
        ++maskIt;
      }

      std::cerr << "  Using specified mask for conversion ..." << std::endl;
      std::cout << "\n  Number of lung voxels: " << noOfMaskVoxels << std::endl;
    }

    //---------------------------------------------------
    // Convert vnl matrix to mha vector.
    //---------------------------------------------------

    std::cerr << "  Converting vnl to mha ..." << std::endl;

    typedef itk::ImageRegionIterator<DisplacementFieldType> FieldIteratorType;
    FieldIteratorType outputIt( outputField, outputField->GetLargestPossibleRegion() );
    outputIt.GoToBegin();

    unsigned int urow = 0;
    DisplacementFieldType::PixelType vectorValue;

    // Default version: no lung mask used to fill the field:
    if( !useMask )
    {
      while( !outputIt.IsAtEnd() )
      {
        vectorValue[0] = outputAsVnlVector[urow][0];
        vectorValue[1] = outputAsVnlVector[urow + 1][0];
        vectorValue[2] = outputAsVnlVector[urow + 2][0];
        urow += 3;
        outputIt.Set( vectorValue );
        ++outputIt;
      }
    }

    // Alternative version: lung mask used to fill the field:
    // If mask == 0, fill (0,0,0) into field. Otherwise
    // fill in vnl vector entries.
    else
    {
      maskIt.GoToBegin();
      while( !outputIt.IsAtEnd() )
      {
        if( maskIt.Get() != 0 )
        {
          vectorValue[0] = outputAsVnlVector[urow][0];
          vectorValue[1] = outputAsVnlVector[urow + 1][0];
          vectorValue[2] = outputAsVnlVector[urow + 2][0];
          urow += 3;
          outputIt.Set( vectorValue );
        }
        else
        {
          vectorValue[0] = 0;
          vectorValue[1] = 0;
          vectorValue[2] = 0;
          outputIt.Set( vectorValue );
        }
        ++outputIt;
        ++maskIt;
      }
    }
    outputField->Modified();
    return true;
  }


  inline bool SaveVnlVectorAsVectorField( vnl_matrix<MatrixValueType>& outputAsVnlVector, std::string outFieldFilename, std::string refFieldname, std::string maskFilename )
  {
    //---------------------------------------------------
    // Load reference field which serves as a template.
    //---------------------------------------------------

    std::cout << "  Loading template field ..." << std::endl;

    typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
    FieldReaderType::Pointer fieldReader = FieldReaderType::New();
    fieldReader->SetFileName( refFieldname.c_str() );
    try
    {
      fieldReader->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr << " Load failed with exception:" << excep << std::endl;
      return EXIT_FAILURE;
    }
    DisplacementFieldPointerType outputField = fieldReader->GetOutput();
    DisplacementFieldType::SizeType fieldSize = outputField->GetLargestPossibleRegion().GetSize();

    //---------------------------------------------------
    // If specified: load mask for conversion:
    //---------------------------------------------------

    bool useMask = false;
    unsigned int noOfMaskVoxels = 0;
    ImagePointerType maskImage;

    typedef itk::ImageRegionConstIterator<ImageType> ConstSegmentationIteratorType;
    ConstSegmentationIteratorType maskIt;

    if( !maskFilename.empty() )
    {
      typedef itk::ImageFileReader<ImageType> MaskReaderType;
      MaskReaderType::Pointer maskReader = MaskReaderType::New();
      maskReader->SetFileName( maskFilename );
      try
      {
        maskReader->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        std::cerr << "  Load failed with exception:" << excep << std::endl;
      }
      useMask = true;
      maskImage = maskReader->GetOutput();
      maskImage->Update();

      maskIt = ConstSegmentationIteratorType( maskImage, maskImage->GetLargestPossibleRegion() );
      maskIt.GoToBegin();
      while( !maskIt.IsAtEnd() )
      {
        if( maskIt.Get() != 0 )
          noOfMaskVoxels++;
        ++maskIt;
      }

      std::cout << "  Using specified mask for conversion ..." << std::endl;
      std::cout << "  Number of lung voxels: " << noOfMaskVoxels << std::endl;
    }

    //---------------------------------------------------
    // Convert vnl matrix to mha vector.
    //---------------------------------------------------

    std::cout <<  "  Converting vnl to mha ..." << std::endl;

    typedef itk::ImageRegionIterator<DisplacementFieldType> FieldIteratorType;
    FieldIteratorType outputIt( outputField, outputField->GetLargestPossibleRegion() );
    outputIt.GoToBegin();

    unsigned int urow = 0;
    DisplacementFieldType::PixelType vectorValue;

    // Default version: no lung mask used to fill the field:
    if( !useMask )
    {
      while( !outputIt.IsAtEnd() )
      {
        vectorValue[0] = outputAsVnlVector[urow][0];
        vectorValue[1] = outputAsVnlVector[urow + 1][0];
        vectorValue[2] = outputAsVnlVector[urow + 2][0];
        urow += 3;
        outputIt.Set( vectorValue );
        ++outputIt;
      }
    }
    
    // Alternative version: lung mask used to fill the field:
    // If mask == 0, fill (0,0,0) into field. Otherwise
    // fill in vnl vector entries.
    else
    {
      maskIt.GoToBegin();
      while( !outputIt.IsAtEnd() )
      {
        if( maskIt.Get() != 0 )
        {
          vectorValue[0] = outputAsVnlVector[urow][0];
          vectorValue[1] = outputAsVnlVector[urow + 1][0];
          vectorValue[2] = outputAsVnlVector[urow + 2][0];
          urow += 3;
          outputIt.Set( vectorValue );
        }
        else
        {
          vectorValue[0] = 0;
          vectorValue[1] = 0;
          vectorValue[2] = 0;
          outputIt.Set( vectorValue );
        }
        ++outputIt;
        ++maskIt;
      }
    }
    outputField->Modified();

    //---------------------------------------------------
    // Save output field as specified.
    //---------------------------------------------------

    std::cout << "  Saving mha ..." << std::endl;

    typedef itk::ImageFileWriter<DisplacementFieldType> FieldWriterType;
    FieldWriterType::Pointer displacementFieldWriter = FieldWriterType::New();

    displacementFieldWriter->SetInput( outputField );
    displacementFieldWriter->SetFileName( outFieldFilename.c_str() );
    try
    {
      displacementFieldWriter->Update();
    }
    catch( itk::ExceptionObject & excep )
    {
      std::cerr <<  " Save failed with exception:" << excep << std::endl;
      return false;
    }

    return true;
  }

} // end of namespace

#endif
