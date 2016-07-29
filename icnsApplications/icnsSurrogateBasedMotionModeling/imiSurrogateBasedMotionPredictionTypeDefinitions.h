/** \file imiSurrogateBasedMotionPredictionTypeDefinitions.h
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2012 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#ifndef __imiSurrogateBasedMotionPredictionTypeDefinitions_h
#define __imiSurrogateBasedMotionPredictionTypeDefinitions_h

#include "imiMacro.h"
#include "imiITKImageTypeDefinitions.h"

#include "vnl/vnl_matlab_filewrite.h"
#include "vnl/vnl_matlab_read.h"

// ITK includes:
#include "itkImageFileReader.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkImageFileWriter.h"

///////////////////////////////////////////////////
//
// Matrix and vector type definitions
//

typedef float imiRealType2;
typedef imiRealType2 MatrixValueType;
typedef vnl_vector<MatrixValueType> VnlVectorType;
typedef vnl_matrix<MatrixValueType> VnlMatrixType;
typedef vnl_diag_matrix<MatrixValueType> VnlDiagMatrixType;


//Helpers. Should be placed somewhere else...
namespace imi
{

  namespace imiSurrogateBasedMotionPredictionHelpers
  {

    inline bool ReadMatrixFromMatlabFileDouble( vnl_matrix<MatrixValueType>& vnlMatrix, std::string filename )
    {
      imiINFO( "  Reading file " << filename << " ...\n" );

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
      //vcl_cerr << vnlMatrix;
      return true;
    }

    inline bool ReadMatrixFromMatlabFileRealType( vnl_matrix<MatrixValueType>& vnlMatrix, std::string filename )
    {
      imiINFO( "  Reading file " << filename << " ...\n" );

      vcl_ifstream fid;
      fid.open( filename.c_str() );
      vnl_matlab_read_or_die( fid, vnlMatrix );
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
          imiERROR( " Load failed with exception:" << excep << "\n" );
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

        imiINFO( "Using specified mask..." );
        std::cout << "\n  Number of lung voxels: " << noOfMaskVoxels << std::endl;
      }

      /* ******************************************************************
       * Now coming up for the core functionality:
       * ******************************************************************/

      typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;

      unsigned int noOfRegressandObservations = fieldFilenames.size();
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
          imiERROR( " Load failed with exception:" << excep << "\n" );
          return -1;
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
          imiINFO( "  Allocate VNL matrix (size: " << noOfRegressandVectorEntries << " x " << noOfRegressandObservations << ") ..." );
          if( !outputMatrix.set_size( noOfRegressandVectorEntries, noOfRegressandObservations ) )
          {
            imiERROR( "  Cannot allocate matrix of requested size !" );
            return -1;
          }
          outputMatrix.fill( 0.0 );
        }
        // Otherwise check size. Has to be the same for all fields.
        else
        {
          if( inputField->GetLargestPossibleRegion().GetSize() != fieldSize )
          {
            imiERROR( "  Loaded fields have different size!" );
            return -1;
          }
          if( useMask && (inputField->GetLargestPossibleRegion().GetSize() != maskSize) )
          {
            imiERROR( "  Loaded field and mask have different size!" );
            return -1;
          }
        }

        // Now, write field information into vnl matrix:
        imiINFO( "  Converting itk field to vnl matrix ..." );

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
          imiERROR( "  Matrix not correctly filled with field values!" );
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

      imiINFO( "  Loading template field ..." );

      typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
      FieldReaderType::Pointer fieldReader = FieldReaderType::New();
      fieldReader->SetFileName( refFieldname.c_str() );
      try
      {
        fieldReader->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        imiERROR( " Load failed with exception:" << excep << "\n" );
        return -1;
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
          imiERROR( "  Load failed with exception:" << excep << "\n" );
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

        imiINFO( "  Using specified mask for conversion ..." );
        std::cout << "\n  Number of lung voxels: " << noOfMaskVoxels << std::endl;
      }

      //---------------------------------------------------
      // Convert vnl matrix to mha vector.
      //---------------------------------------------------

      imiINFO( "  Converting vnl to mha ..." );

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

      imiINFO( "  Loading template field ..." );

      typedef itk::ImageFileReader<DisplacementFieldType> FieldReaderType;
      FieldReaderType::Pointer fieldReader = FieldReaderType::New();
      fieldReader->SetFileName( refFieldname.c_str() );
      try
      {
        fieldReader->Update();
      }
      catch( itk::ExceptionObject & excep )
      {
        imiERROR( " Load failed with exception:" << excep << "\n" );
        return -1;
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
          imiERROR( "  Load failed with exception:" << excep << "\n" );
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

        imiINFO( "  Using specified mask for conversion ..." );
        std::cout << "\n  Number of lung voxels: " << noOfMaskVoxels << std::endl;
      }

      //---------------------------------------------------
      // Convert vnl matrix to mha vector.
      //---------------------------------------------------

      imiINFO( "  Converting vnl to mha ..." );

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

      imiINFO( "  Saving mha ..." );

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
        imiERROR( " Save failed with exception:" << excep << "\n" );
        return false;
      }

      return true;
    }


  }

}

#endif