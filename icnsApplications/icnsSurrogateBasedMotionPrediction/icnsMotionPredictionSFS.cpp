/** \file imiMotionPredictionSFSMethods.cpp
 *
 *  \b Initial \b Author: Matthias Wilms \n\n
 *  \b Copyright (C) 2013 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

#include "icnsMotionPredictionSFS.h"

//#include "imiLinearKernel.h"

#include <limits>

#include <vnl/vnl_tag.h>

namespace imi
{

imiMotionPredictionSFS* imiMotionPredictionSFS::New()
{
  return new imiMotionPredictionSFS();
}

imiMotionPredictionSFS::imiMotionPredictionSFS()
{
  m_predictor = NULL;
  m_ridgeParameterMin = 1;
  m_ridgeParameterMax = 10;
  m_ridgeMultiplicator = 10;
  m_SFSPerformed = false;
  m_regressandMatrix = NULL;
  m_maxNumOfIterations = std::numeric_limits<unsigned int>::max();
  m_optimalSSD = std::numeric_limits<double>::max();
  m_predictors.clear();
}

imiMotionPredictionSFS::~imiMotionPredictionSFS()
{
  //m_predictor->Delete(); //TODO
}

void imiMotionPredictionSFS::PerformSFS()
{
  if( m_predictor )
  {
    //m_predictor->Delete(); // TODO
  }

  m_SFSPerformed=true;

  //TODO
  /*m_predictor = imiMLRMotionPredictionMethods::New();
  m_predictor->SetRegressandTrainingMatrix( *m_regressandMatrix );*/
  unsigned int useEveryNthCol=3;

  VnlMatrixType newRegressand;
  newRegressand.set_size(m_regressandMatrix->rows(),m_regressandMatrix->cols()-1);
  for (unsigned int col=0; col<m_regressandMatrix->cols(); col+=useEveryNthCol)
  {
    icnsMLRMotionPrediction* predictor = icnsMLRMotionPrediction::New();
    unsigned int index=0;
    for (unsigned int j=0; j<m_regressandMatrix->cols();j++)
    {
      if(col!=j)
      {
        newRegressand.set_column(index,m_regressandMatrix->get_column(j));
        index++;
      }
    }
    predictor->SetRegressandTrainingMatrix(newRegressand);
    m_predictors.push_back(predictor);
  }
  

  m_optimalIndices.clear();

  std::set<unsigned int> remainingIndices;

  m_optimalSSD = std::numeric_limits<double>::max();

  for( unsigned int i = 0; i < m_regressorMatrices.size(); i++ )
  {
    remainingIndices.insert( i );
  }

  unsigned int numOfIterations = 0;
  bool ssdImproved = true;
  VnlMatrixType* bestSubset = NULL;
  double bestSubsetSSDValue = std::numeric_limits<double>::max();
  float bestSubsetRidgeParameter = 0;
  unsigned int newIndex;
  while( numOfIterations < m_maxNumOfIterations && remainingIndices.size() > 0 && ssdImproved )
  {
    ssdImproved = false;
    VnlMatrixType* bestTestSubset = NULL;
    VnlMatrixType* oldBestSubset = bestSubset;

    for( std::set<unsigned int>::iterator it = remainingIndices.begin(); it != remainingIndices.end(); ++it )
    {
      VnlMatrixType* testSubset = NULL;
      std::set<unsigned int> testSubsetIndices( m_optimalIndices );
      testSubsetIndices.insert( *it );
      //combine surrogates
      if( numOfIterations > 0 )
      {
        testSubset = CombineMatrices( oldBestSubset, m_regressorMatrices[*it] );
      }
      else
      {
        testSubset = new VnlMatrixType( *m_regressorMatrices[*it] );
      }

      std::cout<<"CombinedMatrix: "<<*testSubset<<std::endl;

//      m_predictor->SetRegressorTrainingMatrix( *testSubset );
      double bestRidgeSSDValue = std::numeric_limits<double>::max();
      float bestRidgeParameter = 0;


      std::cout << "Test indices: ";
      for( std::set<unsigned int>::iterator it2 = testSubsetIndices.begin(); it2 != testSubsetIndices.end(); ++it2 )
      {
        std::cout << *it2 << "," << std::endl;
      }
      std::cout << std::endl;
      std::cout << "Matrix size: " << testSubset->rows() << " x " << testSubset->cols() << std::endl;

      //initialize predictors for cross validation
      VnlMatrixType newRegressor;
      newRegressor.set_size(testSubset->rows(),testSubset->cols()-1);
      unsigned int predictorIndex=0;
      for (unsigned int col=0; col<testSubset->cols(); col+=useEveryNthCol)
      {
        unsigned int index=0;
        for (unsigned int j=0; j<testSubset->cols();j++)
        {
          if(col!=j)
          {
            newRegressor.set_column(index,testSubset->get_column(j));
            index++;
          }
        }
        m_predictors[predictorIndex]->SetRegressorTrainingMatrix(newRegressor);
        m_predictors[predictorIndex]->TrainLSEstimator();
        predictorIndex++;
      }

      //optimize ridge parameter
      double lastSSDValue=std::numeric_limits<double>::max();
      double currSSDValue = 0;
      for( float ridgeParameter = m_ridgeParameterMin; ridgeParameter <= m_ridgeParameterMax && currSSDValue < lastSSDValue; ridgeParameter*=m_ridgeMultiplicator )
      {

        if(ridgeParameter == m_ridgeParameterMin)
        {
          //m_predictor->TrainLSEstimator();
        } else 
        {
            lastSSDValue=currSSDValue;
        }
        currSSDValue = 0;
        unsigned int numOfSSDTests=0;
        //cross validation
          predictorIndex=0;
          for (unsigned int col=0; col<testSubset->cols(); col+=useEveryNthCol)
          {
              std::cout<<"Cross validation column: "<<col<<std::endl;
              m_predictors[predictorIndex]->SetTikhonovRegularizationParameter( ridgeParameter );
              VnlMatrixType output;
              VnlMatrixType column = testSubset->get_n_columns( col, 1 );
              m_predictors[predictorIndex]->PredictOutput( column, output );
              VnlMatrixType diff = output - m_regressandMatrix->get_n_columns( col, 1 );
              VnlMatrixType ssd = diff.transpose() * diff;
              currSSDValue += ssd[0][0];
              numOfSSDTests++;
            predictorIndex++;
          }
        std::cout<<"SSD: "<< currSSDValue<<std::endl;
        long elems=(long)m_regressandMatrix->rows()*numOfSSDTests;
        std::cout<<"Elements: "<< elems<<std::endl;
        std::cout<<"SSD per element: "<< currSSDValue/elems<<std::endl;
/*        for( unsigned int col = 0; col < testSubset->cols(); col+=2 )
        {
          VnlMatrixType output;
          VnlMatrixType column = testSubset->get_n_columns( col, 1 );
          m_predictor->PredictOutput( column, output );
          //double-valued diff computation
          VnlMatrixType diff = output - m_regressandMatrix->get_n_columns( col, 1 );
          VnlMatrixType ssd = diff.transpose() * diff;
          currSSDValue += ssd[0][0];
          numOfSSDTests++;
        }
        std::cout<<"SSD: "<< currSSDValue<<std::endl;
        std::cout<<"SSD per element: "<< currSSDValue/(m_regressandMatrix->rows()*          numOfSSDTests)<<std::endl;*/
        if( currSSDValue < bestRidgeSSDValue )
        {
          bestRidgeParameter = ridgeParameter;
          bestRidgeSSDValue = currSSDValue;
        }
      }
      /*std::cout << "HRidge parameter: " << bestRidgeParameter << std::endl;
      std::cout << "SSD: " << bestRidgeSSDValue << std::endl;*/

      if( bestRidgeSSDValue < bestSubsetSSDValue )
      {
        bestSubset = testSubset;
        bestSubsetSSDValue = bestRidgeSSDValue;
        bestSubsetRidgeParameter = bestRidgeParameter;
        newIndex = *it;
        ssdImproved = true;
      }
      else
      {
        delete testSubset;
      }
    }




    if( ssdImproved )
    {
      m_optimalIndices.insert( newIndex );
      remainingIndices.erase( newIndex );
      m_optimalSSD = bestSubsetSSDValue;
      m_ridgeParameter = bestSubsetRidgeParameter;
      m_optimalSubset = bestSubset;
      /*std::cout << "HIER2" << std::endl;
      std::cout<<"newIndex: "<<newIndex<<std::endl;
      std::cout<<"size: "<<m_optimalIndices.size()<<std::endl;
      std::cout << "Optimal indices: ";
      for( std::set<unsigned int>::iterator it2 = m_optimalIndices.begin(); it2 != m_optimalIndices.end(); ++it2 )
      {
        std::cout << *it2 << "," << std::endl;
      }
      std::cout<<"Optimal SSD: "<<m_optimalSSD<<std::endl;
      std::cout<<"Optimal ridge: "<<m_ridgeParameter<<std::endl;
      std::cout<<"Optimal subset: "<<*m_optimalSubset<<std::endl;*/
    }
    numOfIterations++;
  }
}

VnlMatrixType* imiMotionPredictionSFS::CombineMatrices( VnlMatrixType* matrixA, VnlMatrixType* matrixB )
{
  if(matrixA->cols()!=matrixB->cols() && matrixA->cols() != 0 && matrixB->cols() != 0)
  {
    std::cerr << "Cannot combine matrices! The matrices do not have the same number of columns!" << std::endl;
    return NULL;
  }

  //std::cout<<"MatrixA "<< matrixA->rows()<<"x"<< matrixA->cols()<<std::endl;
  //std::cout<<"MatrixB "<< matrixB->rows()<<"x"<< matrixB->cols()<<std::endl;
  VnlMatrixType* combinedMatrix = new VnlMatrixType( matrixA->rows() + matrixB->rows(), std::max(matrixA->cols(),matrixB->cols()));

  //copy contents of matrixA
  for( unsigned int row = 0; row < matrixA->rows(); row++ )
  {
    for( unsigned int col = 0; col < matrixA->cols(); col++ )
    {
      combinedMatrix->put( row, col, matrixA->get( row, col ) );
    }
  }

  //copy contents of matrixB
  for( unsigned int row = 0; row < matrixB->rows(); row++ )
  {
    for( unsigned int col = 0; col < matrixB->cols(); col++ )
    {
      combinedMatrix->put( row + matrixA->rows(), col, matrixB->get( row, col ) );
    }
  }

  return combinedMatrix;
}

} //endnamespace
