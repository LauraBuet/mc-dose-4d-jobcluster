/** \file imiRPCA.cpp
 *
 *  \b Initial \b Author: Matthias Wilms, Jonas Ortmueller \n\n
 *  \b Copyright (C) 2015 Institute of Medical Informatics,
 *     University of Luebeck
 *
 ****************************************************************************/

// Project includes:
#include "icnsRobustPCA.h"

// ITK includes

// VNL includes:

namespace imi
{

/** \brief CreateimiRPCAse Delete() to destroy. */
imiRPCA* imiRPCA::New()
{
  return new imiRPCA();
}

imiRPCA::imiRPCA()
{
  m_Lambda = 1;
  m_Tolerance = 1e-7;
  m_MaxIterations = 1000;
}

imiRPCA::~imiRPCA()
{
  // TODO Auto-generated destructor stub
}

bool imiRPCA::Execute()
{
  std::cout << "imiRPCA::Execute()" << std::endl;

  //
  // Make some checks
  //
  //assert( 0 < m_D.n_elem );

  //Init
  //-------------------------------------------------------------
  int rows = m_D.n_rows;
  int cols = m_D.n_cols;

  // Low Rank and Sparse Matrix
  ArmaMatrixType A_hat;
  A_hat.set_size( rows, cols );
  A_hat.zeros();

  ArmaMatrixType E_hat;
  E_hat.set_size( rows, cols );
  E_hat.zeros();

  ArmaMatrixType Z;
  Z.set_size( rows, cols );
  Z.zeros();
  //--------------------------------------------------------------

  // Lagrange Multipliers Matrix Y
  //-----------------------------------------------------------

  ArmaMatrixType Y;
  Y = m_D;

  // Calculated the 2-Norm of Y
  double Y_2_Norm = norm( Y, 2 );
  // Calculated the inf-Norm of Y
  double Y_Inf_Norm = norm( vectorise( Y ), "inf" );

  // computes J(Y) = max(||Y||_2,\gamma^{-1}||Y||_{\inf}), cf. [2], eq. (10)
  double J_of_norm = std::max( Y_2_Norm, (Y_Inf_Norm / m_Lambda) );
  // computes line 1 in Alg. 5 of [2]
  Y = Y / J_of_norm;

  // Set Prorgamm Parameters
  //------------------------------------------------------------------

  double m = 1.25 / Y_2_Norm;
  double m_b = m * 1e07;
  double rho = 1.5;
  int iter = 0;
  bool converged = false;

  // Frobenius Norm of the Intput Matrix
  //--------------------------------------------------------------------
  double D_Fro_Norm = arma::norm( m_D, "fro" );

  ArmaMatrixType Temp_Matrix;
  ArmaMatrixType Temp_Matrix1;
  ArmaMatrixType Temp_Matrix2;
  while( converged == false )
  {

    if( iter % 5 == 0 && iter > 0)
    {
      std::cout << "======================================================" << std::endl;
      std::cout << "  Iteration: " << iter << " ..." << std::endl;
      std::cout << "  Error: " << norm(Z,"fro")/D_Fro_Norm << " ..." << std::endl;
      std::cout << "======================================================" << std::endl;
    }

    Temp_Matrix = m_D - A_hat + (Y * (1 / m));

    Temp_Matrix1 = Temp_Matrix - (m_Lambda / m);

    Temp_Matrix2 = Temp_Matrix + (m_Lambda / m);

    // change elements of Temp_Matrix < 0  to 0
    Temp_Matrix1.elem( find( Temp_Matrix1 < 0 ) ).zeros();

    // change elements of Temp_Matrix1 > 0  to 0
    Temp_Matrix2.elem( find( Temp_Matrix2 > 0 ) ).zeros();

    // calculate Matrix E_hat
    E_hat = Temp_Matrix1;
    E_hat = E_hat + Temp_Matrix2;


    // calculate SVD
    ArmaMatrixType U;
    ArmaMatrixType V;
    ArmaVectorType S;
    ArmaVectorType S_temp;
    ArmaMatrixType S1;

    arma::svd_econ( U, S, V, (m_D - E_hat + (Y * (1 / m))) );

    // number of values greater then 1/m
    int num = 0;
    for( int i = 0; i <= (S.n_rows) - 1; i++ )
    {
      if( S( i, 0 ) > (1 / m) )
      {
        num++;
      }
    }


    // Calculate A_hat

    S_temp = S.rows( 0, std::max(0,num - 1) );
    S1.set_size( std::max(0,num - 1), std::max(0,num - 1) );
    S1 = diagmat( (S_temp - (1 / m)) );

    // Transpose because of V*S*V^t
    V = V.t();

    A_hat = U.cols( 0, std::max(0,num - 1) ) * S1 * V.rows( 0, std::max(0,num - 1) );

    // line 8  Algorithms calculate Y
    Z = m_D - A_hat - E_hat;
    Y = Y + (m * Z);

    // new m
    m = std::min( m * rho, m_b );

    // converged?
    iter++;
    if( norm( Z, "fro" ) / D_Fro_Norm < m_Tolerance )
    {
      converged = true;
    }
  }

  m_L=A_hat;
  m_S=E_hat;

  return true;
}

}

