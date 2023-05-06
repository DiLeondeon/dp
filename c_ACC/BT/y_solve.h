#include "header.h"
#include "work_lhs.h"
#include "timers.h"

//---------------------------------------------------------------------
// Performs line solves in Y direction by first factoring
// the block-tridiagonal matrix into an upper triangular matrix, 
// and then performing back substitution to solve for the unknow
// vectors of each line.  
// 
// Make sure we treat elements zero to cell_size in the direction
// of the sweep.
//---------------------------------------------------------------------
void y_solve()
{
  int i, j, k, m, n, isize, jsize, ksize;

  double tmp1, tmp2, tmp3;
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  if (timeron) timer_start(t_ysolve);

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // This function computes the left hand side for the three y-factors   
  //---------------------------------------------------------------------

  jsize = grid_points[1]-1;
  ksize = grid_points[2]-2;
  isize = grid_points[0]-2;
  //---------------------------------------------------------------------
  // Compute the indices for storing the tri-diagonal matrix;
  // determine a (labeled f) and n jacobians for cell c
  //---------------------------------------------------------------------
  #pragma acc parallel loop collapse(2) private(i,j,k,m,n,tmp1,tmp2,tmp3)
  for (k = 1; k <= ksize; k++) {
    //#pragma acc loop
    for (i = 1; i <= isize; i++) {
      //#pragma acc loop
      for (j = 0; j <= jsize; j++) {
        tmp1 = rho_i[k][j][i];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;

        fjac[k][i][j][0][0] = 0.0;
        fjac[k][i][j][1][0] = 0.0;
        fjac[k][i][j][2][0] = 1.0;
        fjac[k][i][j][3][0] = 0.0;
        fjac[k][i][j][4][0] = 0.0;

        fjac[k][i][j][0][1] = - ( u[k][j][i][1]*u[k][j][i][2] ) * tmp2;
        fjac[k][i][j][1][1] = u[k][j][i][2] * tmp1;
        fjac[k][i][j][2][1] = u[k][j][i][1] * tmp1;
        fjac[k][i][j][3][1] = 0.0;
        fjac[k][i][j][4][1] = 0.0;

        fjac[k][i][j][0][2] = - ( u[k][j][i][2]*u[k][j][i][2]*tmp2)
          + c2 * qs[k][j][i];
        fjac[k][i][j][1][2] = - c2 *  u[k][j][i][1] * tmp1;
        fjac[k][i][j][2][2] = ( 2.0 - c2 ) *  u[k][j][i][2] * tmp1;
        fjac[k][i][j][3][2] = - c2 * u[k][j][i][3] * tmp1;
        fjac[k][i][j][4][2] = c2;

        fjac[k][i][j][0][3] = - ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
        fjac[k][i][j][1][3] = 0.0;
        fjac[k][i][j][2][3] = u[k][j][i][3] * tmp1;
        fjac[k][i][j][3][3] = u[k][j][i][2] * tmp1;
        fjac[k][i][j][4][3] = 0.0;

        fjac[k][i][j][0][4] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
          * u[k][j][i][2] * tmp2;
        fjac[k][i][j][1][4] = - c2 * u[k][j][i][1]*u[k][j][i][2] * tmp2;
        fjac[k][i][j][2][4] = c1 * u[k][j][i][4] * tmp1 
          - c2 * ( qs[k][j][i] + u[k][j][i][2]*u[k][j][i][2] * tmp2 );
        fjac[k][i][j][3][4] = - c2 * ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
        fjac[k][i][j][4][4] = c1 * u[k][j][i][2] * tmp1;

        njac[k][i][j][0][0] = 0.0;
        njac[k][i][j][1][0] = 0.0;
        njac[k][i][j][2][0] = 0.0;
        njac[k][i][j][3][0] = 0.0;
        njac[k][i][j][4][0] = 0.0;

        njac[k][i][j][0][1] = - c3c4 * tmp2 * u[k][j][i][1];
        njac[k][i][j][1][1] =   c3c4 * tmp1;
        njac[k][i][j][2][1] =   0.0;
        njac[k][i][j][3][1] =   0.0;
        njac[k][i][j][4][1] =   0.0;

        njac[k][i][j][0][2] = - con43 * c3c4 * tmp2 * u[k][j][i][2];
        njac[k][i][j][1][2] =   0.0;
        njac[k][i][j][2][2] =   con43 * c3c4 * tmp1;
        njac[k][i][j][3][2] =   0.0;
        njac[k][i][j][4][2] =   0.0;

        njac[k][i][j][0][3] = - c3c4 * tmp2 * u[k][j][i][3];
        njac[k][i][j][1][3] =   0.0;
        njac[k][i][j][2][3] =   0.0;
        njac[k][i][j][3][3] =   c3c4 * tmp1;
        njac[k][i][j][4][3] =   0.0;

        njac[k][i][j][0][4] = - (  c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( con43 * c3c4
              - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[k][i][j][1][4] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][1];
        njac[k][i][j][2][4] = ( con43 * c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[k][i][j][3][4] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][3];
        njac[k][i][j][4][4] = ( c1345 ) * tmp1;
      }

      //---------------------------------------------------------------------
      // now joacobians set, so form left hand side in y direction
      //---------------------------------------------------------------------
      //#pragma acc routine (lhsinit) worker
      //lhsinit(k, i, jsize, lhs);//y_lhsinit(lhs[k][i], jsize);
      for (n = 0; n < 5; n++) {
        for (m = 0; m < 5; m++) {
          lhs[k][i][0][0][n][m] = 0.0;
          lhs[k][i][0][1][n][m] = 0.0;
          lhs[k][i][0][2][n][m] = 0.0;
          lhs[k][i][jsize][0][n][m] = 0.0;
          lhs[k][i][jsize][1][n][m] = 0.0;
          lhs[k][i][jsize][2][n][m] = 0.0;
        }
        lhs[k][i][0][1][n][n] = 1.0;
        lhs[k][i][jsize][1][n][n] = 1.0;
      }
      for (j = 1; j <= jsize-1; j++) {
        tmp1 = dt * ty1;
        tmp2 = dt * ty2;

        lhs[k][i][j][AA][0][0] = - tmp2 * fjac[k][i][j-1][0][0]
          - tmp1 * njac[k][i][j-1][0][0]
          - tmp1 * dy1; 
        lhs[k][i][j][AA][1][0] = - tmp2 * fjac[k][i][j-1][1][0]
          - tmp1 * njac[k][i][j-1][1][0];
        lhs[k][i][j][AA][2][0] = - tmp2 * fjac[k][i][j-1][2][0]
          - tmp1 * njac[k][i][j-1][2][0];
        lhs[k][i][j][AA][3][0] = - tmp2 * fjac[k][i][j-1][3][0]
          - tmp1 * njac[k][i][j-1][3][0];
        lhs[k][i][j][AA][4][0] = - tmp2 * fjac[k][i][j-1][4][0]
          - tmp1 * njac[k][i][j-1][4][0];

        lhs[k][i][j][AA][0][1] = - tmp2 * fjac[k][i][j-1][0][1]
          - tmp1 * njac[k][i][j-1][0][1];
        lhs[k][i][j][AA][1][1] = - tmp2 * fjac[k][i][j-1][1][1]
          - tmp1 * njac[k][i][j-1][1][1]
          - tmp1 * dy2;
        lhs[k][i][j][AA][2][1] = - tmp2 * fjac[k][i][j-1][2][1]
          - tmp1 * njac[k][i][j-1][2][1];
        lhs[k][i][j][AA][3][1] = - tmp2 * fjac[k][i][j-1][3][1]
          - tmp1 * njac[k][i][j-1][3][1];
        lhs[k][i][j][AA][4][1] = - tmp2 * fjac[k][i][j-1][4][1]
          - tmp1 * njac[k][i][j-1][4][1];

        lhs[k][i][j][AA][0][2] = - tmp2 * fjac[k][i][j-1][0][2]
          - tmp1 * njac[k][i][j-1][0][2];
        lhs[k][i][j][AA][1][2] = - tmp2 * fjac[k][i][j-1][1][2]
          - tmp1 * njac[k][i][j-1][1][2];
        lhs[k][i][j][AA][2][2] = - tmp2 * fjac[k][i][j-1][2][2]
          - tmp1 * njac[k][i][j-1][2][2]
          - tmp1 * dy3;
        lhs[k][i][j][AA][3][2] = - tmp2 * fjac[k][i][j-1][3][2]
          - tmp1 * njac[k][i][j-1][3][2];
        lhs[k][i][j][AA][4][2] = - tmp2 * fjac[k][i][j-1][4][2]
          - tmp1 * njac[k][i][j-1][4][2];

        lhs[k][i][j][AA][0][3] = - tmp2 * fjac[k][i][j-1][0][3]
          - tmp1 * njac[k][i][j-1][0][3];
        lhs[k][i][j][AA][1][3] = - tmp2 * fjac[k][i][j-1][1][3]
          - tmp1 * njac[k][i][j-1][1][3];
        lhs[k][i][j][AA][2][3] = - tmp2 * fjac[k][i][j-1][2][3]
          - tmp1 * njac[k][i][j-1][2][3];
        lhs[k][i][j][AA][3][3] = - tmp2 * fjac[k][i][j-1][3][3]
          - tmp1 * njac[k][i][j-1][3][3]
          - tmp1 * dy4;
        lhs[k][i][j][AA][4][3] = - tmp2 * fjac[k][i][j-1][4][3]
          - tmp1 * njac[k][i][j-1][4][3];

        lhs[k][i][j][AA][0][4] = - tmp2 * fjac[k][i][j-1][0][4]
          - tmp1 * njac[k][i][j-1][0][4];
        lhs[k][i][j][AA][1][4] = - tmp2 * fjac[k][i][j-1][1][4]
          - tmp1 * njac[k][i][j-1][1][4];
        lhs[k][i][j][AA][2][4] = - tmp2 * fjac[k][i][j-1][2][4]
          - tmp1 * njac[k][i][j-1][2][4];
        lhs[k][i][j][AA][3][4] = - tmp2 * fjac[k][i][j-1][3][4]
          - tmp1 * njac[k][i][j-1][3][4];
        lhs[k][i][j][AA][4][4] = - tmp2 * fjac[k][i][j-1][4][4]
          - tmp1 * njac[k][i][j-1][4][4]
          - tmp1 * dy5;

        lhs[k][i][j][BB][0][0] = 1.0
          + tmp1 * 2.0 * njac[k][i][j][0][0]
          + tmp1 * 2.0 * dy1;
        lhs[k][i][j][BB][1][0] = tmp1 * 2.0 * njac[k][i][j][1][0];
        lhs[k][i][j][BB][2][0] = tmp1 * 2.0 * njac[k][i][j][2][0];
        lhs[k][i][j][BB][3][0] = tmp1 * 2.0 * njac[k][i][j][3][0];
        lhs[k][i][j][BB][4][0] = tmp1 * 2.0 * njac[k][i][j][4][0];

        lhs[k][i][j][BB][0][1] = tmp1 * 2.0 * njac[k][i][j][0][1];
        lhs[k][i][j][BB][1][1] = 1.0
          + tmp1 * 2.0 * njac[k][i][j][1][1]
          + tmp1 * 2.0 * dy2;
        lhs[k][i][j][BB][2][1] = tmp1 * 2.0 * njac[k][i][j][2][1];
        lhs[k][i][j][BB][3][1] = tmp1 * 2.0 * njac[k][i][j][3][1];
        lhs[k][i][j][BB][4][1] = tmp1 * 2.0 * njac[k][i][j][4][1];

        lhs[k][i][j][BB][0][2] = tmp1 * 2.0 * njac[k][i][j][0][2];
        lhs[k][i][j][BB][1][2] = tmp1 * 2.0 * njac[k][i][j][1][2];
        lhs[k][i][j][BB][2][2] = 1.0
          + tmp1 * 2.0 * njac[k][i][j][2][2]
          + tmp1 * 2.0 * dy3;
        lhs[k][i][j][BB][3][2] = tmp1 * 2.0 * njac[k][i][j][3][2];
        lhs[k][i][j][BB][4][2] = tmp1 * 2.0 * njac[k][i][j][4][2];

        lhs[k][i][j][BB][0][3] = tmp1 * 2.0 * njac[k][i][j][0][3];
        lhs[k][i][j][BB][1][3] = tmp1 * 2.0 * njac[k][i][j][1][3];
        lhs[k][i][j][BB][2][3] = tmp1 * 2.0 * njac[k][i][j][2][3];
        lhs[k][i][j][BB][3][3] = 1.0
          + tmp1 * 2.0 * njac[k][i][j][3][3]
          + tmp1 * 2.0 * dy4;
        lhs[k][i][j][BB][4][3] = tmp1 * 2.0 * njac[k][i][j][4][3];

        lhs[k][i][j][BB][0][4] = tmp1 * 2.0 * njac[k][i][j][0][4];
        lhs[k][i][j][BB][1][4] = tmp1 * 2.0 * njac[k][i][j][1][4];
        lhs[k][i][j][BB][2][4] = tmp1 * 2.0 * njac[k][i][j][2][4];
        lhs[k][i][j][BB][3][4] = tmp1 * 2.0 * njac[k][i][j][3][4];
        lhs[k][i][j][BB][4][4] = 1.0
          + tmp1 * 2.0 * njac[k][i][j][4][4] 
          + tmp1 * 2.0 * dy5;

        lhs[k][i][j][CC][0][0] =  tmp2 * fjac[k][i][j+1][0][0]
          - tmp1 * njac[k][i][j+1][0][0]
          - tmp1 * dy1;
        lhs[k][i][j][CC][1][0] =  tmp2 * fjac[k][i][j+1][1][0]
          - tmp1 * njac[k][i][j+1][1][0];
        lhs[k][i][j][CC][2][0] =  tmp2 * fjac[k][i][j+1][2][0]
          - tmp1 * njac[k][i][j+1][2][0];
        lhs[k][i][j][CC][3][0] =  tmp2 * fjac[k][i][j+1][3][0]
          - tmp1 * njac[k][i][j+1][3][0];
        lhs[k][i][j][CC][4][0] =  tmp2 * fjac[k][i][j+1][4][0]
          - tmp1 * njac[k][i][j+1][4][0];

        lhs[k][i][j][CC][0][1] =  tmp2 * fjac[k][i][j+1][0][1]
          - tmp1 * njac[k][i][j+1][0][1];
        lhs[k][i][j][CC][1][1] =  tmp2 * fjac[k][i][j+1][1][1]
          - tmp1 * njac[k][i][j+1][1][1]
          - tmp1 * dy2;
        lhs[k][i][j][CC][2][1] =  tmp2 * fjac[k][i][j+1][2][1]
          - tmp1 * njac[k][i][j+1][2][1];
        lhs[k][i][j][CC][3][1] =  tmp2 * fjac[k][i][j+1][3][1]
          - tmp1 * njac[k][i][j+1][3][1];
        lhs[k][i][j][CC][4][1] =  tmp2 * fjac[k][i][j+1][4][1]
          - tmp1 * njac[k][i][j+1][4][1];

        lhs[k][i][j][CC][0][2] =  tmp2 * fjac[k][i][j+1][0][2]
          - tmp1 * njac[k][i][j+1][0][2];
        lhs[k][i][j][CC][1][2] =  tmp2 * fjac[k][i][j+1][1][2]
          - tmp1 * njac[k][i][j+1][1][2];
        lhs[k][i][j][CC][2][2] =  tmp2 * fjac[k][i][j+1][2][2]
          - tmp1 * njac[k][i][j+1][2][2]
          - tmp1 * dy3;
        lhs[k][i][j][CC][3][2] =  tmp2 * fjac[k][i][j+1][3][2]
          - tmp1 * njac[k][i][j+1][3][2];
        lhs[k][i][j][CC][4][2] =  tmp2 * fjac[k][i][j+1][4][2]
          - tmp1 * njac[k][i][j+1][4][2];

        lhs[k][i][j][CC][0][3] =  tmp2 * fjac[k][i][j+1][0][3]
          - tmp1 * njac[k][i][j+1][0][3];
        lhs[k][i][j][CC][1][3] =  tmp2 * fjac[k][i][j+1][1][3]
          - tmp1 * njac[k][i][j+1][1][3];
        lhs[k][i][j][CC][2][3] =  tmp2 * fjac[k][i][j+1][2][3]
          - tmp1 * njac[k][i][j+1][2][3];
        lhs[k][i][j][CC][3][3] =  tmp2 * fjac[k][i][j+1][3][3]
          - tmp1 * njac[k][i][j+1][3][3]
          - tmp1 * dy4;
        lhs[k][i][j][CC][4][3] =  tmp2 * fjac[k][i][j+1][4][3]
          - tmp1 * njac[k][i][j+1][4][3];

        lhs[k][i][j][CC][0][4] =  tmp2 * fjac[k][i][j+1][0][4]
          - tmp1 * njac[k][i][j+1][0][4];
        lhs[k][i][j][CC][1][4] =  tmp2 * fjac[k][i][j+1][1][4]
          - tmp1 * njac[k][i][j+1][1][4];
        lhs[k][i][j][CC][2][4] =  tmp2 * fjac[k][i][j+1][2][4]
          - tmp1 * njac[k][i][j+1][2][4];
        lhs[k][i][j][CC][3][4] =  tmp2 * fjac[k][i][j+1][3][4]
          - tmp1 * njac[k][i][j+1][3][4];
        lhs[k][i][j][CC][4][4] =  tmp2 * fjac[k][i][j+1][4][4]
          - tmp1 * njac[k][i][j+1][4][4]
          - tmp1 * dy5;
      }

      //---------------------------------------------------------------------
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // performs guaussian elimination on this cell.
      // 
      // assumes that unpacking routines for non-first cells 
      // preload C' and rhs' from previous cell.
      // 
      // assumed send happens outside this routine, but that
      // c'(JMAX) and rhs'(JMAX) will be sent to next cell
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // multiply c[k][0][i] by b_inverse and copy back to c
      // multiply rhs(0) by b_inverse(0) and copy to rhs
      //---------------------------------------------------------------------
      //#pragma acc routine (binvcrhs) worker
      binvcrhs( k, i, 0, BB, k, i, 0, CC, k, 0, i, lhs, rhs );//y_binvcrhs( lhs[k][i][0][BB], lhs[k][i][0][CC], rhs[k][0][i] );

      //---------------------------------------------------------------------
      // begin inner most do loop
      // do all the elements of the cell unless last 
      //---------------------------------------------------------------------
      for (j = 1; j <= jsize-1; j++) {
        //-------------------------------------------------------------------
        // subtract A*lhs_vector(j-1) from lhs_vector(j)
        // 
        // rhs(j) = rhs(j) - A*rhs(j-1)
        //-------------------------------------------------------------------
        //#pragma acc routine (matvec_sub) worker
        //matvec_sub(k, i, j, AA, k, j-1, i, k, j, i, lhs, rhs);//y_matvec_sub(lhs[k][i][j][AA], rhs[k][j-1][i], rhs[k][j][i]);
        rhs[k][j][i][0] = rhs[k][j][i][0] - lhs[k][i][j][AA][0][0]*rhs[k][j-1][i][0]
                    - lhs[k][i][j][AA][1][0]*rhs[k][j-1][i][1]
                    - lhs[k][i][j][AA][2][0]*rhs[k][j-1][i][2]
                    - lhs[k][i][j][AA][3][0]*rhs[k][j-1][i][3]
                    - lhs[k][i][j][AA][4][0]*rhs[k][j-1][i][4];
        rhs[k][j][i][1] = rhs[k][j][i][1] - lhs[k][i][j][AA][0][1]*rhs[k][j-1][i][0]
                    - lhs[k][i][j][AA][1][1]*rhs[k][j-1][i][1]
                    - lhs[k][i][j][AA][2][1]*rhs[k][j-1][i][2]
                    - lhs[k][i][j][AA][3][1]*rhs[k][j-1][i][3]
                    - lhs[k][i][j][AA][4][1]*rhs[k][j-1][i][4];
        rhs[k][j][i][2] = rhs[k][j][i][2] - lhs[k][i][j][AA][0][2]*rhs[k][j-1][i][0]
                    - lhs[k][i][j][AA][1][2]*rhs[k][j-1][i][1]
                    - lhs[k][i][j][AA][2][2]*rhs[k][j-1][i][2]
                    - lhs[k][i][j][AA][3][2]*rhs[k][j-1][i][3]
                    - lhs[k][i][j][AA][4][2]*rhs[k][j-1][i][4];
        rhs[k][j][i][3] = rhs[k][j][i][3] - lhs[k][i][j][AA][0][3]*rhs[k][j-1][i][0]
                    - lhs[k][i][j][AA][1][3]*rhs[k][j-1][i][1]
                    - lhs[k][i][j][AA][2][3]*rhs[k][j-1][i][2]
                    - lhs[k][i][j][AA][3][3]*rhs[k][j-1][i][3]
                    - lhs[k][i][j][AA][4][3]*rhs[k][j-1][i][4];
        rhs[k][j][i][4] = rhs[k][j][i][4] - lhs[k][i][j][AA][0][4]*rhs[k][j-1][i][0]
                    - lhs[k][i][j][AA][1][4]*rhs[k][j-1][i][1]
                    - lhs[k][i][j][AA][2][4]*rhs[k][j-1][i][2]
                    - lhs[k][i][j][AA][3][4]*rhs[k][j-1][i][3]
                    - lhs[k][i][j][AA][4][4]*rhs[k][j-1][i][4];

        //-------------------------------------------------------------------
        // B(j) = B(j) - C(j-1)*A(j)
        //-------------------------------------------------------------------
        //#pragma acc routine (matmul_sub) worker
        matmul_sub(k, i, j, AA, k, i, j-1, CC, k, i, j, BB, lhs);//y_matmul_sub(lhs[k][i][j][AA], lhs[k][i][j-1][CC], lhs[k][i][j][BB]);

        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[k][0][i] by b_inverse[k][0][i] and copy to rhs
        //-------------------------------------------------------------------
        //#pragma acc routine (binvcrhs) worker
        binvcrhs( k, i, j, BB, k, i, j, CC, k, j, i, lhs, rhs );//y_binvcrhs( lhs[k][i][j][BB], lhs[k][i][j][CC], rhs[k][j][i] );
      }

      //---------------------------------------------------------------------
      // rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
      //---------------------------------------------------------------------
      //#pragma acc routine (matvec_sub) worker
      //matvec_sub(k, i, jsize, AA, k, jsize-1, i, k, jsize, i, lhs, rhs);//y_matvec_sub(lhs[k][i][jsize][AA], rhs[k][jsize-1][i], rhs[k][jsize][i]);
      rhs[k][jsize][i][0] = rhs[k][jsize][i][0] - lhs[k][i][jsize][AA][0][0]*rhs[k][jsize-1][i][0]
                    - lhs[k][i][jsize][AA][1][0]*rhs[k][jsize-1][i][1]
                    - lhs[k][i][jsize][AA][2][0]*rhs[k][jsize-1][i][2]
                    - lhs[k][i][jsize][AA][3][0]*rhs[k][jsize-1][i][3]
                    - lhs[k][i][jsize][AA][4][0]*rhs[k][jsize-1][i][4];
      rhs[k][jsize][i][1] = rhs[k][jsize][i][1] - lhs[k][i][jsize][AA][0][1]*rhs[k][jsize-1][i][0]
                    - lhs[k][i][jsize][AA][1][1]*rhs[k][jsize-1][i][1]
                    - lhs[k][i][jsize][AA][2][1]*rhs[k][jsize-1][i][2]
                    - lhs[k][i][jsize][AA][3][1]*rhs[k][jsize-1][i][3]
                    - lhs[k][i][jsize][AA][4][1]*rhs[k][jsize-1][i][4];
      rhs[k][jsize][i][2] = rhs[k][jsize][i][2] - lhs[k][i][jsize][AA][0][2]*rhs[k][jsize-1][i][0]
                    - lhs[k][i][jsize][AA][1][2]*rhs[k][jsize-1][i][1]
                    - lhs[k][i][jsize][AA][2][2]*rhs[k][jsize-1][i][2]
                    - lhs[k][i][jsize][AA][3][2]*rhs[k][jsize-1][i][3]
                    - lhs[k][i][jsize][AA][4][2]*rhs[k][jsize-1][i][4];
      rhs[k][jsize][i][3] = rhs[k][jsize][i][3] - lhs[k][i][jsize][AA][0][3]*rhs[k][jsize-1][i][0]
                    - lhs[k][i][jsize][AA][1][3]*rhs[k][jsize-1][i][1]
                    - lhs[k][i][jsize][AA][2][3]*rhs[k][jsize-1][i][2]
                    - lhs[k][i][jsize][AA][3][3]*rhs[k][jsize-1][i][3]
                    - lhs[k][i][jsize][AA][4][3]*rhs[k][jsize-1][i][4];
      rhs[k][jsize][i][4] = rhs[k][jsize][i][4] - lhs[k][i][jsize][AA][0][4]*rhs[k][jsize-1][i][0]
                    - lhs[k][i][jsize][AA][1][4]*rhs[k][jsize-1][i][1]
                    - lhs[k][i][jsize][AA][2][4]*rhs[k][jsize-1][i][2]
                    - lhs[k][i][jsize][AA][3][4]*rhs[k][jsize-1][i][3]
                    - lhs[k][i][jsize][AA][4][4]*rhs[k][jsize-1][i][4];

      //---------------------------------------------------------------------
      // B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
      // matmul_sub(AA,i,jsize,k,c,
      // $              CC,i,jsize-1,k,c,BB,i,jsize,k)
      //---------------------------------------------------------------------
      //#pragma acc routine (matmul_sub) worker
      matmul_sub(k, i, jsize, AA, k, i, jsize-1, CC, k, i, jsize, BB, lhs);//y_matmul_sub(lhs[k][i][jsize][AA], lhs[k][i][jsize-1][CC], lhs[k][i][jsize][BB]);

      //---------------------------------------------------------------------
      // multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
      //---------------------------------------------------------------------
      //#pragma acc routine (binvrhs) worker
      binvrhs( k, i, jsize, BB, k, jsize, i, lhs, rhs );//y_binvrhs( lhs[k][i][jsize][BB], rhs[k][jsize][i] );

      //---------------------------------------------------------------------
      // back solve: if last cell, then generate U(jsize)=rhs(jsize)
      // else assume U(jsize) is loaded in un pack backsub_info
      // so just use it
      // after u(jstart) will be sent to next cell
      //---------------------------------------------------------------------
      for (j = jsize-1; j >= 0; j--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] 
              - lhs[k][i][j][CC][n][m]*rhs[k][j+1][i][n];
          }
        }
      }
    }
  }
  if (timeron) timer_stop(t_ysolve);
}
