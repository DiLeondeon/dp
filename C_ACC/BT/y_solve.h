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
  int i, j, k, m, n, jsize;

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

  //---------------------------------------------------------------------
  // Compute the indices for storing the tri-diagonal matrix;
  // determine a (labeled f) and n jacobians for cell c
  //---------------------------------------------------------------------
  #pragma acc parallel loop collapse(2) private(i,j,k,m,n,tmp1,tmp2,tmp3)
  for (k = 1; k <= grid_points[2]-2; k++) {
    //#pragma acc loop
    for (i = 1; i <= grid_points[0]-2; i++) {
      //#pragma acc loop
      for (j = 0; j <= jsize; j++) {
        tmp1 = rho_i[k][j][i];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;

        fjac[0][0][k][i][j] = 0.0; //fjac[0][0][k][i][j]
        fjac[1][0][k][i][j] = 0.0; //fjac[1][0][k][i][j]
        fjac[2][0][k][i][j] = 1.0; //fjac[2][0][k][i][j]
        fjac[3][0][k][i][j] = 0.0; //fjac[3][0][k][i][j]
        fjac[4][0][k][i][j] = 0.0; //fjac[4][0][k][i][j]

        fjac[0][1][k][i][j] = - ( u[k][j][i][1]*u[k][j][i][2] ) * tmp2; //fjac[0][1][k][i][j]
        fjac[1][1][k][i][j] = u[k][j][i][2] * tmp1; //fjac[1][1][k][i][j]
        fjac[2][1][k][i][j] = u[k][j][i][1] * tmp1; //fjac[2][1][k][i][j]
        fjac[3][1][k][i][j] = 0.0; //fjac[3][1][k][i][j]
        fjac[4][1][k][i][j] = 0.0; //fjac[4][1][k][i][j]

        fjac[0][2][k][i][j] = - ( u[k][j][i][2]*u[k][j][i][2]*tmp2) //fjac[0][2][k][i][j]
          + c2 * qs[k][j][i];
        fjac[1][2][k][i][j] = - c2 *  u[k][j][i][1] * tmp1; //fjac[1][2][k][i][j]
        fjac[2][2][k][i][j] = ( 2.0 - c2 ) *  u[k][j][i][2] * tmp1; //fjac[2][2][k][i][j]
        fjac[3][2][k][i][j] = - c2 * u[k][j][i][3] * tmp1; //fjac[3][2][k][i][j]
        fjac[4][2][k][i][j] = c2; //fjac[4][2][k][i][j]

        fjac[0][3][k][i][j] = - ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2; //fjac[0][3][k][i][j]
        fjac[1][3][k][i][j] = 0.0; //fjac[1][3][k][i][j]
        fjac[2][3][k][i][j] = u[k][j][i][3] * tmp1; //fjac[2][3][k][i][j]
        fjac[3][3][k][i][j] = u[k][j][i][2] * tmp1; //fjac[3][3][k][i][j]
        fjac[4][3][k][i][j] = 0.0; //fjac[4][3][k][i][j]

        fjac[0][4][k][i][j] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] ) //fjac[0][4][k][i][j]
          * u[k][j][i][2] * tmp2;
        fjac[1][4][k][i][j] = - c2 * u[k][j][i][1]*u[k][j][i][2] * tmp2; //fjac[1][4][k][i][j]
        fjac[2][4][k][i][j] = c1 * u[k][j][i][4] * tmp1  //fjac[2][4][k][i][j]
          - c2 * ( qs[k][j][i] + u[k][j][i][2]*u[k][j][i][2] * tmp2 );
        fjac[3][4][k][i][j] = - c2 * ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2; //fjac[3][4][k][i][j]
        fjac[4][4][k][i][j] = c1 * u[k][j][i][2] * tmp1; //fjac[4][4][k][i][j]

        njac[0][0][k][i][j] = 0.0;
        njac[1][0][k][i][j] = 0.0;
        njac[2][0][k][i][j] = 0.0;
        njac[3][0][k][i][j] = 0.0;
        njac[4][0][k][i][j] = 0.0;

        njac[0][1][k][i][j] = - c3c4 * tmp2 * u[k][j][i][1];
        njac[1][1][k][i][j] =   c3c4 * tmp1;
        njac[2][1][k][i][j] =   0.0;
        njac[3][1][k][i][j] =   0.0;
        njac[4][1][k][i][j] =   0.0;

        njac[0][2][k][i][j] = - con43 * c3c4 * tmp2 * u[k][j][i][2];
        njac[1][2][k][i][j] =   0.0;
        njac[2][2][k][i][j] =   con43 * c3c4 * tmp1;
        njac[3][2][k][i][j] =   0.0;
        njac[4][2][k][i][j] =   0.0;

        njac[0][3][k][i][j] = - c3c4 * tmp2 * u[k][j][i][3];
        njac[1][3][k][i][j] =   0.0;
        njac[2][3][k][i][j] =   0.0;
        njac[3][3][k][i][j] =   c3c4 * tmp1;
        njac[4][3][k][i][j] =   0.0;

        njac[0][4][k][i][j] = - (  c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( con43 * c3c4
              - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[1][4][k][i][j] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][1];
        njac[2][4][k][i][j] = ( con43 * c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[3][4][k][i][j] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][3];
        njac[4][4][k][i][j] = ( c1345 ) * tmp1;
      }

      //---------------------------------------------------------------------
      // now joacobians set, so form left hand side in y direction
      //---------------------------------------------------------------------
      //#pragma acc routine (lhsinit) worker
      lhsinit(k, i, jsize, lhs);//y_lhsinit(lhs[k][i], jsize);
      for (j = 1; j <= jsize-1; j++) {
        tmp1 = dt * ty1;
        tmp2 = dt * ty2;

        lhs[AA][0][0][k][i][j] = - tmp2 * fjac[0][0][k][i][j-1] //lhs[AA][0][0][k][i][j] fjac[0][0][k][i][j-1]
          - tmp1 * njac[0][0][k][i][j-1]
          - tmp1 * dy1; 
        lhs[AA][1][0][k][i][j] = - tmp2 * fjac[1][0][k][i][j-1]
          - tmp1 * njac[1][0][k][i][j-1];
        lhs[AA][2][0][k][i][j] = - tmp2 * fjac[2][0][k][i][j-1]
          - tmp1 * njac[2][0][k][i][j-1];
        lhs[AA][3][0][k][i][j] = - tmp2 * fjac[3][0][k][i][j-1]
          - tmp1 * njac[3][0][k][i][j-1];
        lhs[AA][4][0][k][i][j] = - tmp2 * fjac[4][0][k][i][j-1]
          - tmp1 * njac[4][0][k][i][j-1];

        lhs[AA][0][1][k][i][j] = - tmp2 * fjac[0][1][k][i][j-1]
          - tmp1 * njac[0][1][k][i][j-1];
        lhs[AA][1][1][k][i][j] = - tmp2 * fjac[1][1][k][i][j-1]
          - tmp1 * njac[1][1][k][i][j-1]
          - tmp1 * dy2;
        lhs[AA][2][1][k][i][j] = - tmp2 * fjac[2][1][k][i][j-1]
          - tmp1 * njac[2][1][k][i][j-1];
        lhs[AA][3][1][k][i][j] = - tmp2 * fjac[3][1][k][i][j-1]
          - tmp1 * njac[3][1][k][i][j-1];
        lhs[AA][4][1][k][i][j] = - tmp2 * fjac[4][1][k][i][j-1]
          - tmp1 * njac[4][1][k][i][j-1];

        lhs[AA][0][2][k][i][j] = - tmp2 * fjac[0][2][k][i][j-1]
          - tmp1 * njac[0][2][k][i][j-1];
        lhs[AA][1][2][k][i][j] = - tmp2 * fjac[1][2][k][i][j-1]
          - tmp1 * njac[1][2][k][i][j-1];
        lhs[AA][2][2][k][i][j] = - tmp2 * fjac[2][2][k][i][j-1]
          - tmp1 * njac[2][2][k][i][j-1]
          - tmp1 * dy3;
        lhs[AA][3][2][k][i][j] = - tmp2 * fjac[3][2][k][i][j-1]
          - tmp1 * njac[3][2][k][i][j-1];
        lhs[AA][4][2][k][i][j] = - tmp2 * fjac[4][2][k][i][j-1]
          - tmp1 * njac[4][2][k][i][j-1];

        lhs[AA][0][3][k][i][j] = - tmp2 * fjac[0][3][k][i][j-1]
          - tmp1 * njac[0][3][k][i][j-1];
        lhs[AA][1][3][k][i][j] = - tmp2 * fjac[1][3][k][i][j-1]
          - tmp1 * njac[1][3][k][i][j-1];
        lhs[AA][2][3][k][i][j] = - tmp2 * fjac[2][3][k][i][j-1]
          - tmp1 * njac[2][3][k][i][j-1];
        lhs[AA][3][3][k][i][j] = - tmp2 * fjac[3][3][k][i][j-1]
          - tmp1 * njac[3][3][k][i][j-1]
          - tmp1 * dy4;
        lhs[AA][4][3][k][i][j] = - tmp2 * fjac[4][3][k][i][j-1]
          - tmp1 * njac[4][3][k][i][j-1];

        lhs[AA][0][4][k][i][j] = - tmp2 * fjac[0][4][k][i][j-1]
          - tmp1 * njac[0][4][k][i][j-1];
        lhs[AA][1][4][k][i][j] = - tmp2 * fjac[1][4][k][i][j-1]
          - tmp1 * njac[1][4][k][i][j-1];
        lhs[AA][2][4][k][i][j] = - tmp2 * fjac[2][4][k][i][j-1]
          - tmp1 * njac[2][4][k][i][j-1];
        lhs[AA][3][4][k][i][j] = - tmp2 * fjac[3][4][k][i][j-1]
          - tmp1 * njac[3][4][k][i][j-1];
        lhs[AA][4][4][k][i][j] = - tmp2 * fjac[4][4][k][i][j-1]
          - tmp1 * njac[4][4][k][i][j-1]
          - tmp1 * dy5;

        lhs[BB][0][0][k][i][j] = 1.0
          + tmp1 * 2.0 * njac[0][0][k][i][j]
          + tmp1 * 2.0 * dy1;
        lhs[BB][1][0][k][i][j] = tmp1 * 2.0 * njac[1][0][k][i][j];
        lhs[BB][2][0][k][i][j] = tmp1 * 2.0 * njac[2][0][k][i][j];
        lhs[BB][3][0][k][i][j] = tmp1 * 2.0 * njac[3][0][k][i][j];
        lhs[BB][4][0][k][i][j] = tmp1 * 2.0 * njac[4][0][k][i][j];

        lhs[BB][0][1][k][i][j] = tmp1 * 2.0 * njac[0][1][k][i][j];
        lhs[BB][1][1][k][i][j] = 1.0
          + tmp1 * 2.0 * njac[1][1][k][i][j]
          + tmp1 * 2.0 * dy2;
        lhs[BB][2][1][k][i][j] = tmp1 * 2.0 * njac[2][1][k][i][j];
        lhs[BB][3][1][k][i][j] = tmp1 * 2.0 * njac[3][1][k][i][j];
        lhs[BB][4][1][k][i][j] = tmp1 * 2.0 * njac[4][1][k][i][j];

        lhs[BB][0][2][k][i][j] = tmp1 * 2.0 * njac[0][2][k][i][j];
        lhs[BB][1][2][k][i][j] = tmp1 * 2.0 * njac[1][2][k][i][j];
        lhs[BB][2][2][k][i][j] = 1.0
          + tmp1 * 2.0 * njac[2][2][k][i][j]
          + tmp1 * 2.0 * dy3;
        lhs[BB][3][2][k][i][j] = tmp1 * 2.0 * njac[3][2][k][i][j];
        lhs[BB][4][2][k][i][j] = tmp1 * 2.0 * njac[4][2][k][i][j];

        lhs[BB][0][3][k][i][j] = tmp1 * 2.0 * njac[0][3][k][i][j];
        lhs[BB][1][3][k][i][j] = tmp1 * 2.0 * njac[1][3][k][i][j];
        lhs[BB][2][3][k][i][j] = tmp1 * 2.0 * njac[2][3][k][i][j];
        lhs[BB][3][3][k][i][j] = 1.0
          + tmp1 * 2.0 * njac[3][3][k][i][j]
          + tmp1 * 2.0 * dy4;
        lhs[BB][4][3][k][i][j] = tmp1 * 2.0 * njac[4][3][k][i][j];

        lhs[BB][0][4][k][i][j] = tmp1 * 2.0 * njac[0][4][k][i][j];
        lhs[BB][1][4][k][i][j] = tmp1 * 2.0 * njac[1][4][k][i][j];
        lhs[BB][2][4][k][i][j] = tmp1 * 2.0 * njac[2][4][k][i][j];
        lhs[BB][3][4][k][i][j] = tmp1 * 2.0 * njac[3][4][k][i][j];
        lhs[BB][4][4][k][i][j] = 1.0
          + tmp1 * 2.0 * njac[4][4][k][i][j] 
          + tmp1 * 2.0 * dy5;

        lhs[CC][0][0][k][i][j] =  tmp2 * fjac[0][0][k][i][j+1]
          - tmp1 * njac[0][0][k][i][j+1]
          - tmp1 * dy1;
        lhs[CC][1][0][k][i][j] =  tmp2 * fjac[1][0][k][i][j+1]
          - tmp1 * njac[1][0][k][i][j+1];
        lhs[CC][2][0][k][i][j] =  tmp2 * fjac[2][0][k][i][j+1]
          - tmp1 * njac[2][0][k][i][j+1];
        lhs[CC][3][0][k][i][j] =  tmp2 * fjac[3][0][k][i][j+1]
          - tmp1 * njac[3][0][k][i][j+1];
        lhs[CC][4][0][k][i][j] =  tmp2 * fjac[4][0][k][i][j+1]
          - tmp1 * njac[4][0][k][i][j+1];

        lhs[CC][0][1][k][i][j] =  tmp2 * fjac[0][1][k][i][j+1]
          - tmp1 * njac[0][1][k][i][j+1];
        lhs[CC][1][1][k][i][j] =  tmp2 * fjac[1][1][k][i][j+1]
          - tmp1 * njac[1][1][k][i][j+1]
          - tmp1 * dy2;
        lhs[CC][2][1][k][i][j] =  tmp2 * fjac[2][1][k][i][j+1]
          - tmp1 * njac[2][1][k][i][j+1];
        lhs[CC][3][1][k][i][j] =  tmp2 * fjac[3][1][k][i][j+1]
          - tmp1 * njac[3][1][k][i][j+1];
        lhs[CC][4][1][k][i][j] =  tmp2 * fjac[4][1][k][i][j+1]
          - tmp1 * njac[4][1][k][i][j+1];

        lhs[CC][0][2][k][i][j] =  tmp2 * fjac[0][2][k][i][j+1]
          - tmp1 * njac[0][2][k][i][j+1];
        lhs[CC][1][2][k][i][j] =  tmp2 * fjac[1][2][k][i][j+1]
          - tmp1 * njac[1][2][k][i][j+1];
        lhs[CC][2][2][k][i][j] =  tmp2 * fjac[2][2][k][i][j+1]
          - tmp1 * njac[2][2][k][i][j+1]
          - tmp1 * dy3;
        lhs[CC][3][2][k][i][j] =  tmp2 * fjac[3][2][k][i][j+1]
          - tmp1 * njac[3][2][k][i][j+1];
        lhs[CC][4][2][k][i][j] =  tmp2 * fjac[4][2][k][i][j+1]
          - tmp1 * njac[4][2][k][i][j+1];

        lhs[CC][0][3][k][i][j] =  tmp2 * fjac[0][3][k][i][j+1]
          - tmp1 * njac[0][3][k][i][j+1];
        lhs[CC][1][3][k][i][j] =  tmp2 * fjac[1][3][k][i][j+1]
          - tmp1 * njac[1][3][k][i][j+1];
        lhs[CC][2][3][k][i][j] =  tmp2 * fjac[2][3][k][i][j+1]
          - tmp1 * njac[2][3][k][i][j+1];
        lhs[CC][3][3][k][i][j] =  tmp2 * fjac[3][3][k][i][j+1]
          - tmp1 * njac[3][3][k][i][j+1]
          - tmp1 * dy4;
        lhs[CC][4][3][k][i][j] =  tmp2 * fjac[4][3][k][i][j+1]
          - tmp1 * njac[4][3][k][i][j+1];

        lhs[CC][0][4][k][i][j] =  tmp2 * fjac[0][4][k][i][j+1]
          - tmp1 * njac[0][4][k][i][j+1];
        lhs[CC][1][4][k][i][j] =  tmp2 * fjac[1][4][k][i][j+1]
          - tmp1 * njac[1][4][k][i][j+1];
        lhs[CC][2][4][k][i][j] =  tmp2 * fjac[2][4][k][i][j+1]
          - tmp1 * njac[2][4][k][i][j+1];
        lhs[CC][3][4][k][i][j] =  tmp2 * fjac[3][4][k][i][j+1]
          - tmp1 * njac[3][4][k][i][j+1];
        lhs[CC][4][4][k][i][j] =  tmp2 * fjac[4][4][k][i][j+1]
          - tmp1 * njac[4][4][k][i][j+1]
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
        matvec_sub(k, i, j, AA, k, j-1, i, k, j, i, lhs, rhs);//y_matvec_sub(lhs[k][i][j][AA], rhs[k][j-1][i], rhs[k][j][i]);

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
      matvec_sub(k, i, jsize, AA, k, jsize-1, i, k, jsize, i, lhs, rhs);//y_matvec_sub(lhs[k][i][jsize][AA], rhs[k][jsize-1][i], rhs[k][jsize][i]);

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
              - lhs[CC][n][m][k][i][j]*rhs[k][j+1][i][n];
          }
        }
      }
    }
  }
  if (timeron) timer_stop(t_ysolve);
}
