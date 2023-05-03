#include "header.h"
#include "work_lhs.h"
#include "timers.h"

//---------------------------------------------------------------------
// Performs line solves in Z direction by first factoring
// the block-tridiagonal matrix into an upper triangular matrix, 
// and then performing back substitution to solve for the unknow
// vectors of each line.  
// 
// Make sure we treat elements zero to cell_size in the direction
// of the sweep.
//---------------------------------------------------------------------
void z_solve()
{
  int i, j, k, m, n, ksize;

  double tmp1, tmp2, tmp3;
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  if (timeron) timer_start(t_zsolve);

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // This function computes the left hand side for the three z-factors   
  //---------------------------------------------------------------------

  ksize = grid_points[2]-1;

  //---------------------------------------------------------------------
  // Compute the indices for storing the block-diagonal matrix;
  // determine c (labeled f) and s jacobians
  //---------------------------------------------------------------------
  #pragma acc parallel loop private(i,j,k,m,n,tmp1,tmp2,tmp3)
  for (j = 1; j <= grid_points[1]-2; j++) {
    //#pragma acc loop
    for (i = 1; i <= grid_points[0]-2; i++) {
      //#pragma acc loop
      for (k = 0; k <= ksize; k++) {
        tmp1 = 1.0 / u[k][j][i][0];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;

        fjac[0][0][j][i][k] = 0.0;
        fjac[1][0][j][i][k] = 0.0;
        fjac[2][0][j][i][k] = 0.0;
        fjac[3][0][j][i][k] = 1.0;
        fjac[4][0][j][i][k] = 0.0;

        fjac[0][1][j][i][k] = - ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[1][1][j][i][k] = u[k][j][i][3] * tmp1;
        fjac[2][1][j][i][k] = 0.0;
        fjac[3][1][j][i][k] = u[k][j][i][1] * tmp1;
        fjac[4][1][j][i][k] = 0.0;

        fjac[0][2][j][i][k] = - ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
        fjac[1][2][j][i][k] = 0.0;
        fjac[2][2][j][i][k] = u[k][j][i][3] * tmp1;
        fjac[3][2][j][i][k] = u[k][j][i][2] * tmp1;
        fjac[4][2][j][i][k] = 0.0;

        fjac[0][3][j][i][k] = - (u[k][j][i][3]*u[k][j][i][3] * tmp2 ) 
          + c2 * qs[k][j][i];
        fjac[1][3][j][i][k] = - c2 *  u[k][j][i][1] * tmp1;
        fjac[2][3][j][i][k] = - c2 *  u[k][j][i][2] * tmp1;
        fjac[3][3][j][i][k] = ( 2.0 - c2 ) *  u[k][j][i][3] * tmp1;
        fjac[4][3][j][i][k] = c2;

        fjac[0][4][j][i][k] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
          * u[k][j][i][3] * tmp2;
        fjac[1][4][j][i][k] = - c2 * ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[2][4][j][i][k] = - c2 * ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
        fjac[3][4][j][i][k] = c1 * ( u[k][j][i][4] * tmp1 )
          - c2 * ( qs[k][j][i] + u[k][j][i][3]*u[k][j][i][3] * tmp2 );
        fjac[4][4][j][i][k] = c1 * u[k][j][i][3] * tmp1;

        njac[0][0][j][i][k] = 0.0;
        njac[1][0][j][i][k] = 0.0;
        njac[2][0][j][i][k] = 0.0;
        njac[3][0][j][i][k] = 0.0;
        njac[4][0][j][i][k] = 0.0;

        njac[0][1][j][i][k] = - c3c4 * tmp2 * u[k][j][i][1];
        njac[1][1][j][i][k] =   c3c4 * tmp1;
        njac[2][1][j][i][k] =   0.0;
        njac[3][1][j][i][k] =   0.0;
        njac[4][1][j][i][k] =   0.0;

        njac[0][2][j][i][k] = - c3c4 * tmp2 * u[k][j][i][2];
        njac[1][2][j][i][k] =   0.0;
        njac[2][2][j][i][k] =   c3c4 * tmp1;
        njac[3][2][j][i][k] =   0.0;
        njac[4][2][j][i][k] =   0.0;

        njac[0][3][j][i][k] = - con43 * c3c4 * tmp2 * u[k][j][i][3];
        njac[1][3][j][i][k] =   0.0;
        njac[2][3][j][i][k] =   0.0;
        njac[3][3][j][i][k] =   con43 * c3 * c4 * tmp1;
        njac[4][3][j][i][k] =   0.0;

        njac[0][4][j][i][k] = - (  c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( con43 * c3c4
              - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[1][4][j][i][k] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][1];
        njac[2][4][j][i][k] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[3][4][j][i][k] = ( con43 * c3c4
            - c1345 ) * tmp2 * u[k][j][i][3];
        njac[4][4][j][i][k] = ( c1345 )* tmp1;
      }

      //---------------------------------------------------------------------
      // now jacobians set, so form left hand side in z direction
      //---------------------------------------------------------------------
      //#pragma acc routine (lhsinit) worker
      lhsinit(j, i, ksize, lhs);//z_lhsinit(lhs[j][i], ksize);
      for (k = 1; k <= ksize-1; k++) {
        tmp1 = dt * tz1;
        tmp2 = dt * tz2;

        lhs[AA][0][0][j][i][k] = - tmp2 * fjac[0][0][j][i][k-1]
          - tmp1 * njac[0][0][j][i][k-1]
          - tmp1 * dz1; 
        lhs[AA][1][0][j][i][k] = - tmp2 * fjac[1][0][j][i][k-1]
          - tmp1 * njac[1][0][j][i][k-1];
        lhs[AA][2][0][j][i][k] = - tmp2 * fjac[2][0][j][i][k-1]
          - tmp1 * njac[2][0][j][i][k-1];
        lhs[AA][3][0][j][i][k] = - tmp2 * fjac[3][0][j][i][k-1]
          - tmp1 * njac[3][0][j][i][k-1];
        lhs[AA][4][0][j][i][k] = - tmp2 * fjac[4][0][j][i][k-1]
          - tmp1 * njac[4][0][j][i][k-1];

        lhs[AA][0][1][j][i][k] = - tmp2 * fjac[0][1][j][i][k-1]
          - tmp1 * njac[0][1][j][i][k-1];
        lhs[AA][1][1][j][i][k] = - tmp2 * fjac[1][1][j][i][k-1]
          - tmp1 * njac[1][1][j][i][k-1]
          - tmp1 * dz2;
        lhs[AA][2][1][j][i][k] = - tmp2 * fjac[2][1][j][i][k-1]
          - tmp1 * njac[2][1][j][i][k-1];
        lhs[AA][3][1][j][i][k] = - tmp2 * fjac[3][1][j][i][k-1]
          - tmp1 * njac[3][1][j][i][k-1];
        lhs[AA][4][1][j][i][k] = - tmp2 * fjac[4][1][j][i][k-1]
          - tmp1 * njac[4][1][j][i][k-1];

        lhs[AA][0][2][j][i][k] = - tmp2 * fjac[0][2][j][i][k-1]
          - tmp1 * njac[0][2][j][i][k-1];
        lhs[AA][1][2][j][i][k] = - tmp2 * fjac[1][2][j][i][k-1]
          - tmp1 * njac[1][2][j][i][k-1];
        lhs[AA][2][2][j][i][k] = - tmp2 * fjac[2][2][j][i][k-1]
          - tmp1 * njac[2][2][j][i][k-1]
          - tmp1 * dz3;
        lhs[AA][3][2][j][i][k] = - tmp2 * fjac[3][2][j][i][k-1]
          - tmp1 * njac[3][2][j][i][k-1];
        lhs[AA][4][2][j][i][k] = - tmp2 * fjac[4][2][j][i][k-1]
          - tmp1 * njac[4][2][j][i][k-1];

        lhs[AA][0][3][j][i][k] = - tmp2 * fjac[0][3][j][i][k-1]
          - tmp1 * njac[0][3][j][i][k-1];
        lhs[AA][1][3][j][i][k] = - tmp2 * fjac[1][3][j][i][k-1]
          - tmp1 * njac[1][3][j][i][k-1];
        lhs[AA][2][3][j][i][k] = - tmp2 * fjac[2][3][j][i][k-1]
          - tmp1 * njac[2][3][j][i][k-1];
        lhs[AA][3][3][j][i][k] = - tmp2 * fjac[3][3][j][i][k-1]
          - tmp1 * njac[3][3][j][i][k-1]
          - tmp1 * dz4;
        lhs[AA][4][3][j][i][k] = - tmp2 * fjac[4][3][j][i][k-1]
          - tmp1 * njac[4][3][j][i][k-1];

        lhs[AA][0][4][j][i][k] = - tmp2 * fjac[0][4][j][i][k-1]
          - tmp1 * njac[0][4][j][i][k-1];
        lhs[AA][1][4][j][i][k] = - tmp2 * fjac[1][4][j][i][k-1]
          - tmp1 * njac[1][4][j][i][k-1];
        lhs[AA][2][4][j][i][k] = - tmp2 * fjac[2][4][j][i][k-1]
          - tmp1 * njac[2][4][j][i][k-1];
        lhs[AA][3][4][j][i][k] = - tmp2 * fjac[3][4][j][i][k-1]
          - tmp1 * njac[3][4][j][i][k-1];
        lhs[AA][4][4][j][i][k] = - tmp2 * fjac[4][4][j][i][k-1]
          - tmp1 * njac[4][4][j][i][k-1]
          - tmp1 * dz5;

        lhs[BB][0][0][j][i][k] = 1.0
          + tmp1 * 2.0 * njac[0][0][j][i][k]
          + tmp1 * 2.0 * dz1;
        lhs[BB][1][0][j][i][k] = tmp1 * 2.0 * njac[1][0][j][i][k];
        lhs[BB][2][0][j][i][k] = tmp1 * 2.0 * njac[2][0][j][i][k];
        lhs[BB][3][0][j][i][k] = tmp1 * 2.0 * njac[3][0][j][i][k];
        lhs[BB][4][0][j][i][k] = tmp1 * 2.0 * njac[4][0][j][i][k];

        lhs[BB][0][1][j][i][k] = tmp1 * 2.0 * njac[0][1][j][i][k];
        lhs[BB][1][1][j][i][k] = 1.0
          + tmp1 * 2.0 * njac[1][1][j][i][k]
          + tmp1 * 2.0 * dz2;
        lhs[BB][2][1][j][i][k] = tmp1 * 2.0 * njac[2][1][j][i][k];
        lhs[BB][3][1][j][i][k] = tmp1 * 2.0 * njac[3][1][j][i][k];
        lhs[BB][4][1][j][i][k] = tmp1 * 2.0 * njac[4][1][j][i][k];

        lhs[BB][0][2][j][i][k] = tmp1 * 2.0 * njac[0][2][j][i][k];
        lhs[BB][1][2][j][i][k] = tmp1 * 2.0 * njac[1][2][j][i][k];
        lhs[BB][2][2][j][i][k] = 1.0
          + tmp1 * 2.0 * njac[2][2][j][i][k]
          + tmp1 * 2.0 * dz3;
        lhs[BB][3][2][j][i][k] = tmp1 * 2.0 * njac[3][2][j][i][k];
        lhs[BB][4][2][j][i][k] = tmp1 * 2.0 * njac[4][2][j][i][k];

        lhs[BB][0][3][j][i][k] = tmp1 * 2.0 * njac[0][3][j][i][k];
        lhs[BB][1][3][j][i][k] = tmp1 * 2.0 * njac[1][3][j][i][k];
        lhs[BB][2][3][j][i][k] = tmp1 * 2.0 * njac[2][3][j][i][k];
        lhs[BB][3][3][j][i][k] = 1.0
          + tmp1 * 2.0 * njac[3][3][j][i][k]
          + tmp1 * 2.0 * dz4;
        lhs[BB][4][3][j][i][k] = tmp1 * 2.0 * njac[4][3][j][i][k];

        lhs[BB][0][4][j][i][k] = tmp1 * 2.0 * njac[0][4][j][i][k];
        lhs[BB][1][4][j][i][k] = tmp1 * 2.0 * njac[1][4][j][i][k];
        lhs[BB][2][4][j][i][k] = tmp1 * 2.0 * njac[2][4][j][i][k];
        lhs[BB][3][4][j][i][k] = tmp1 * 2.0 * njac[3][4][j][i][k];
        lhs[BB][4][4][j][i][k] = 1.0
          + tmp1 * 2.0 * njac[4][4][j][i][k] 
          + tmp1 * 2.0 * dz5;

        lhs[CC][0][0][j][i][k] =  tmp2 * fjac[0][0][j][i][k+1]
          - tmp1 * njac[0][0][j][i][k+1]
          - tmp1 * dz1;
        lhs[CC][1][0][j][i][k] =  tmp2 * fjac[1][0][j][i][k+1]
          - tmp1 * njac[1][0][j][i][k+1];
        lhs[CC][2][0][j][i][k] =  tmp2 * fjac[2][0][j][i][k+1]
          - tmp1 * njac[2][0][j][i][k+1];
        lhs[CC][3][0][j][i][k] =  tmp2 * fjac[3][0][j][i][k+1]
          - tmp1 * njac[3][0][j][i][k+1];
        lhs[CC][4][0][j][i][k] =  tmp2 * fjac[4][0][j][i][k+1]
          - tmp1 * njac[4][0][j][i][k+1];

        lhs[CC][0][1][j][i][k] =  tmp2 * fjac[0][1][j][i][k+1]
          - tmp1 * njac[0][1][j][i][k+1];
        lhs[CC][1][1][j][i][k] =  tmp2 * fjac[1][1][j][i][k+1]
          - tmp1 * njac[1][1][j][i][k+1]
          - tmp1 * dz2;
        lhs[CC][2][1][j][i][k] =  tmp2 * fjac[2][1][j][i][k+1]
          - tmp1 * njac[2][1][j][i][k+1];
        lhs[CC][3][1][j][i][k] =  tmp2 * fjac[3][1][j][i][k+1]
          - tmp1 * njac[3][1][j][i][k+1];
        lhs[CC][4][1][j][i][k] =  tmp2 * fjac[4][1][j][i][k+1]
          - tmp1 * njac[4][1][j][i][k+1];

        lhs[CC][0][2][j][i][k] =  tmp2 * fjac[0][2][j][i][k+1]
          - tmp1 * njac[0][2][j][i][k+1];
        lhs[CC][1][2][j][i][k] =  tmp2 * fjac[1][2][j][i][k+1]
          - tmp1 * njac[1][2][j][i][k+1];
        lhs[CC][2][2][j][i][k] =  tmp2 * fjac[2][2][j][i][k+1]
          - tmp1 * njac[2][2][j][i][k+1]
          - tmp1 * dz3;
        lhs[CC][3][2][j][i][k] =  tmp2 * fjac[3][2][j][i][k+1]
          - tmp1 * njac[3][2][j][i][k+1];
        lhs[CC][4][2][j][i][k] =  tmp2 * fjac[4][2][j][i][k+1]
          - tmp1 * njac[4][2][j][i][k+1];

        lhs[CC][0][3][j][i][k] =  tmp2 * fjac[0][3][j][i][k+1]
          - tmp1 * njac[0][3][j][i][k+1];
        lhs[CC][1][3][j][i][k] =  tmp2 * fjac[1][3][j][i][k+1]
          - tmp1 * njac[1][3][j][i][k+1];
        lhs[CC][2][3][j][i][k] =  tmp2 * fjac[2][3][j][i][k+1]
          - tmp1 * njac[2][3][j][i][k+1];
        lhs[CC][3][3][j][i][k] =  tmp2 * fjac[3][3][j][i][k+1]
          - tmp1 * njac[3][3][j][i][k+1]
          - tmp1 * dz4;
        lhs[CC][4][3][j][i][k] =  tmp2 * fjac[4][3][j][i][k+1]
          - tmp1 * njac[4][3][j][i][k+1];

        lhs[CC][0][4][j][i][k] =  tmp2 * fjac[0][4][j][i][k+1]
          - tmp1 * njac[0][4][j][i][k+1];
        lhs[CC][1][4][j][i][k] =  tmp2 * fjac[1][4][j][i][k+1]
          - tmp1 * njac[1][4][j][i][k+1];
        lhs[CC][2][4][j][i][k] =  tmp2 * fjac[2][4][j][i][k+1]
          - tmp1 * njac[2][4][j][i][k+1];
        lhs[CC][3][4][j][i][k] =  tmp2 * fjac[3][4][j][i][k+1]
          - tmp1 * njac[3][4][j][i][k+1];
        lhs[CC][4][4][j][i][k] =  tmp2 * fjac[4][4][j][i][k+1]
          - tmp1 * njac[4][4][j][i][k+1]
          - tmp1 * dz5;
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
      // c'(KMAX) and rhs'(KMAX) will be sent to next cell.
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // outer most do loops - sweeping in i direction
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // multiply c[0][j][i] by b_inverse and copy back to c
      // multiply rhs(0) by b_inverse(0) and copy to rhs
      //---------------------------------------------------------------------
      //#pragma acc routine (binvcrhs) worker
      binvcrhs( j, i, 0, BB, j, i, 0, CC, 0, j, i, lhs, rhs );//z_binvcrhs( lhs[j][i][0][BB], lhs[j][i][0][CC], rhs[0][j][i] );

      //---------------------------------------------------------------------
      // begin inner most do loop
      // do all the elements of the cell unless last 
      //---------------------------------------------------------------------
      for (k = 1; k <= ksize-1; k++) {
        //-------------------------------------------------------------------
        // subtract A*lhs_vector(k-1) from lhs_vector(k)
        // 
        // rhs(k) = rhs(k) - A*rhs(k-1)
        //-------------------------------------------------------------------
        //#pragma acc routine (matvec_sub) worker
        matvec_sub(j, i, k, AA, k-1, j, i, k, j, i, lhs, rhs);//z_matvec_sub(lhs[j][i][k][AA], rhs[k-1][j][i], rhs[k][j][i]);

        //-------------------------------------------------------------------
        // B(k) = B(k) - C(k-1)*A(k)
        // matmul_sub(AA,i,j,k,c,CC,i,j,k-1,c,BB,i,j,k)
        //-------------------------------------------------------------------
        //#pragma acc routine (matmul_sub) worker
        matmul_sub(j, i, k, AA, j, i, k-1, CC, j, i, k, BB, lhs);//z_matmul_sub(lhs[j][i][k][AA], lhs[j][i][k-1][CC], lhs[j][i][k][BB]);

        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[0][j][i] by b_inverse[0][j][i] and copy to rhs
        //-------------------------------------------------------------------
        //#pragma acc routine (binvcrhs) worker
        binvcrhs( j, i, k, BB, j, i, k, CC, k, j, i, lhs, rhs );//z_binvcrhs( lhs[j][i][k][BB], lhs[j][i][k][CC], rhs[k][j][i] );
      }

      //---------------------------------------------------------------------
      // Now finish up special cases for last cell
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
      //---------------------------------------------------------------------
      //#pragma acc routine (matvec_sub) worker
      matvec_sub(j, i, ksize, AA, ksize-1, j, i, ksize, j, i, lhs, rhs);//z_matvec_sub(lhs[j][i][ksize][AA], rhs[ksize-1][j][i], rhs[ksize][j][i]);

      //---------------------------------------------------------------------
      // B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
      // matmul_sub(AA,i,j,ksize,c,
      // $              CC,i,j,ksize-1,c,BB,i,j,ksize)
      //---------------------------------------------------------------------
      //#pragma acc routine (matmul_sub) worker
      matmul_sub(j, i, ksize, AA, j, i, ksize-1, CC, j, i, ksize, BB, lhs);//z_matmul_sub(lhs[j][i][ksize][AA], lhs[j][i][ksize-1][CC], lhs[j][i][ksize][BB]);

      //---------------------------------------------------------------------
      // multiply rhs(ksize) by b_inverse(ksize) and copy to rhs
      //---------------------------------------------------------------------
      //#pragma acc routine (binvrhs) worker
      binvrhs( j, i, ksize, BB, ksize, j, i, lhs, rhs );//z_binvrhs( lhs[j][i][ksize][BB], rhs[ksize][j][i] );

      //---------------------------------------------------------------------
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // back solve: if last cell, then generate U(ksize)=rhs(ksize)
      // else assume U(ksize) is loaded in un pack backsub_info
      // so just use it
      // after u(kstart) will be sent to next cell
      //---------------------------------------------------------------------

      for (k = ksize-1; k >= 0; k--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] 
              - lhs[CC][n][m][j][i][k]*rhs[k+1][j][i][n];
          }
        }
      }
    }
  }
  if (timeron) timer_stop(t_zsolve);
}
