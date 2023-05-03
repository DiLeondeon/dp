#include "header.h"
#include "work_lhs.h"
#include "timers.h"
#include <stdio.h>

//---------------------------------------------------------------------
// 
// Performs line solves in X direction by first factoring
// the block-tridiagonal matrix into an upper triangular matrix, 
// and then performing back substitution to solve for the unknow
// vectors of each line.  
// 
// Make sure we treat elements zero to cell_size in the direction
// of the sweep.
// 
//---------------------------------------------------------------------
void x_solve()
{
  int i, j, k, m, n, isize;

  double tmp1, tmp2, tmp3;
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  if (timeron) timer_start(t_xsolve);

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // This function computes the left hand side in the xi-direction
  //---------------------------------------------------------------------

  isize = grid_points[0]-1;

  //---------------------------------------------------------------------
  // determine a (labeled f) and n jacobians
  //---------------------------------------------------------------------
  #pragma acc parallel loop collapse(2) private(i,j,k,m,n,tmp1,tmp2,tmp3)
  for (k = 1; k <= grid_points[2]-2; k++) {
    //#pragma acc loop 
    for (j = 1; j <= grid_points[1]-2; j++) {
      //#pragma acc loop
      for (i = 0; i <= isize; i++) {
        tmp1 = rho_i[k][j][i];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;
        //-------------------------------------------------------------------
        // 
        //-------------------------------------------------------------------
        fjac[0][0][k][j][i] = 0.0;
        fjac[1][0][k][j][i] = 1.0;
        fjac[2][0][k][j][i] = 0.0;
        fjac[3][0][k][j][i] = 0.0;
        fjac[4][0][k][j][i] = 0.0;

        fjac[0][1][k][j][i] = -(u[k][j][i][1] * tmp2 * u[k][j][i][1])
          + c2 * qs[k][j][i];
        fjac[1][1][k][j][i] = ( 2.0 - c2 ) * ( u[k][j][i][1] / u[k][j][i][0] );
        fjac[2][1][k][j][i] = - c2 * ( u[k][j][i][2] * tmp1 );
        fjac[3][1][k][j][i] = - c2 * ( u[k][j][i][3] * tmp1 );
        fjac[4][1][k][j][i] = c2;

        fjac[0][2][k][j][i] = - ( u[k][j][i][1]*u[k][j][i][2] ) * tmp2;
        fjac[1][2][k][j][i] = u[k][j][i][2] * tmp1;
        fjac[2][2][k][j][i] = u[k][j][i][1] * tmp1;
        fjac[3][2][k][j][i] = 0.0;
        fjac[4][2][k][j][i] = 0.0;

        fjac[0][3][k][j][i] = - ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[1][3][k][j][i] = u[k][j][i][3] * tmp1;
        fjac[2][3][k][j][i] = 0.0;
        fjac[3][3][k][j][i] = u[k][j][i][1] * tmp1;
        fjac[4][3][k][j][i] = 0.0;

        fjac[0][4][k][j][i] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
          * ( u[k][j][i][1] * tmp2 );
        fjac[1][4][k][j][i] = c1 *  u[k][j][i][4] * tmp1 
          - c2 * ( u[k][j][i][1]*u[k][j][i][1] * tmp2 + qs[k][j][i] );
        fjac[2][4][k][j][i] = - c2 * ( u[k][j][i][2]*u[k][j][i][1] ) * tmp2;
        fjac[3][4][k][j][i] = - c2 * ( u[k][j][i][3]*u[k][j][i][1] ) * tmp2;
        fjac[4][4][k][j][i] = c1 * ( u[k][j][i][1] * tmp1 );

        njac[0][0][k][j][i] = 0.0;
        njac[1][0][k][j][i] = 0.0;
        njac[2][0][k][j][i] = 0.0;
        njac[3][0][k][j][i] = 0.0;
        njac[4][0][k][j][i] = 0.0;

        njac[0][1][k][j][i] = - con43 * c3c4 * tmp2 * u[k][j][i][1];
        njac[1][1][k][j][i] =   con43 * c3c4 * tmp1;
        njac[2][1][k][j][i] =   0.0;
        njac[3][1][k][j][i] =   0.0;
        njac[4][1][k][j][i] =   0.0;

        njac[0][2][k][j][i] = - c3c4 * tmp2 * u[k][j][i][2];
        njac[1][2][k][j][i] =   0.0;
        njac[2][2][k][j][i] =   c3c4 * tmp1;
        njac[3][2][k][j][i] =   0.0;
        njac[4][2][k][j][i] =   0.0;

        njac[0][3][k][j][i] = - c3c4 * tmp2 * u[k][j][i][3];
        njac[1][3][k][j][i] =   0.0;
        njac[2][3][k][j][i] =   0.0;
        njac[3][3][k][j][i] =   c3c4 * tmp1;
        njac[4][3][k][j][i] =   0.0;

        njac[0][4][k][j][i] = - ( con43 * c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[1][4][k][j][i] = ( con43 * c3c4
            - c1345 ) * tmp2 * u[k][j][i][1];
        njac[2][4][k][j][i] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[3][4][k][j][i] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][3];
        njac[4][4][k][j][i] = ( c1345 ) * tmp1;
      }
      //---------------------------------------------------------------------
      // now jacobians set, so form left hand side in x direction
      //---------------------------------------------------------------------
      //#pragma acc routine (lhsinit) worker
      lhsinit(k, j, isize, lhs);//x_lhsinit(lhs[k][j], isize);
      for (i = 1; i <= isize-1; i++) {
        tmp1 = dt * tx1;
        tmp2 = dt * tx2;

        lhs[AA][0][0][k][j][i] = - tmp2 * fjac[0][0][k][j][i-1]
          - tmp1 * njac[0][0][k][j][i-1]
          - tmp1 * dx1; 
        lhs[AA][1][0][k][j][i] = - tmp2 * fjac[1][0][k][j][i-1]
          - tmp1 * njac[1][0][k][j][i-1];
        lhs[AA][2][0][k][j][i] = - tmp2 * fjac[2][0][k][j][i-1]
          - tmp1 * njac[2][0][k][j][i-1];
        lhs[AA][3][0][k][j][i] = - tmp2 * fjac[3][0][k][j][i-1]
          - tmp1 * njac[3][0][k][j][i-1];
        lhs[AA][4][0][k][j][i] = - tmp2 * fjac[4][0][k][j][i-1]
          - tmp1 * njac[4][0][k][j][i-1];

        lhs[AA][0][1][k][j][i] = - tmp2 * fjac[0][1][k][j][i-1]
          - tmp1 * njac[0][1][k][j][i-1];
        lhs[AA][1][1][k][j][i] = - tmp2 * fjac[1][1][k][j][i-1]
          - tmp1 * njac[1][1][k][j][i-1]
          - tmp1 * dx2;
        lhs[AA][2][1][k][j][i] = - tmp2 * fjac[2][1][k][j][i-1]
          - tmp1 * njac[2][1][k][j][i-1];
        lhs[AA][3][1][k][j][i] = - tmp2 * fjac[3][1][k][j][i-1]
          - tmp1 * njac[3][1][k][j][i-1];
        lhs[AA][4][1][k][j][i] = - tmp2 * fjac[4][1][k][j][i-1]
          - tmp1 * njac[4][1][k][j][i-1];

        lhs[AA][0][2][k][j][i] = - tmp2 * fjac[0][2][k][j][i-1]
          - tmp1 * njac[0][2][k][j][i-1];
        lhs[AA][1][2][k][j][i] = - tmp2 * fjac[1][2][k][j][i-1]
          - tmp1 * njac[1][2][k][j][i-1];
        lhs[AA][2][2][k][j][i] = - tmp2 * fjac[2][2][k][j][i-1]
          - tmp1 * njac[2][2][k][j][i-1]
          - tmp1 * dx3;
        lhs[AA][3][2][k][j][i] = - tmp2 * fjac[3][2][k][j][i-1]
          - tmp1 * njac[3][2][k][j][i-1];
        lhs[AA][4][2][k][j][i] = - tmp2 * fjac[4][2][k][j][i-1]
          - tmp1 * njac[4][2][k][j][i-1];

        lhs[AA][0][3][k][j][i] = - tmp2 * fjac[0][3][k][j][i-1]
          - tmp1 * njac[0][3][k][j][i-1];
        lhs[AA][1][3][k][j][i] = - tmp2 * fjac[1][3][k][j][i-1]
          - tmp1 * njac[1][3][k][j][i-1];
        lhs[AA][2][3][k][j][i] = - tmp2 * fjac[2][3][k][j][i-1]
          - tmp1 * njac[2][3][k][j][i-1];
        lhs[AA][3][3][k][j][i] = - tmp2 * fjac[3][3][k][j][i-1]
          - tmp1 * njac[3][3][k][j][i-1]
          - tmp1 * dx4;
        lhs[AA][4][3][k][j][i] = - tmp2 * fjac[4][3][k][j][i-1]
          - tmp1 * njac[4][3][k][j][i-1];

        lhs[AA][0][4][k][j][i] = - tmp2 * fjac[0][4][k][j][i-1]
          - tmp1 * njac[0][4][k][j][i-1];
        lhs[AA][1][4][k][j][i] = - tmp2 * fjac[1][4][k][j][i-1]
          - tmp1 * njac[1][4][k][j][i-1];
        lhs[AA][2][4][k][j][i] = - tmp2 * fjac[2][4][k][j][i-1]
          - tmp1 * njac[2][4][k][j][i-1];
        lhs[AA][3][4][k][j][i] = - tmp2 * fjac[3][4][k][j][i-1]
          - tmp1 * njac[3][4][k][j][i-1];
        lhs[AA][4][4][k][j][i] = - tmp2 * fjac[4][4][k][j][i-1]
          - tmp1 * njac[4][4][k][j][i-1]
          - tmp1 * dx5;

        lhs[BB][0][0][k][j][i] = 1.0
          + tmp1 * 2.0 * njac[0][0][k][j][i]
          + tmp1 * 2.0 * dx1;
        lhs[BB][1][0][k][j][i] = tmp1 * 2.0 * njac[1][0][k][j][i];
        lhs[BB][2][0][k][j][i] = tmp1 * 2.0 * njac[2][0][k][j][i];
        lhs[BB][3][0][k][j][i] = tmp1 * 2.0 * njac[3][0][k][j][i];
        lhs[BB][4][0][k][j][i] = tmp1 * 2.0 * njac[4][0][k][j][i];

        lhs[BB][0][1][k][j][i] = tmp1 * 2.0 * njac[0][1][k][j][i];
        lhs[BB][1][1][k][j][i] = 1.0
          + tmp1 * 2.0 * njac[1][1][k][j][i]
          + tmp1 * 2.0 * dx2;
        lhs[BB][2][1][k][j][i] = tmp1 * 2.0 * njac[2][1][k][j][i];
        lhs[BB][3][1][k][j][i] = tmp1 * 2.0 * njac[3][1][k][j][i];
        lhs[BB][4][1][k][j][i] = tmp1 * 2.0 * njac[4][1][k][j][i];

        lhs[BB][0][2][k][j][i] = tmp1 * 2.0 * njac[0][2][k][j][i];
        lhs[BB][1][2][k][j][i] = tmp1 * 2.0 * njac[1][2][k][j][i];
        lhs[BB][2][2][k][j][i] = 1.0
          + tmp1 * 2.0 * njac[2][2][k][j][i]
          + tmp1 * 2.0 * dx3;
        lhs[BB][3][2][k][j][i] = tmp1 * 2.0 * njac[3][2][k][j][i];
        lhs[BB][4][2][k][j][i] = tmp1 * 2.0 * njac[4][2][k][j][i];

        lhs[BB][0][3][k][j][i] = tmp1 * 2.0 * njac[0][3][k][j][i];
        lhs[BB][1][3][k][j][i] = tmp1 * 2.0 * njac[1][3][k][j][i];
        lhs[BB][2][3][k][j][i] = tmp1 * 2.0 * njac[2][3][k][j][i];
        lhs[BB][3][3][k][j][i] = 1.0
          + tmp1 * 2.0 * njac[3][3][k][j][i]
          + tmp1 * 2.0 * dx4;
        lhs[BB][4][3][k][j][i] = tmp1 * 2.0 * njac[4][3][k][j][i];

        lhs[BB][0][4][k][j][i] = tmp1 * 2.0 * njac[0][4][k][j][i];
        lhs[BB][1][4][k][j][i] = tmp1 * 2.0 * njac[1][4][k][j][i];
        lhs[BB][2][4][k][j][i] = tmp1 * 2.0 * njac[2][4][k][j][i];
        lhs[BB][3][4][k][j][i] = tmp1 * 2.0 * njac[3][4][k][j][i];
        lhs[BB][4][4][k][j][i] = 1.0
          + tmp1 * 2.0 * njac[4][4][k][j][i]
          + tmp1 * 2.0 * dx5;

        lhs[CC][0][0][k][j][i] =  tmp2 * fjac[0][0][k][j][i+1]
          - tmp1 * njac[0][0][k][j][i+1]
          - tmp1 * dx1;
        lhs[CC][1][0][k][j][i] =  tmp2 * fjac[1][0][k][j][i+1]
          - tmp1 * njac[1][0][k][j][i+1];
        lhs[CC][2][0][k][j][i] =  tmp2 * fjac[2][0][k][j][i+1]
          - tmp1 * njac[2][0][k][j][i+1];
        lhs[CC][3][0][k][j][i] =  tmp2 * fjac[3][0][k][j][i+1]
          - tmp1 * njac[3][0][k][j][i+1];
        lhs[CC][4][0][k][j][i] =  tmp2 * fjac[4][0][k][j][i+1]
          - tmp1 * njac[4][0][k][j][i+1];

        lhs[CC][0][1][k][j][i] =  tmp2 * fjac[0][1][k][j][i+1]
          - tmp1 * njac[0][1][k][j][i+1];
        lhs[CC][1][1][k][j][i] =  tmp2 * fjac[1][1][k][j][i+1]
          - tmp1 * njac[1][1][k][j][i+1]
          - tmp1 * dx2;
        lhs[CC][2][1][k][j][i] =  tmp2 * fjac[2][1][k][j][i+1]
          - tmp1 * njac[2][1][k][j][i+1];
        lhs[CC][3][1][k][j][i] =  tmp2 * fjac[3][1][k][j][i+1]
          - tmp1 * njac[3][1][k][j][i+1];
        lhs[CC][4][1][k][j][i] =  tmp2 * fjac[4][1][k][j][i+1]
          - tmp1 * njac[4][1][k][j][i+1];

        lhs[CC][0][2][k][j][i] =  tmp2 * fjac[0][2][k][j][i+1]
          - tmp1 * njac[0][2][k][j][i+1];
        lhs[CC][1][2][k][j][i] =  tmp2 * fjac[1][2][k][j][i+1]
          - tmp1 * njac[1][2][k][j][i+1];
        lhs[CC][2][2][k][j][i] =  tmp2 * fjac[2][2][k][j][i+1]
          - tmp1 * njac[2][2][k][j][i+1]
          - tmp1 * dx3;
        lhs[CC][3][2][k][j][i] =  tmp2 * fjac[3][2][k][j][i+1]
          - tmp1 * njac[3][2][k][j][i+1];
        lhs[CC][4][2][k][j][i] =  tmp2 * fjac[4][2][k][j][i+1]
          - tmp1 * njac[4][2][k][j][i+1];

        lhs[CC][0][3][k][j][i] =  tmp2 * fjac[0][3][k][j][i+1]
          - tmp1 * njac[0][3][k][j][i+1];
        lhs[CC][1][3][k][j][i] =  tmp2 * fjac[1][3][k][j][i+1]
          - tmp1 * njac[1][3][k][j][i+1];
        lhs[CC][2][3][k][j][i] =  tmp2 * fjac[2][3][k][j][i+1]
          - tmp1 * njac[2][3][k][j][i+1];
        lhs[CC][3][3][k][j][i] =  tmp2 * fjac[3][3][k][j][i+1]
          - tmp1 * njac[3][3][k][j][i+1]
          - tmp1 * dx4;
        lhs[CC][4][3][k][j][i] =  tmp2 * fjac[4][3][k][j][i+1]
          - tmp1 * njac[4][3][k][j][i+1];

        lhs[CC][0][4][k][j][i] =  tmp2 * fjac[0][4][k][j][i+1]
          - tmp1 * njac[0][4][k][j][i+1];
        lhs[CC][1][4][k][j][i] =  tmp2 * fjac[1][4][k][j][i+1]
          - tmp1 * njac[1][4][k][j][i+1];
        lhs[CC][2][4][k][j][i] =  tmp2 * fjac[2][4][k][j][i+1]
          - tmp1 * njac[2][4][k][j][i+1];
        lhs[CC][3][4][k][j][i] =  tmp2 * fjac[3][4][k][j][i+1]
          - tmp1 * njac[3][4][k][j][i+1];
        lhs[CC][4][4][k][j][i] =  tmp2 * fjac[4][4][k][j][i+1]
          - tmp1 * njac[4][4][k][j][i+1]
          - tmp1 * dx5;
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
      // c'(IMAX) and rhs'(IMAX) will be sent to next cell
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // outer most do loops - sweeping in i direction
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // multiply c[k][j][0] by b_inverse and copy back to c
      // multiply rhs(0) by b_inverse(0) and copy to rhs
      //---------------------------------------------------------------------
      //#pragma acc routine (binvcrhs) worker
      binvcrhs( k, j, 0, BB, k, j, 0, CC, k, j, 0, lhs, rhs );//x_binvcrhs( lhs[k][j][0][BB], lhs[k][j][0][CC], rhs[k][j][0] );

      //---------------------------------------------------------------------
      // begin inner most do loop
      // do all the elements of the cell unless last 
      //---------------------------------------------------------------------
      for (i = 1; i <= isize-1; i++) {
        //-------------------------------------------------------------------
        // rhs(i) = rhs(i) - A*rhs(i-1)
        //-------------------------------------------------------------------
        //#pragma acc routine (matvec_sub) worker
        matvec_sub(k, j, i, AA, k, j, i-1, k, j, i, lhs, rhs);//x_matvec_sub(lhs[k][j][i][AA], rhs[k][j][i-1], rhs[k][j][i]);

        //-------------------------------------------------------------------
        // B(i) = B(i) - C(i-1)*A(i)
        //-------------------------------------------------------------------
        //#pragma acc routine (matmul_sub) worker
        matmul_sub(k, j, i, AA, k, j, i-1, CC, k, j, i, BB, lhs);//x_matmul_sub(lhs[k][j][i][AA], lhs[k][j][i-1][CC], lhs[k][j][i][BB]);


        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[k][j][0] by b_inverse[k][j][0] and copy to rhs
        //-------------------------------------------------------------------
        //#pragma acc routine (binvcrhs) worker
        binvcrhs( k, j, i, BB, k, j, i, CC, k, j, i, lhs, rhs );//x_binvcrhs( lhs[k][j][i][BB], lhs[k][j][i][CC], rhs[k][j][i] );
      }

      //---------------------------------------------------------------------
      // rhs(isize) = rhs(isize) - A*rhs(isize-1)
      //---------------------------------------------------------------------
      //#pragma acc routine (matvec_sub) worker
      matvec_sub(k, j, isize, AA, k, j, isize-1, k, j, isize, lhs, rhs);//x_matvec_sub(lhs[k][j][isize][AA], rhs[k][j][isize-1], rhs[k][j][isize]);

      //---------------------------------------------------------------------
      // B(isize) = B(isize) - C(isize-1)*A(isize)
      //---------------------------------------------------------------------
      //#pragma acc routine (matmul_sub) worker
      matmul_sub(k, j, isize, AA, k, j, isize-1, CC, k, j, isize, BB, lhs);//x_matmul_sub(lhs[k][j][isize][AA], lhs[k][j][isize-1][CC], lhs[k][j][isize][BB]);

      //---------------------------------------------------------------------
      // multiply rhs() by b_inverse() and copy to rhs
      //---------------------------------------------------------------------
      //#pragma acc routine (binvrhs) worker
      binvrhs( k, j, isize, BB, k, j, isize, lhs, rhs );//x_binvrhs( lhs[k][j][isize][BB], rhs[k][j][isize] );

      //---------------------------------------------------------------------
      // back solve: if last cell, then generate U(isize)=rhs(isize)
      // else assume U(isize) is loaded in un pack backsub_info
      // so just use it
      // after u(istart) will be sent to next cell
      //---------------------------------------------------------------------
      for (i = isize-1; i >=0; i--) {
        for (m = 0; m < BLOCK_SIZE; m++) {
          for (n = 0; n < BLOCK_SIZE; n++) {
            rhs[k][j][i][m] = rhs[k][j][i][m] 
              - lhs[CC][n][m][k][j][i]*rhs[k][j][i+1][n];
          }
        }
      }
    }
  }
  if (timeron) timer_stop(t_xsolve);
}
