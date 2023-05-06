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
  int i, j, k, m, n, m1, isize, jsize, ksize;

  double tmp1, tmp2, tmp3;
  double pivot, coeff;
  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  if (timeron) timer_start(t_zsolve);

  //---------------------------------------------------------------------
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // This function computes the left hand side for the three z-factors   
  //---------------------------------------------------------------------
  isize = grid_points[0]-2;
  ksize = grid_points[2]-1;
  jsize = grid_points[1]-2;
  //---------------------------------------------------------------------
  // Compute the indices for storing the block-diagonal matrix;
  // determine c (labeled f) and s jacobians
  //---------------------------------------------------------------------
  #pragma acc parallel loop collapse(2) private(i,j,k,m,m1,n,tmp1,tmp2,tmp3,pivot,coeff)
  for (j = 1; j <= jsize; j++) {
    //#pragma acc loop
    for (i = 1; i <= isize; i++) {
      //#pragma acc loop
      for (k = 0; k <= ksize; k++) {
        tmp1 = 1.0 / u[k][j][i][0];
        tmp2 = tmp1 * tmp1;
        tmp3 = tmp1 * tmp2;

        fjac[j][i][k][0][0] = 0.0;
        fjac[j][i][k][1][0] = 0.0;
        fjac[j][i][k][2][0] = 0.0;
        fjac[j][i][k][3][0] = 1.0;
        fjac[j][i][k][4][0] = 0.0;

        fjac[j][i][k][0][1] = - ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[j][i][k][1][1] = u[k][j][i][3] * tmp1;
        fjac[j][i][k][2][1] = 0.0;
        fjac[j][i][k][3][1] = u[k][j][i][1] * tmp1;
        fjac[j][i][k][4][1] = 0.0;

        fjac[j][i][k][0][2] = - ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
        fjac[j][i][k][1][2] = 0.0;
        fjac[j][i][k][2][2] = u[k][j][i][3] * tmp1;
        fjac[j][i][k][3][2] = u[k][j][i][2] * tmp1;
        fjac[j][i][k][4][2] = 0.0;

        fjac[j][i][k][0][3] = - (u[k][j][i][3]*u[k][j][i][3] * tmp2 ) 
          + c2 * qs[k][j][i];
        fjac[j][i][k][1][3] = - c2 *  u[k][j][i][1] * tmp1;
        fjac[j][i][k][2][3] = - c2 *  u[k][j][i][2] * tmp1;
        fjac[j][i][k][3][3] = ( 2.0 - c2 ) *  u[k][j][i][3] * tmp1;
        fjac[j][i][k][4][3] = c2;

        fjac[j][i][k][0][4] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
          * u[k][j][i][3] * tmp2;
        fjac[j][i][k][1][4] = - c2 * ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[j][i][k][2][4] = - c2 * ( u[k][j][i][2]*u[k][j][i][3] ) * tmp2;
        fjac[j][i][k][3][4] = c1 * ( u[k][j][i][4] * tmp1 )
          - c2 * ( qs[k][j][i] + u[k][j][i][3]*u[k][j][i][3] * tmp2 );
        fjac[j][i][k][4][4] = c1 * u[k][j][i][3] * tmp1;

        njac[j][i][k][0][0] = 0.0;
        njac[j][i][k][1][0] = 0.0;
        njac[j][i][k][2][0] = 0.0;
        njac[j][i][k][3][0] = 0.0;
        njac[j][i][k][4][0] = 0.0;

        njac[j][i][k][0][1] = - c3c4 * tmp2 * u[k][j][i][1];
        njac[j][i][k][1][1] =   c3c4 * tmp1;
        njac[j][i][k][2][1] =   0.0;
        njac[j][i][k][3][1] =   0.0;
        njac[j][i][k][4][1] =   0.0;

        njac[j][i][k][0][2] = - c3c4 * tmp2 * u[k][j][i][2];
        njac[j][i][k][1][2] =   0.0;
        njac[j][i][k][2][2] =   c3c4 * tmp1;
        njac[j][i][k][3][2] =   0.0;
        njac[j][i][k][4][2] =   0.0;

        njac[j][i][k][0][3] = - con43 * c3c4 * tmp2 * u[k][j][i][3];
        njac[j][i][k][1][3] =   0.0;
        njac[j][i][k][2][3] =   0.0;
        njac[j][i][k][3][3] =   con43 * c3 * c4 * tmp1;
        njac[j][i][k][4][3] =   0.0;

        njac[j][i][k][0][4] = - (  c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( con43 * c3c4
              - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[j][i][k][1][4] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][1];
        njac[j][i][k][2][4] = (  c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[j][i][k][3][4] = ( con43 * c3c4
            - c1345 ) * tmp2 * u[k][j][i][3];
        njac[j][i][k][4][4] = ( c1345 )* tmp1;
      }

      //---------------------------------------------------------------------
      // now jacobians set, so form left hand side in z direction
      //---------------------------------------------------------------------
      //#pragma acc routine (lhsinit) worker
      //lhsinit(j, i, ksize, lhs);//z_lhsinit(lhs[j][i], ksize);
      for (n = 0; n < 5; n++) {
        for (m = 0; m < 5; m++) {
          lhs[j][i][0][0][n][m] = 0.0;
          lhs[j][i][0][1][n][m] = 0.0;
          lhs[j][i][0][2][n][m] = 0.0;
          lhs[j][i][ksize][0][n][m] = 0.0;
          lhs[j][i][ksize][1][n][m] = 0.0;
          lhs[j][i][ksize][2][n][m] = 0.0;
        }
        lhs[j][i][0][1][n][n] = 1.0;
        lhs[j][i][ksize][1][n][n] = 1.0;
      }
      for (k = 1; k <= ksize-1; k++) {
        tmp1 = dt * tz1;
        tmp2 = dt * tz2;

        lhs[j][i][k][AA][0][0] = - tmp2 * fjac[j][i][k-1][0][0]
          - tmp1 * njac[j][i][k-1][0][0]
          - tmp1 * dz1; 
        lhs[j][i][k][AA][1][0] = - tmp2 * fjac[j][i][k-1][1][0]
          - tmp1 * njac[j][i][k-1][1][0];
        lhs[j][i][k][AA][2][0] = - tmp2 * fjac[j][i][k-1][2][0]
          - tmp1 * njac[j][i][k-1][2][0];
        lhs[j][i][k][AA][3][0] = - tmp2 * fjac[j][i][k-1][3][0]
          - tmp1 * njac[j][i][k-1][3][0];
        lhs[j][i][k][AA][4][0] = - tmp2 * fjac[j][i][k-1][4][0]
          - tmp1 * njac[j][i][k-1][4][0];

        lhs[j][i][k][AA][0][1] = - tmp2 * fjac[j][i][k-1][0][1]
          - tmp1 * njac[j][i][k-1][0][1];
        lhs[j][i][k][AA][1][1] = - tmp2 * fjac[j][i][k-1][1][1]
          - tmp1 * njac[j][i][k-1][1][1]
          - tmp1 * dz2;
        lhs[j][i][k][AA][2][1] = - tmp2 * fjac[j][i][k-1][2][1]
          - tmp1 * njac[j][i][k-1][2][1];
        lhs[j][i][k][AA][3][1] = - tmp2 * fjac[j][i][k-1][3][1]
          - tmp1 * njac[j][i][k-1][3][1];
        lhs[j][i][k][AA][4][1] = - tmp2 * fjac[j][i][k-1][4][1]
          - tmp1 * njac[j][i][k-1][4][1];

        lhs[j][i][k][AA][0][2] = - tmp2 * fjac[j][i][k-1][0][2]
          - tmp1 * njac[j][i][k-1][0][2];
        lhs[j][i][k][AA][1][2] = - tmp2 * fjac[j][i][k-1][1][2]
          - tmp1 * njac[j][i][k-1][1][2];
        lhs[j][i][k][AA][2][2] = - tmp2 * fjac[j][i][k-1][2][2]
          - tmp1 * njac[j][i][k-1][2][2]
          - tmp1 * dz3;
        lhs[j][i][k][AA][3][2] = - tmp2 * fjac[j][i][k-1][3][2]
          - tmp1 * njac[j][i][k-1][3][2];
        lhs[j][i][k][AA][4][2] = - tmp2 * fjac[j][i][k-1][4][2]
          - tmp1 * njac[j][i][k-1][4][2];

        lhs[j][i][k][AA][0][3] = - tmp2 * fjac[j][i][k-1][0][3]
          - tmp1 * njac[j][i][k-1][0][3];
        lhs[j][i][k][AA][1][3] = - tmp2 * fjac[j][i][k-1][1][3]
          - tmp1 * njac[j][i][k-1][1][3];
        lhs[j][i][k][AA][2][3] = - tmp2 * fjac[j][i][k-1][2][3]
          - tmp1 * njac[j][i][k-1][2][3];
        lhs[j][i][k][AA][3][3] = - tmp2 * fjac[j][i][k-1][3][3]
          - tmp1 * njac[j][i][k-1][3][3]
          - tmp1 * dz4;
        lhs[j][i][k][AA][4][3] = - tmp2 * fjac[j][i][k-1][4][3]
          - tmp1 * njac[j][i][k-1][4][3];

        lhs[j][i][k][AA][0][4] = - tmp2 * fjac[j][i][k-1][0][4]
          - tmp1 * njac[j][i][k-1][0][4];
        lhs[j][i][k][AA][1][4] = - tmp2 * fjac[j][i][k-1][1][4]
          - tmp1 * njac[j][i][k-1][1][4];
        lhs[j][i][k][AA][2][4] = - tmp2 * fjac[j][i][k-1][2][4]
          - tmp1 * njac[j][i][k-1][2][4];
        lhs[j][i][k][AA][3][4] = - tmp2 * fjac[j][i][k-1][3][4]
          - tmp1 * njac[j][i][k-1][3][4];
        lhs[j][i][k][AA][4][4] = - tmp2 * fjac[j][i][k-1][4][4]
          - tmp1 * njac[j][i][k-1][4][4]
          - tmp1 * dz5;

        lhs[j][i][k][BB][0][0] = 1.0
          + tmp1 * 2.0 * njac[j][i][k][0][0]
          + tmp1 * 2.0 * dz1;
        lhs[j][i][k][BB][1][0] = tmp1 * 2.0 * njac[j][i][k][1][0];
        lhs[j][i][k][BB][2][0] = tmp1 * 2.0 * njac[j][i][k][2][0];
        lhs[j][i][k][BB][3][0] = tmp1 * 2.0 * njac[j][i][k][3][0];
        lhs[j][i][k][BB][4][0] = tmp1 * 2.0 * njac[j][i][k][4][0];

        lhs[j][i][k][BB][0][1] = tmp1 * 2.0 * njac[j][i][k][0][1];
        lhs[j][i][k][BB][1][1] = 1.0
          + tmp1 * 2.0 * njac[j][i][k][1][1]
          + tmp1 * 2.0 * dz2;
        lhs[j][i][k][BB][2][1] = tmp1 * 2.0 * njac[j][i][k][2][1];
        lhs[j][i][k][BB][3][1] = tmp1 * 2.0 * njac[j][i][k][3][1];
        lhs[j][i][k][BB][4][1] = tmp1 * 2.0 * njac[j][i][k][4][1];

        lhs[j][i][k][BB][0][2] = tmp1 * 2.0 * njac[j][i][k][0][2];
        lhs[j][i][k][BB][1][2] = tmp1 * 2.0 * njac[j][i][k][1][2];
        lhs[j][i][k][BB][2][2] = 1.0
          + tmp1 * 2.0 * njac[j][i][k][2][2]
          + tmp1 * 2.0 * dz3;
        lhs[j][i][k][BB][3][2] = tmp1 * 2.0 * njac[j][i][k][3][2];
        lhs[j][i][k][BB][4][2] = tmp1 * 2.0 * njac[j][i][k][4][2];

        lhs[j][i][k][BB][0][3] = tmp1 * 2.0 * njac[j][i][k][0][3];
        lhs[j][i][k][BB][1][3] = tmp1 * 2.0 * njac[j][i][k][1][3];
        lhs[j][i][k][BB][2][3] = tmp1 * 2.0 * njac[j][i][k][2][3];
        lhs[j][i][k][BB][3][3] = 1.0
          + tmp1 * 2.0 * njac[j][i][k][3][3]
          + tmp1 * 2.0 * dz4;
        lhs[j][i][k][BB][4][3] = tmp1 * 2.0 * njac[j][i][k][4][3];

        lhs[j][i][k][BB][0][4] = tmp1 * 2.0 * njac[j][i][k][0][4];
        lhs[j][i][k][BB][1][4] = tmp1 * 2.0 * njac[j][i][k][1][4];
        lhs[j][i][k][BB][2][4] = tmp1 * 2.0 * njac[j][i][k][2][4];
        lhs[j][i][k][BB][3][4] = tmp1 * 2.0 * njac[j][i][k][3][4];
        lhs[j][i][k][BB][4][4] = 1.0
          + tmp1 * 2.0 * njac[j][i][k][4][4] 
          + tmp1 * 2.0 * dz5;

        lhs[j][i][k][CC][0][0] =  tmp2 * fjac[j][i][k+1][0][0]
          - tmp1 * njac[j][i][k+1][0][0]
          - tmp1 * dz1;
        lhs[j][i][k][CC][1][0] =  tmp2 * fjac[j][i][k+1][1][0]
          - tmp1 * njac[j][i][k+1][1][0];
        lhs[j][i][k][CC][2][0] =  tmp2 * fjac[j][i][k+1][2][0]
          - tmp1 * njac[j][i][k+1][2][0];
        lhs[j][i][k][CC][3][0] =  tmp2 * fjac[j][i][k+1][3][0]
          - tmp1 * njac[j][i][k+1][3][0];
        lhs[j][i][k][CC][4][0] =  tmp2 * fjac[j][i][k+1][4][0]
          - tmp1 * njac[j][i][k+1][4][0];

        lhs[j][i][k][CC][0][1] =  tmp2 * fjac[j][i][k+1][0][1]
          - tmp1 * njac[j][i][k+1][0][1];
        lhs[j][i][k][CC][1][1] =  tmp2 * fjac[j][i][k+1][1][1]
          - tmp1 * njac[j][i][k+1][1][1]
          - tmp1 * dz2;
        lhs[j][i][k][CC][2][1] =  tmp2 * fjac[j][i][k+1][2][1]
          - tmp1 * njac[j][i][k+1][2][1];
        lhs[j][i][k][CC][3][1] =  tmp2 * fjac[j][i][k+1][3][1]
          - tmp1 * njac[j][i][k+1][3][1];
        lhs[j][i][k][CC][4][1] =  tmp2 * fjac[j][i][k+1][4][1]
          - tmp1 * njac[j][i][k+1][4][1];

        lhs[j][i][k][CC][0][2] =  tmp2 * fjac[j][i][k+1][0][2]
          - tmp1 * njac[j][i][k+1][0][2];
        lhs[j][i][k][CC][1][2] =  tmp2 * fjac[j][i][k+1][1][2]
          - tmp1 * njac[j][i][k+1][1][2];
        lhs[j][i][k][CC][2][2] =  tmp2 * fjac[j][i][k+1][2][2]
          - tmp1 * njac[j][i][k+1][2][2]
          - tmp1 * dz3;
        lhs[j][i][k][CC][3][2] =  tmp2 * fjac[j][i][k+1][3][2]
          - tmp1 * njac[j][i][k+1][3][2];
        lhs[j][i][k][CC][4][2] =  tmp2 * fjac[j][i][k+1][4][2]
          - tmp1 * njac[j][i][k+1][4][2];

        lhs[j][i][k][CC][0][3] =  tmp2 * fjac[j][i][k+1][0][3]
          - tmp1 * njac[j][i][k+1][0][3];
        lhs[j][i][k][CC][1][3] =  tmp2 * fjac[j][i][k+1][1][3]
          - tmp1 * njac[j][i][k+1][1][3];
        lhs[j][i][k][CC][2][3] =  tmp2 * fjac[j][i][k+1][2][3]
          - tmp1 * njac[j][i][k+1][2][3];
        lhs[j][i][k][CC][3][3] =  tmp2 * fjac[j][i][k+1][3][3]
          - tmp1 * njac[j][i][k+1][3][3]
          - tmp1 * dz4;
        lhs[j][i][k][CC][4][3] =  tmp2 * fjac[j][i][k+1][4][3]
          - tmp1 * njac[j][i][k+1][4][3];

        lhs[j][i][k][CC][0][4] =  tmp2 * fjac[j][i][k+1][0][4]
          - tmp1 * njac[j][i][k+1][0][4];
        lhs[j][i][k][CC][1][4] =  tmp2 * fjac[j][i][k+1][1][4]
          - tmp1 * njac[j][i][k+1][1][4];
        lhs[j][i][k][CC][2][4] =  tmp2 * fjac[j][i][k+1][2][4]
          - tmp1 * njac[j][i][k+1][2][4];
        lhs[j][i][k][CC][3][4] =  tmp2 * fjac[j][i][k+1][3][4]
          - tmp1 * njac[j][i][k+1][3][4];
        lhs[j][i][k][CC][4][4] =  tmp2 * fjac[j][i][k+1][4][4]
          - tmp1 * njac[j][i][k+1][4][4]
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
      //binvcrhs( j, i, 0, BB, j, i, 0, CC, 0, j, i, lhs, rhs );//z_binvcrhs( lhs[j][i][0][BB], lhs[j][i][0][CC], rhs[0][j][i] );
      for (m=0;m<5;m++){
        pivot = 1.00/lhs[j][i][0][BB][m][m];
        for (n=m+1;n<5;n++){
          lhs[j][i][0][BB][n][m] = lhs[j][i][0][BB][n][m]*pivot;
        }
        for (n=0;n<5;n++){
          lhs[j][i][0][CC][n][m] = lhs[j][i][0][CC][n][m]*pivot;
        }
        rhs[0][j][i][m]   = rhs[0][j][i][m]  *pivot;
        for (m1=0;m1<5;m1++){
          if (m1 != m){
            coeff = lhs[j][i][0][BB][m][m1];
            for (n=m+1;n<5;n++){
              lhs[j][i][0][BB][n][m1]= lhs[j][i][0][BB][n][m1] - coeff*lhs[j][i][0][BB][n][m];
            }
            for (n=0;n<5;n++){
              lhs[j][i][0][CC][n][m1] = lhs[j][i][0][CC][n][m1] - coeff*lhs[j][i][0][CC][n][m];
            }
            rhs[0][j][i][m1]   = rhs[0][j][i][m1]   - coeff*rhs[0][j][i][m];
          }
        }
      }

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
        //matvec_sub(j, i, k, AA, k-1, j, i, k, j, i, lhs, rhs);//z_matvec_sub(lhs[j][i][k][AA], rhs[k-1][j][i], rhs[k][j][i]);
        rhs[k][j][i][0] = rhs[k][j][i][0] - lhs[j][i][k][AA][0][0]*rhs[k-1][j][i][0]
                    - lhs[j][i][k][AA][1][0]*rhs[k-1][j][i][1]
                    - lhs[j][i][k][AA][2][0]*rhs[k-1][j][i][2]
                    - lhs[j][i][k][AA][3][0]*rhs[k-1][j][i][3]
                    - lhs[j][i][k][AA][4][0]*rhs[k-1][j][i][4];
        rhs[k][j][i][1] = rhs[k][j][i][1] - lhs[j][i][k][AA][0][1]*rhs[k-1][j][i][0]
                    - lhs[j][i][k][AA][1][1]*rhs[k-1][j][i][1]
                    - lhs[j][i][k][AA][2][1]*rhs[k-1][j][i][2]
                    - lhs[j][i][k][AA][3][1]*rhs[k-1][j][i][3]
                    - lhs[j][i][k][AA][4][1]*rhs[k-1][j][i][4];
        rhs[k][j][i][2] = rhs[k][j][i][2] - lhs[j][i][k][AA][0][2]*rhs[k-1][j][i][0]
                    - lhs[j][i][k][AA][1][2]*rhs[k-1][j][i][1]
                    - lhs[j][i][k][AA][2][2]*rhs[k-1][j][i][2]
                    - lhs[j][i][k][AA][3][2]*rhs[k-1][j][i][3]
                    - lhs[j][i][k][AA][4][2]*rhs[k-1][j][i][4];
        rhs[k][j][i][3] = rhs[k][j][i][3] - lhs[j][i][k][AA][0][3]*rhs[k-1][j][i][0]
                    - lhs[j][i][k][AA][1][3]*rhs[k-1][j][i][1]
                    - lhs[j][i][k][AA][2][3]*rhs[k-1][j][i][2]
                    - lhs[j][i][k][AA][3][3]*rhs[k-1][j][i][3]
                    - lhs[j][i][k][AA][4][3]*rhs[k-1][j][i][4];
        rhs[k][j][i][4] = rhs[k][j][i][4] - lhs[j][i][k][AA][0][4]*rhs[k-1][j][i][0]
                    - lhs[j][i][k][AA][1][4]*rhs[k-1][j][i][1]
                    - lhs[j][i][k][AA][2][4]*rhs[k-1][j][i][2]
                    - lhs[j][i][k][AA][3][4]*rhs[k-1][j][i][3]
                    - lhs[j][i][k][AA][4][4]*rhs[k-1][j][i][4];

        //-------------------------------------------------------------------
        // B(k) = B(k) - C(k-1)*A(k)
        // matmul_sub(AA,i,j,k,c,CC,i,j,k-1,c,BB,i,j,k)
        //-------------------------------------------------------------------
        //#pragma acc routine (matmul_sub) worker
        //matmul_sub(j, i, k, AA, j, i, k-1, CC, j, i, k, BB, lhs);//z_matmul_sub(lhs[j][i][k][AA], lhs[j][i][k-1][CC], lhs[j][i][k][BB]);
                  //j, i, ksize, AA, j, i, ksize-1, CC, j, i, ksize, BB, lhs
        lhs[j][i][k][BB][0][0] = lhs[j][i][k][BB][0][0] - lhs[j][i][k][AA][0][0]*lhs[j][i][k-1][CC][0][0]
                              - lhs[j][i][k][AA][1][0]*lhs[j][i][k-1][CC][0][1]
                              - lhs[j][i][k][AA][2][0]*lhs[j][i][k-1][CC][0][2]
                              - lhs[j][i][k][AA][3][0]*lhs[j][i][k-1][CC][0][3]
                              - lhs[j][i][k][AA][4][0]*lhs[j][i][k-1][CC][0][4];
        lhs[j][i][k][BB][0][1] = lhs[j][i][k][BB][0][1] - lhs[j][i][k][AA][0][1]*lhs[j][i][k-1][CC][0][0]
                              - lhs[j][i][k][AA][1][1]*lhs[j][i][k-1][CC][0][1]
                              - lhs[j][i][k][AA][2][1]*lhs[j][i][k-1][CC][0][2]
                              - lhs[j][i][k][AA][3][1]*lhs[j][i][k-1][CC][0][3]
                              - lhs[j][i][k][AA][4][1]*lhs[j][i][k-1][CC][0][4];
        lhs[j][i][k][BB][0][2] = lhs[j][i][k][BB][0][2] - lhs[j][i][k][AA][0][2]*lhs[j][i][k-1][CC][0][0]
                              - lhs[j][i][k][AA][1][2]*lhs[j][i][k-1][CC][0][1]
                              - lhs[j][i][k][AA][2][2]*lhs[j][i][k-1][CC][0][2]
                              - lhs[j][i][k][AA][3][2]*lhs[j][i][k-1][CC][0][3]
                              - lhs[j][i][k][AA][4][2]*lhs[j][i][k-1][CC][0][4];
        lhs[j][i][k][BB][0][3] = lhs[j][i][k][BB][0][3] - lhs[j][i][k][AA][0][3]*lhs[j][i][k-1][CC][0][0]
                              - lhs[j][i][k][AA][1][3]*lhs[j][i][k-1][CC][0][1]
                              - lhs[j][i][k][AA][2][3]*lhs[j][i][k-1][CC][0][2]
                              - lhs[j][i][k][AA][3][3]*lhs[j][i][k-1][CC][0][3]
                              - lhs[j][i][k][AA][4][3]*lhs[j][i][k-1][CC][0][4];
        lhs[j][i][k][BB][0][4] = lhs[j][i][k][BB][0][4] - lhs[j][i][k][AA][0][4]*lhs[j][i][k-1][CC][0][0]
                              - lhs[j][i][k][AA][1][4]*lhs[j][i][k-1][CC][0][1]
                              - lhs[j][i][k][AA][2][4]*lhs[j][i][k-1][CC][0][2]
                              - lhs[j][i][k][AA][3][4]*lhs[j][i][k-1][CC][0][3]
                              - lhs[j][i][k][AA][4][4]*lhs[j][i][k-1][CC][0][4];
        lhs[j][i][k][BB][1][0] = lhs[j][i][k][BB][1][0] - lhs[j][i][k][AA][0][0]*lhs[j][i][k-1][CC][1][0]
                              - lhs[j][i][k][AA][1][0]*lhs[j][i][k-1][CC][1][1]
                              - lhs[j][i][k][AA][2][0]*lhs[j][i][k-1][CC][1][2]
                              - lhs[j][i][k][AA][3][0]*lhs[j][i][k-1][CC][1][3]
                              - lhs[j][i][k][AA][4][0]*lhs[j][i][k-1][CC][1][4];
        lhs[j][i][k][BB][1][1] = lhs[j][i][k][BB][1][1] - lhs[j][i][k][AA][0][1]*lhs[j][i][k-1][CC][1][0]
                              - lhs[j][i][k][AA][1][1]*lhs[j][i][k-1][CC][1][1]
                              - lhs[j][i][k][AA][2][1]*lhs[j][i][k-1][CC][1][2]
                              - lhs[j][i][k][AA][3][1]*lhs[j][i][k-1][CC][1][3]
                              - lhs[j][i][k][AA][4][1]*lhs[j][i][k-1][CC][1][4];
        lhs[j][i][k][BB][1][2] = lhs[j][i][k][BB][1][2] - lhs[j][i][k][AA][0][2]*lhs[j][i][k-1][CC][1][0]
                              - lhs[j][i][k][AA][1][2]*lhs[j][i][k-1][CC][1][1]
                              - lhs[j][i][k][AA][2][2]*lhs[j][i][k-1][CC][1][2]
                              - lhs[j][i][k][AA][3][2]*lhs[j][i][k-1][CC][1][3]
                              - lhs[j][i][k][AA][4][2]*lhs[j][i][k-1][CC][1][4];
        lhs[j][i][k][BB][1][3] = lhs[j][i][k][BB][1][3] - lhs[j][i][k][AA][0][3]*lhs[j][i][k-1][CC][1][0]
                              - lhs[j][i][k][AA][1][3]*lhs[j][i][k-1][CC][1][1]
                              - lhs[j][i][k][AA][2][3]*lhs[j][i][k-1][CC][1][2]
                              - lhs[j][i][k][AA][3][3]*lhs[j][i][k-1][CC][1][3]
                              - lhs[j][i][k][AA][4][3]*lhs[j][i][k-1][CC][1][4];
        lhs[j][i][k][BB][1][4] = lhs[j][i][k][BB][1][4] - lhs[j][i][k][AA][0][4]*lhs[j][i][k-1][CC][1][0]
                              - lhs[j][i][k][AA][1][4]*lhs[j][i][k-1][CC][1][1]
                              - lhs[j][i][k][AA][2][4]*lhs[j][i][k-1][CC][1][2]
                              - lhs[j][i][k][AA][3][4]*lhs[j][i][k-1][CC][1][3]
                              - lhs[j][i][k][AA][4][4]*lhs[j][i][k-1][CC][1][4];
        lhs[j][i][k][BB][2][0] = lhs[j][i][k][BB][2][0] - lhs[j][i][k][AA][0][0]*lhs[j][i][k-1][CC][2][0]
                              - lhs[j][i][k][AA][1][0]*lhs[j][i][k-1][CC][2][1]
                              - lhs[j][i][k][AA][2][0]*lhs[j][i][k-1][CC][2][2]
                              - lhs[j][i][k][AA][3][0]*lhs[j][i][k-1][CC][2][3]
                              - lhs[j][i][k][AA][4][0]*lhs[j][i][k-1][CC][2][4];
        lhs[j][i][k][BB][2][1] = lhs[j][i][k][BB][2][1] - lhs[j][i][k][AA][0][1]*lhs[j][i][k-1][CC][2][0]
                              - lhs[j][i][k][AA][1][1]*lhs[j][i][k-1][CC][2][1]
                              - lhs[j][i][k][AA][2][1]*lhs[j][i][k-1][CC][2][2]
                              - lhs[j][i][k][AA][3][1]*lhs[j][i][k-1][CC][2][3]
                              - lhs[j][i][k][AA][4][1]*lhs[j][i][k-1][CC][2][4];
        lhs[j][i][k][BB][2][2] = lhs[j][i][k][BB][2][2] - lhs[j][i][k][AA][0][2]*lhs[j][i][k-1][CC][2][0]
                              - lhs[j][i][k][AA][1][2]*lhs[j][i][k-1][CC][2][1]
                              - lhs[j][i][k][AA][2][2]*lhs[j][i][k-1][CC][2][2]
                              - lhs[j][i][k][AA][3][2]*lhs[j][i][k-1][CC][2][3]
                              - lhs[j][i][k][AA][4][2]*lhs[j][i][k-1][CC][2][4];
        lhs[j][i][k][BB][2][3] = lhs[j][i][k][BB][2][3] - lhs[j][i][k][AA][0][3]*lhs[j][i][k-1][CC][2][0]
                              - lhs[j][i][k][AA][1][3]*lhs[j][i][k-1][CC][2][1]
                              - lhs[j][i][k][AA][2][3]*lhs[j][i][k-1][CC][2][2]
                              - lhs[j][i][k][AA][3][3]*lhs[j][i][k-1][CC][2][3]
                              - lhs[j][i][k][AA][4][3]*lhs[j][i][k-1][CC][2][4];
        lhs[j][i][k][BB][2][4] = lhs[j][i][k][BB][2][4] - lhs[j][i][k][AA][0][4]*lhs[j][i][k-1][CC][2][0]
                              - lhs[j][i][k][AA][1][4]*lhs[j][i][k-1][CC][2][1]
                              - lhs[j][i][k][AA][2][4]*lhs[j][i][k-1][CC][2][2]
                              - lhs[j][i][k][AA][3][4]*lhs[j][i][k-1][CC][2][3]
                              - lhs[j][i][k][AA][4][4]*lhs[j][i][k-1][CC][2][4];
        lhs[j][i][k][BB][3][0] = lhs[j][i][k][BB][3][0] - lhs[j][i][k][AA][0][0]*lhs[j][i][k-1][CC][3][0]
                              - lhs[j][i][k][AA][1][0]*lhs[j][i][k-1][CC][3][1]
                              - lhs[j][i][k][AA][2][0]*lhs[j][i][k-1][CC][3][2]
                              - lhs[j][i][k][AA][3][0]*lhs[j][i][k-1][CC][3][3]
                              - lhs[j][i][k][AA][4][0]*lhs[j][i][k-1][CC][3][4];
        lhs[j][i][k][BB][3][1] = lhs[j][i][k][BB][3][1] - lhs[j][i][k][AA][0][1]*lhs[j][i][k-1][CC][3][0]
                              - lhs[j][i][k][AA][1][1]*lhs[j][i][k-1][CC][3][1]
                              - lhs[j][i][k][AA][2][1]*lhs[j][i][k-1][CC][3][2]
                              - lhs[j][i][k][AA][3][1]*lhs[j][i][k-1][CC][3][3]
                              - lhs[j][i][k][AA][4][1]*lhs[j][i][k-1][CC][3][4];
        lhs[j][i][k][BB][3][2] = lhs[j][i][k][BB][3][2] - lhs[j][i][k][AA][0][2]*lhs[j][i][k-1][CC][3][0]
                              - lhs[j][i][k][AA][1][2]*lhs[j][i][k-1][CC][3][1]
                              - lhs[j][i][k][AA][2][2]*lhs[j][i][k-1][CC][3][2]
                              - lhs[j][i][k][AA][3][2]*lhs[j][i][k-1][CC][3][3]
                              - lhs[j][i][k][AA][4][2]*lhs[j][i][k-1][CC][3][4];
        lhs[j][i][k][BB][3][3] = lhs[j][i][k][BB][3][3] - lhs[j][i][k][AA][0][3]*lhs[j][i][k-1][CC][3][0]
                              - lhs[j][i][k][AA][1][3]*lhs[j][i][k-1][CC][3][1]
                              - lhs[j][i][k][AA][2][3]*lhs[j][i][k-1][CC][3][2]
                              - lhs[j][i][k][AA][3][3]*lhs[j][i][k-1][CC][3][3]
                              - lhs[j][i][k][AA][4][3]*lhs[j][i][k-1][CC][3][4];
        lhs[j][i][k][BB][3][4] = lhs[j][i][k][BB][3][4] - lhs[j][i][k][AA][0][4]*lhs[j][i][k-1][CC][3][0]
                              - lhs[j][i][k][AA][1][4]*lhs[j][i][k-1][CC][3][1]
                              - lhs[j][i][k][AA][2][4]*lhs[j][i][k-1][CC][3][2]
                              - lhs[j][i][k][AA][3][4]*lhs[j][i][k-1][CC][3][3]
                              - lhs[j][i][k][AA][4][4]*lhs[j][i][k-1][CC][3][4];
        lhs[j][i][k][BB][4][0] = lhs[j][i][k][BB][4][0] - lhs[j][i][k][AA][0][0]*lhs[j][i][k-1][CC][4][0]
                              - lhs[j][i][k][AA][1][0]*lhs[j][i][k-1][CC][4][1]
                              - lhs[j][i][k][AA][2][0]*lhs[j][i][k-1][CC][4][2]
                              - lhs[j][i][k][AA][3][0]*lhs[j][i][k-1][CC][4][3]
                              - lhs[j][i][k][AA][4][0]*lhs[j][i][k-1][CC][4][4];
        lhs[j][i][k][BB][4][1] = lhs[j][i][k][BB][4][1] - lhs[j][i][k][AA][0][1]*lhs[j][i][k-1][CC][4][0]
                              - lhs[j][i][k][AA][1][1]*lhs[j][i][k-1][CC][4][1]
                              - lhs[j][i][k][AA][2][1]*lhs[j][i][k-1][CC][4][2]
                              - lhs[j][i][k][AA][3][1]*lhs[j][i][k-1][CC][4][3]
                              - lhs[j][i][k][AA][4][1]*lhs[j][i][k-1][CC][4][4];
        lhs[j][i][k][BB][4][2] = lhs[j][i][k][BB][4][2] - lhs[j][i][k][AA][0][2]*lhs[j][i][k-1][CC][4][0]
                              - lhs[j][i][k][AA][1][2]*lhs[j][i][k-1][CC][4][1]
                              - lhs[j][i][k][AA][2][2]*lhs[j][i][k-1][CC][4][2]
                              - lhs[j][i][k][AA][3][2]*lhs[j][i][k-1][CC][4][3]
                              - lhs[j][i][k][AA][4][2]*lhs[j][i][k-1][CC][4][4];
        lhs[j][i][k][BB][4][3] = lhs[j][i][k][BB][4][3] - lhs[j][i][k][AA][0][3]*lhs[j][i][k-1][CC][4][0]
                              - lhs[j][i][k][AA][1][3]*lhs[j][i][k-1][CC][4][1]
                              - lhs[j][i][k][AA][2][3]*lhs[j][i][k-1][CC][4][2]
                              - lhs[j][i][k][AA][3][3]*lhs[j][i][k-1][CC][4][3]
                              - lhs[j][i][k][AA][4][3]*lhs[j][i][k-1][CC][4][4];
        lhs[j][i][k][BB][4][4] = lhs[j][i][k][BB][4][4] - lhs[j][i][k][AA][0][4]*lhs[j][i][k-1][CC][4][0]
                              - lhs[j][i][k][AA][1][4]*lhs[j][i][k-1][CC][4][1]
                              - lhs[j][i][k][AA][2][4]*lhs[j][i][k-1][CC][4][2]
                              - lhs[j][i][k][AA][3][4]*lhs[j][i][k-1][CC][4][3]
                              - lhs[j][i][k][AA][4][4]*lhs[j][i][k-1][CC][4][4];

        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[0][j][i] by b_inverse[0][j][i] and copy to rhs
        //-------------------------------------------------------------------
        //#pragma acc routine (binvcrhs) worker
        //binvcrhs( j, i, k, BB, j, i, k, CC, k, j, i, lhs, rhs );//z_binvcrhs( lhs[j][i][k][BB], lhs[j][i][k][CC], rhs[k][j][i] );
        for (m=0;m<5;m++){
          pivot = 1.00/lhs[j][i][k][BB][m][m];
          for (n=m+1;n<5;n++){
            lhs[j][i][k][BB][n][m] = lhs[j][i][k][BB][n][m]*pivot;
          }
          for (n=0;n<5;n++){
            lhs[j][i][k][CC][n][m] = lhs[j][i][k][CC][n][m]*pivot;
          }
          rhs[k][j][i][m]   = rhs[k][j][i][m]  *pivot;
          for (m1=0;m1<5;m1++){
            if (m1 != m){
              coeff = lhs[j][i][k][BB][m][m1];
              for (n=m+1;n<5;n++){
                lhs[j][i][k][BB][n][m1]= lhs[j][i][k][BB][n][m1] - coeff*lhs[j][i][k][BB][n][m];
              }
              for (n=0;n<5;n++){
                lhs[j][i][k][CC][n][m1] = lhs[j][i][k][CC][n][m1] - coeff*lhs[j][i][k][CC][n][m];
              }
              rhs[k][j][i][m1]   = rhs[k][j][i][m1]   - coeff*rhs[k][j][i][m];
            }
          }
        }                      
      }

      //---------------------------------------------------------------------
      // Now finish up special cases for last cell
      //---------------------------------------------------------------------

      //---------------------------------------------------------------------
      // rhs(ksize) = rhs(ksize) - A*rhs(ksize-1)
      //---------------------------------------------------------------------
      //#pragma acc routine (matvec_sub) worker
      //matvec_sub(j, i, ksize, AA, ksize-1, j, i, ksize, j, i, lhs, rhs);//z_matvec_sub(lhs[j][i][ksize][AA], rhs[ksize-1][j][i], rhs[ksize][j][i]);
      rhs[ksize][j][i][0] = rhs[ksize][j][i][0] - lhs[j][i][ksize][AA][0][0]*rhs[ksize-1][j][i][0]
                    - lhs[j][i][ksize][AA][1][0]*rhs[ksize-1][j][i][1]
                    - lhs[j][i][ksize][AA][2][0]*rhs[ksize-1][j][i][2]
                    - lhs[j][i][ksize][AA][3][0]*rhs[ksize-1][j][i][3]
                    - lhs[j][i][ksize][AA][4][0]*rhs[ksize-1][j][i][4];
      rhs[ksize][j][i][1] = rhs[ksize][j][i][1] - lhs[j][i][ksize][AA][0][1]*rhs[ksize-1][j][i][0]
                    - lhs[j][i][ksize][AA][1][1]*rhs[ksize-1][j][i][1]
                    - lhs[j][i][ksize][AA][2][1]*rhs[ksize-1][j][i][2]
                    - lhs[j][i][ksize][AA][3][1]*rhs[ksize-1][j][i][3]
                    - lhs[j][i][ksize][AA][4][1]*rhs[ksize-1][j][i][4];
      rhs[ksize][j][i][2] = rhs[ksize][j][i][2] - lhs[j][i][ksize][AA][0][2]*rhs[ksize-1][j][i][0]
                    - lhs[j][i][ksize][AA][1][2]*rhs[ksize-1][j][i][1]
                    - lhs[j][i][ksize][AA][2][2]*rhs[ksize-1][j][i][2]
                    - lhs[j][i][ksize][AA][3][2]*rhs[ksize-1][j][i][3]
                    - lhs[j][i][ksize][AA][4][2]*rhs[ksize-1][j][i][4];
      rhs[ksize][j][i][3] = rhs[ksize][j][i][3] - lhs[j][i][ksize][AA][0][3]*rhs[ksize-1][j][i][0]
                    - lhs[j][i][ksize][AA][1][3]*rhs[ksize-1][j][i][1]
                    - lhs[j][i][ksize][AA][2][3]*rhs[ksize-1][j][i][2]
                    - lhs[j][i][ksize][AA][3][3]*rhs[ksize-1][j][i][3]
                    - lhs[j][i][ksize][AA][4][3]*rhs[ksize-1][j][i][4];
      rhs[ksize][j][i][4] = rhs[ksize][j][i][4] - lhs[j][i][ksize][AA][0][4]*rhs[ksize-1][j][i][0]
                    - lhs[j][i][ksize][AA][1][4]*rhs[ksize-1][j][i][1]
                    - lhs[j][i][ksize][AA][2][4]*rhs[ksize-1][j][i][2]
                    - lhs[j][i][ksize][AA][3][4]*rhs[ksize-1][j][i][3]
                    - lhs[j][i][ksize][AA][4][4]*rhs[ksize-1][j][i][4];

      //---------------------------------------------------------------------
      // B(ksize) = B(ksize) - C(ksize-1)*A(ksize)
      // matmul_sub(AA,i,j,ksize,c,
      // $              CC,i,j,ksize-1,c,BB,i,j,ksize)
      //---------------------------------------------------------------------
      //#pragma acc routine (matmul_sub) worker
      //matmul_sub(j, i, ksize, AA, j, i, ksize-1, CC, j, i, ksize, BB, lhs);//z_matmul_sub(lhs[j][i][ksize][AA], lhs[j][i][ksize-1][CC], lhs[j][i][ksize][BB]);
      lhs[j][i][ksize][BB][0][0] = lhs[j][i][ksize][BB][0][0] - lhs[j][i][ksize][AA][0][0]*lhs[j][i][ksize-1][CC][0][0]
                              - lhs[j][i][ksize][AA][1][0]*lhs[j][i][ksize-1][CC][0][1]
                              - lhs[j][i][ksize][AA][2][0]*lhs[j][i][ksize-1][CC][0][2]
                              - lhs[j][i][ksize][AA][3][0]*lhs[j][i][ksize-1][CC][0][3]
                              - lhs[j][i][ksize][AA][4][0]*lhs[j][i][ksize-1][CC][0][4];
      lhs[j][i][ksize][BB][0][1] = lhs[j][i][ksize][BB][0][1] - lhs[j][i][ksize][AA][0][1]*lhs[j][i][ksize-1][CC][0][0]
                              - lhs[j][i][ksize][AA][1][1]*lhs[j][i][ksize-1][CC][0][1]
                              - lhs[j][i][ksize][AA][2][1]*lhs[j][i][ksize-1][CC][0][2]
                              - lhs[j][i][ksize][AA][3][1]*lhs[j][i][ksize-1][CC][0][3]
                              - lhs[j][i][ksize][AA][4][1]*lhs[j][i][ksize-1][CC][0][4];
      lhs[j][i][ksize][BB][0][2] = lhs[j][i][ksize][BB][0][2] - lhs[j][i][ksize][AA][0][2]*lhs[j][i][ksize-1][CC][0][0]
                              - lhs[j][i][ksize][AA][1][2]*lhs[j][i][ksize-1][CC][0][1]
                              - lhs[j][i][ksize][AA][2][2]*lhs[j][i][ksize-1][CC][0][2]
                              - lhs[j][i][ksize][AA][3][2]*lhs[j][i][ksize-1][CC][0][3]
                              - lhs[j][i][ksize][AA][4][2]*lhs[j][i][ksize-1][CC][0][4];
      lhs[j][i][ksize][BB][0][3] = lhs[j][i][ksize][BB][0][3] - lhs[j][i][ksize][AA][0][3]*lhs[j][i][ksize-1][CC][0][0]
                              - lhs[j][i][ksize][AA][1][3]*lhs[j][i][ksize-1][CC][0][1]
                              - lhs[j][i][ksize][AA][2][3]*lhs[j][i][ksize-1][CC][0][2]
                              - lhs[j][i][ksize][AA][3][3]*lhs[j][i][ksize-1][CC][0][3]
                              - lhs[j][i][ksize][AA][4][3]*lhs[j][i][ksize-1][CC][0][4];
      lhs[j][i][ksize][BB][0][4] = lhs[j][i][ksize][BB][0][4] - lhs[j][i][ksize][AA][0][4]*lhs[j][i][ksize-1][CC][0][0]
                              - lhs[j][i][ksize][AA][1][4]*lhs[j][i][ksize-1][CC][0][1]
                              - lhs[j][i][ksize][AA][2][4]*lhs[j][i][ksize-1][CC][0][2]
                              - lhs[j][i][ksize][AA][3][4]*lhs[j][i][ksize-1][CC][0][3]
                              - lhs[j][i][ksize][AA][4][4]*lhs[j][i][ksize-1][CC][0][4];
      lhs[j][i][ksize][BB][1][0] = lhs[j][i][ksize][BB][1][0] - lhs[j][i][ksize][AA][0][0]*lhs[j][i][ksize-1][CC][1][0]
                              - lhs[j][i][ksize][AA][1][0]*lhs[j][i][ksize-1][CC][1][1]
                              - lhs[j][i][ksize][AA][2][0]*lhs[j][i][ksize-1][CC][1][2]
                              - lhs[j][i][ksize][AA][3][0]*lhs[j][i][ksize-1][CC][1][3]
                              - lhs[j][i][ksize][AA][4][0]*lhs[j][i][ksize-1][CC][1][4];
      lhs[j][i][ksize][BB][1][1] = lhs[j][i][ksize][BB][1][1] - lhs[j][i][ksize][AA][0][1]*lhs[j][i][ksize-1][CC][1][0]
                              - lhs[j][i][ksize][AA][1][1]*lhs[j][i][ksize-1][CC][1][1]
                              - lhs[j][i][ksize][AA][2][1]*lhs[j][i][ksize-1][CC][1][2]
                              - lhs[j][i][ksize][AA][3][1]*lhs[j][i][ksize-1][CC][1][3]
                              - lhs[j][i][ksize][AA][4][1]*lhs[j][i][ksize-1][CC][1][4];
      lhs[j][i][ksize][BB][1][2] = lhs[j][i][ksize][BB][1][2] - lhs[j][i][ksize][AA][0][2]*lhs[j][i][ksize-1][CC][1][0]
                              - lhs[j][i][ksize][AA][1][2]*lhs[j][i][ksize-1][CC][1][1]
                              - lhs[j][i][ksize][AA][2][2]*lhs[j][i][ksize-1][CC][1][2]
                              - lhs[j][i][ksize][AA][3][2]*lhs[j][i][ksize-1][CC][1][3]
                              - lhs[j][i][ksize][AA][4][2]*lhs[j][i][ksize-1][CC][1][4];
      lhs[j][i][ksize][BB][1][3] = lhs[j][i][ksize][BB][1][3] - lhs[j][i][ksize][AA][0][3]*lhs[j][i][ksize-1][CC][1][0]
                              - lhs[j][i][ksize][AA][1][3]*lhs[j][i][ksize-1][CC][1][1]
                              - lhs[j][i][ksize][AA][2][3]*lhs[j][i][ksize-1][CC][1][2]
                              - lhs[j][i][ksize][AA][3][3]*lhs[j][i][ksize-1][CC][1][3]
                              - lhs[j][i][ksize][AA][4][3]*lhs[j][i][ksize-1][CC][1][4];
      lhs[j][i][ksize][BB][1][4] = lhs[j][i][ksize][BB][1][4] - lhs[j][i][ksize][AA][0][4]*lhs[j][i][ksize-1][CC][1][0]
                              - lhs[j][i][ksize][AA][1][4]*lhs[j][i][ksize-1][CC][1][1]
                              - lhs[j][i][ksize][AA][2][4]*lhs[j][i][ksize-1][CC][1][2]
                              - lhs[j][i][ksize][AA][3][4]*lhs[j][i][ksize-1][CC][1][3]
                              - lhs[j][i][ksize][AA][4][4]*lhs[j][i][ksize-1][CC][1][4];
      lhs[j][i][ksize][BB][2][0] = lhs[j][i][ksize][BB][2][0] - lhs[j][i][ksize][AA][0][0]*lhs[j][i][ksize-1][CC][2][0]
                              - lhs[j][i][ksize][AA][1][0]*lhs[j][i][ksize-1][CC][2][1]
                              - lhs[j][i][ksize][AA][2][0]*lhs[j][i][ksize-1][CC][2][2]
                              - lhs[j][i][ksize][AA][3][0]*lhs[j][i][ksize-1][CC][2][3]
                              - lhs[j][i][ksize][AA][4][0]*lhs[j][i][ksize-1][CC][2][4];
      lhs[j][i][ksize][BB][2][1] = lhs[j][i][ksize][BB][2][1] - lhs[j][i][ksize][AA][0][1]*lhs[j][i][ksize-1][CC][2][0]
                              - lhs[j][i][ksize][AA][1][1]*lhs[j][i][ksize-1][CC][2][1]
                              - lhs[j][i][ksize][AA][2][1]*lhs[j][i][ksize-1][CC][2][2]
                              - lhs[j][i][ksize][AA][3][1]*lhs[j][i][ksize-1][CC][2][3]
                              - lhs[j][i][ksize][AA][4][1]*lhs[j][i][ksize-1][CC][2][4];
      lhs[j][i][ksize][BB][2][2] = lhs[j][i][ksize][BB][2][2] - lhs[j][i][ksize][AA][0][2]*lhs[j][i][ksize-1][CC][2][0]
                              - lhs[j][i][ksize][AA][1][2]*lhs[j][i][ksize-1][CC][2][1]
                              - lhs[j][i][ksize][AA][2][2]*lhs[j][i][ksize-1][CC][2][2]
                              - lhs[j][i][ksize][AA][3][2]*lhs[j][i][ksize-1][CC][2][3]
                              - lhs[j][i][ksize][AA][4][2]*lhs[j][i][ksize-1][CC][2][4];
      lhs[j][i][ksize][BB][2][3] = lhs[j][i][ksize][BB][2][3] - lhs[j][i][ksize][AA][0][3]*lhs[j][i][ksize-1][CC][2][0]
                              - lhs[j][i][ksize][AA][1][3]*lhs[j][i][ksize-1][CC][2][1]
                              - lhs[j][i][ksize][AA][2][3]*lhs[j][i][ksize-1][CC][2][2]
                              - lhs[j][i][ksize][AA][3][3]*lhs[j][i][ksize-1][CC][2][3]
                              - lhs[j][i][ksize][AA][4][3]*lhs[j][i][ksize-1][CC][2][4];
      lhs[j][i][ksize][BB][2][4] = lhs[j][i][ksize][BB][2][4] - lhs[j][i][ksize][AA][0][4]*lhs[j][i][ksize-1][CC][2][0]
                              - lhs[j][i][ksize][AA][1][4]*lhs[j][i][ksize-1][CC][2][1]
                              - lhs[j][i][ksize][AA][2][4]*lhs[j][i][ksize-1][CC][2][2]
                              - lhs[j][i][ksize][AA][3][4]*lhs[j][i][ksize-1][CC][2][3]
                              - lhs[j][i][ksize][AA][4][4]*lhs[j][i][ksize-1][CC][2][4];
      lhs[j][i][ksize][BB][3][0] = lhs[j][i][ksize][BB][3][0] - lhs[j][i][ksize][AA][0][0]*lhs[j][i][ksize-1][CC][3][0]
                              - lhs[j][i][ksize][AA][1][0]*lhs[j][i][ksize-1][CC][3][1]
                              - lhs[j][i][ksize][AA][2][0]*lhs[j][i][ksize-1][CC][3][2]
                              - lhs[j][i][ksize][AA][3][0]*lhs[j][i][ksize-1][CC][3][3]
                              - lhs[j][i][ksize][AA][4][0]*lhs[j][i][ksize-1][CC][3][4];
      lhs[j][i][ksize][BB][3][1] = lhs[j][i][ksize][BB][3][1] - lhs[j][i][ksize][AA][0][1]*lhs[j][i][ksize-1][CC][3][0]
                              - lhs[j][i][ksize][AA][1][1]*lhs[j][i][ksize-1][CC][3][1]
                              - lhs[j][i][ksize][AA][2][1]*lhs[j][i][ksize-1][CC][3][2]
                              - lhs[j][i][ksize][AA][3][1]*lhs[j][i][ksize-1][CC][3][3]
                              - lhs[j][i][ksize][AA][4][1]*lhs[j][i][ksize-1][CC][3][4];
      lhs[j][i][ksize][BB][3][2] = lhs[j][i][ksize][BB][3][2] - lhs[j][i][ksize][AA][0][2]*lhs[j][i][ksize-1][CC][3][0]
                              - lhs[j][i][ksize][AA][1][2]*lhs[j][i][ksize-1][CC][3][1]
                              - lhs[j][i][ksize][AA][2][2]*lhs[j][i][ksize-1][CC][3][2]
                              - lhs[j][i][ksize][AA][3][2]*lhs[j][i][ksize-1][CC][3][3]
                              - lhs[j][i][ksize][AA][4][2]*lhs[j][i][ksize-1][CC][3][4];
      lhs[j][i][ksize][BB][3][3] = lhs[j][i][ksize][BB][3][3] - lhs[j][i][ksize][AA][0][3]*lhs[j][i][ksize-1][CC][3][0]
                              - lhs[j][i][ksize][AA][1][3]*lhs[j][i][ksize-1][CC][3][1]
                              - lhs[j][i][ksize][AA][2][3]*lhs[j][i][ksize-1][CC][3][2]
                              - lhs[j][i][ksize][AA][3][3]*lhs[j][i][ksize-1][CC][3][3]
                              - lhs[j][i][ksize][AA][4][3]*lhs[j][i][ksize-1][CC][3][4];
      lhs[j][i][ksize][BB][3][4] = lhs[j][i][ksize][BB][3][4] - lhs[j][i][ksize][AA][0][4]*lhs[j][i][ksize-1][CC][3][0]
                              - lhs[j][i][ksize][AA][1][4]*lhs[j][i][ksize-1][CC][3][1]
                              - lhs[j][i][ksize][AA][2][4]*lhs[j][i][ksize-1][CC][3][2]
                              - lhs[j][i][ksize][AA][3][4]*lhs[j][i][ksize-1][CC][3][3]
                              - lhs[j][i][ksize][AA][4][4]*lhs[j][i][ksize-1][CC][3][4];
      lhs[j][i][ksize][BB][4][0] = lhs[j][i][ksize][BB][4][0] - lhs[j][i][ksize][AA][0][0]*lhs[j][i][ksize-1][CC][4][0]
                              - lhs[j][i][ksize][AA][1][0]*lhs[j][i][ksize-1][CC][4][1]
                              - lhs[j][i][ksize][AA][2][0]*lhs[j][i][ksize-1][CC][4][2]
                              - lhs[j][i][ksize][AA][3][0]*lhs[j][i][ksize-1][CC][4][3]
                              - lhs[j][i][ksize][AA][4][0]*lhs[j][i][ksize-1][CC][4][4];
      lhs[j][i][ksize][BB][4][1] = lhs[j][i][ksize][BB][4][1] - lhs[j][i][ksize][AA][0][1]*lhs[j][i][ksize-1][CC][4][0]
                              - lhs[j][i][ksize][AA][1][1]*lhs[j][i][ksize-1][CC][4][1]
                              - lhs[j][i][ksize][AA][2][1]*lhs[j][i][ksize-1][CC][4][2]
                              - lhs[j][i][ksize][AA][3][1]*lhs[j][i][ksize-1][CC][4][3]
                              - lhs[j][i][ksize][AA][4][1]*lhs[j][i][ksize-1][CC][4][4];
      lhs[j][i][ksize][BB][4][2] = lhs[j][i][ksize][BB][4][2] - lhs[j][i][ksize][AA][0][2]*lhs[j][i][ksize-1][CC][4][0]
                              - lhs[j][i][ksize][AA][1][2]*lhs[j][i][ksize-1][CC][4][1]
                              - lhs[j][i][ksize][AA][2][2]*lhs[j][i][ksize-1][CC][4][2]
                              - lhs[j][i][ksize][AA][3][2]*lhs[j][i][ksize-1][CC][4][3]
                              - lhs[j][i][ksize][AA][4][2]*lhs[j][i][ksize-1][CC][4][4];
      lhs[j][i][ksize][BB][4][3] = lhs[j][i][ksize][BB][4][3] - lhs[j][i][ksize][AA][0][3]*lhs[j][i][ksize-1][CC][4][0]
                              - lhs[j][i][ksize][AA][1][3]*lhs[j][i][ksize-1][CC][4][1]
                              - lhs[j][i][ksize][AA][2][3]*lhs[j][i][ksize-1][CC][4][2]
                              - lhs[j][i][ksize][AA][3][3]*lhs[j][i][ksize-1][CC][4][3]
                              - lhs[j][i][ksize][AA][4][3]*lhs[j][i][ksize-1][CC][4][4];
      lhs[j][i][ksize][BB][4][4] = lhs[j][i][ksize][BB][4][4] - lhs[j][i][ksize][AA][0][4]*lhs[j][i][ksize-1][CC][4][0]
                              - lhs[j][i][ksize][AA][1][4]*lhs[j][i][ksize-1][CC][4][1]
                              - lhs[j][i][ksize][AA][2][4]*lhs[j][i][ksize-1][CC][4][2]
                              - lhs[j][i][ksize][AA][3][4]*lhs[j][i][ksize-1][CC][4][3]
                              - lhs[j][i][ksize][AA][4][4]*lhs[j][i][ksize-1][CC][4][4];

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
              - lhs[j][i][k][CC][n][m]*rhs[k+1][j][i][n];
          }
        }
      }
    }
  }
  if (timeron) timer_stop(t_zsolve);
}
