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
/*void x_lhsinit(int k, int j, int ni, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5])//(double lhs[][3][5][5], int ni)
{
  int i, m, n;

  //---------------------------------------------------------------------
  // zero the whole left hand side for starters
  // set all diagonal values to 1. This is overkill, but convenient
  //---------------------------------------------------------------------
  i = 0;
  for (n = 0; n < 5; n++) {
    for (m = 0; m < 5; m++) {
      lhs[k][j][i][0][n][m] = 0.0;
      lhs[k][j][i][1][n][m] = 0.0;
      lhs[k][j][i][2][n][m] = 0.0;
    }
    lhs[k][j][i][1][n][n] = 1.0;
  }
  i = ni;
  for (n = 0; n < 5; n++) {
    for (m = 0; m < 5; m++) {
      lhs[k][j][i][0][n][m] = 0.0;
      lhs[k][j][i][1][n][m] = 0.0;
      lhs[k][j][i][2][n][m] = 0.0;
    }
    lhs[k][j][i][1][n][n] = 1.0;
  }
}

//#pragma acc routine
//x_matvec_sub(k, j, i, AA, k, j, i-1, k, j, i, lhs, rhs);//x_matvec_sub(lhs[k][j][i][AA], rhs[k][j][i-1], rhs[k][j][i]);
//(double ablock[5][5], double avec[5], double bvec[5])//x_matvec_sub(k, j, i, AA, i-1);
//x_matvec_sub(lhs[k][j][i][first], rhs[k][j][second], rhs[k][j][i]);     lhs[k1][j1][i1][first]    rhs[k2][j2][i2]   rhs[k3][j3][i3]
void x_matvec_sub(int k1, int j1, int i1, int first, int k2, int j2, int i2, int k3, int j3, int i3, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5], double rhs[KMAX][JMAXP+1][IMAXP+1][5])//(double ablock[5][5], double avec[5], double bvec[5])//x_matvec_sub(k, j, i, AA, i-1);//x_matvec_sub(lhs[k][j][i][first], rhs[k][j][second], rhs[k][j][i]);
{
  //---------------------------------------------------------------------
  // rhs[kc][jc][ic][i] = rhs[kc][jc][ic][i] 
  // $                  - lhs[ia][ablock][0][i]*
  //---------------------------------------------------------------------
  rhs[k3][j3][i3][0] = rhs[k3][j3][i3][0] - lhs[k1][j1][i1][first][0][0]*rhs[k2][j2][i2][0]
                    - lhs[k1][j1][i1][first][1][0]*rhs[k2][j2][i2][1]
                    - lhs[k1][j1][i1][first][2][0]*rhs[k2][j2][i2][2]
                    - lhs[k1][j1][i1][first][3][0]*rhs[k2][j2][i2][3]
                    - lhs[k1][j1][i1][first][4][0]*rhs[k2][j2][i2][4];
  rhs[k3][j3][i3][1] = rhs[k3][j3][i3][1] - lhs[k1][j1][i1][first][0][1]*rhs[k2][j2][i2][0]
                    - lhs[k1][j1][i1][first][1][1]*rhs[k2][j2][i2][1]
                    - lhs[k1][j1][i1][first][2][1]*rhs[k2][j2][i2][2]
                    - lhs[k1][j1][i1][first][3][1]*rhs[k2][j2][i2][3]
                    - lhs[k1][j1][i1][first][4][1]*rhs[k2][j2][i2][4];
  rhs[k3][j3][i3][2] = rhs[k3][j3][i3][2] - lhs[k1][j1][i1][first][0][2]*rhs[k2][j2][i2][0]
                    - lhs[k1][j1][i1][first][1][2]*rhs[k2][j2][i2][1]
                    - lhs[k1][j1][i1][first][2][2]*rhs[k2][j2][i2][2]
                    - lhs[k1][j1][i1][first][3][2]*rhs[k2][j2][i2][3]
                    - lhs[k1][j1][i1][first][4][2]*rhs[k2][j2][i2][4];
  rhs[k3][j3][i3][3] = rhs[k3][j3][i3][3] - lhs[k1][j1][i1][first][0][3]*rhs[k2][j2][i2][0]
                    - lhs[k1][j1][i1][first][1][3]*rhs[k2][j2][i2][1]
                    - lhs[k1][j1][i1][first][2][3]*rhs[k2][j2][i2][2]
                    - lhs[k1][j1][i1][first][3][3]*rhs[k2][j2][i2][3]
                    - lhs[k1][j1][i1][first][4][3]*rhs[k2][j2][i2][4];
  rhs[k3][j3][i3][4] = rhs[k3][j3][i3][4] - lhs[k1][j1][i1][first][0][4]*rhs[k2][j2][i2][0]
                    - lhs[k1][j1][i1][first][1][4]*rhs[k2][j2][i2][1]
                    - lhs[k1][j1][i1][first][2][4]*rhs[k2][j2][i2][2]
                    - lhs[k1][j1][i1][first][3][4]*rhs[k2][j2][i2][3]
                    - lhs[k1][j1][i1][first][4][4]*rhs[k2][j2][i2][4];
}


//---------------------------------------------------------------------
// subtracts a(i,j,k) X b(i,j,k) from c(i,j,k)
//---------------------------------------------------------------------
//#pragma acc routine
//x_matmul_sub(k, j, i, AA, k, j, i-1, CC, k, j, i, BB, lhs);
//(double ablock[5][5], double bblock[5][5], double cblock[5][5])
//x_matmul_sub(lhs[k][j][i][first], lhs[k][j][second][third], lhs[k][j][i][forth]);    lhs[k1][j1][i1][first]  lhs[k2][j2][i2][second] lhs[k3][j3][i3][third]
void x_matmul_sub(int k1, int j1, int i1, int first, int k2, int j2, int i2, int second, int k3, int j3, int i3, int third, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5])//(double ablock[5][5], double bblock[5][5], double cblock[5][5])
{
  lhs[k3][j3][i3][third][0][0] = lhs[k3][j3][i3][third][0][0] - lhs[k1][j1][i1][first][0][0]*lhs[k2][j2][i2][second][0][0]
                              - lhs[k1][j1][i1][first][1][0]*lhs[k2][j2][i2][second][0][1]
                              - lhs[k1][j1][i1][first][2][0]*lhs[k2][j2][i2][second][0][2]
                              - lhs[k1][j1][i1][first][3][0]*lhs[k2][j2][i2][second][0][3]
                              - lhs[k1][j1][i1][first][4][0]*lhs[k2][j2][i2][second][0][4];
  lhs[k3][j3][i3][third][0][1] = lhs[k3][j3][i3][third][0][1] - lhs[k1][j1][i1][first][0][1]*lhs[k2][j2][i2][second][0][0]
                              - lhs[k1][j1][i1][first][1][1]*lhs[k2][j2][i2][second][0][1]
                              - lhs[k1][j1][i1][first][2][1]*lhs[k2][j2][i2][second][0][2]
                              - lhs[k1][j1][i1][first][3][1]*lhs[k2][j2][i2][second][0][3]
                              - lhs[k1][j1][i1][first][4][1]*lhs[k2][j2][i2][second][0][4];
  lhs[k3][j3][i3][third][0][2] = lhs[k3][j3][i3][third][0][2] - lhs[k1][j1][i1][first][0][2]*lhs[k2][j2][i2][second][0][0]
                              - lhs[k1][j1][i1][first][1][2]*lhs[k2][j2][i2][second][0][1]
                              - lhs[k1][j1][i1][first][2][2]*lhs[k2][j2][i2][second][0][2]
                              - lhs[k1][j1][i1][first][3][2]*lhs[k2][j2][i2][second][0][3]
                              - lhs[k1][j1][i1][first][4][2]*lhs[k2][j2][i2][second][0][4];
  lhs[k3][j3][i3][third][0][3] = lhs[k3][j3][i3][third][0][3] - lhs[k1][j1][i1][first][0][3]*lhs[k2][j2][i2][second][0][0]
                              - lhs[k1][j1][i1][first][1][3]*lhs[k2][j2][i2][second][0][1]
                              - lhs[k1][j1][i1][first][2][3]*lhs[k2][j2][i2][second][0][2]
                              - lhs[k1][j1][i1][first][3][3]*lhs[k2][j2][i2][second][0][3]
                              - lhs[k1][j1][i1][first][4][3]*lhs[k2][j2][i2][second][0][4];
  lhs[k3][j3][i3][third][0][4] = lhs[k3][j3][i3][third][0][4] - lhs[k1][j1][i1][first][0][4]*lhs[k2][j2][i2][second][0][0]
                              - lhs[k1][j1][i1][first][1][4]*lhs[k2][j2][i2][second][0][1]
                              - lhs[k1][j1][i1][first][2][4]*lhs[k2][j2][i2][second][0][2]
                              - lhs[k1][j1][i1][first][3][4]*lhs[k2][j2][i2][second][0][3]
                              - lhs[k1][j1][i1][first][4][4]*lhs[k2][j2][i2][second][0][4];
  lhs[k3][j3][i3][third][1][0] = lhs[k3][j3][i3][third][1][0] - lhs[k1][j1][i1][first][0][0]*lhs[k2][j2][i2][second][1][0]
                              - lhs[k1][j1][i1][first][1][0]*lhs[k2][j2][i2][second][1][1]
                              - lhs[k1][j1][i1][first][2][0]*lhs[k2][j2][i2][second][1][2]
                              - lhs[k1][j1][i1][first][3][0]*lhs[k2][j2][i2][second][1][3]
                              - lhs[k1][j1][i1][first][4][0]*lhs[k2][j2][i2][second][1][4];
  lhs[k3][j3][i3][third][1][1] = lhs[k3][j3][i3][third][1][1] - lhs[k1][j1][i1][first][0][1]*lhs[k2][j2][i2][second][1][0]
                              - lhs[k1][j1][i1][first][1][1]*lhs[k2][j2][i2][second][1][1]
                              - lhs[k1][j1][i1][first][2][1]*lhs[k2][j2][i2][second][1][2]
                              - lhs[k1][j1][i1][first][3][1]*lhs[k2][j2][i2][second][1][3]
                              - lhs[k1][j1][i1][first][4][1]*lhs[k2][j2][i2][second][1][4];
  lhs[k3][j3][i3][third][1][2] = lhs[k3][j3][i3][third][1][2] - lhs[k1][j1][i1][first][0][2]*lhs[k2][j2][i2][second][1][0]
                              - lhs[k1][j1][i1][first][1][2]*lhs[k2][j2][i2][second][1][1]
                              - lhs[k1][j1][i1][first][2][2]*lhs[k2][j2][i2][second][1][2]
                              - lhs[k1][j1][i1][first][3][2]*lhs[k2][j2][i2][second][1][3]
                              - lhs[k1][j1][i1][first][4][2]*lhs[k2][j2][i2][second][1][4];
  lhs[k3][j3][i3][third][1][3] = lhs[k3][j3][i3][third][1][3] - lhs[k1][j1][i1][first][0][3]*lhs[k2][j2][i2][second][1][0]
                              - lhs[k1][j1][i1][first][1][3]*lhs[k2][j2][i2][second][1][1]
                              - lhs[k1][j1][i1][first][2][3]*lhs[k2][j2][i2][second][1][2]
                              - lhs[k1][j1][i1][first][3][3]*lhs[k2][j2][i2][second][1][3]
                              - lhs[k1][j1][i1][first][4][3]*lhs[k2][j2][i2][second][1][4];
  lhs[k3][j3][i3][third][1][4] = lhs[k3][j3][i3][third][1][4] - lhs[k1][j1][i1][first][0][4]*lhs[k2][j2][i2][second][1][0]
                              - lhs[k1][j1][i1][first][1][4]*lhs[k2][j2][i2][second][1][1]
                              - lhs[k1][j1][i1][first][2][4]*lhs[k2][j2][i2][second][1][2]
                              - lhs[k1][j1][i1][first][3][4]*lhs[k2][j2][i2][second][1][3]
                              - lhs[k1][j1][i1][first][4][4]*lhs[k2][j2][i2][second][1][4];
  lhs[k3][j3][i3][third][2][0] = lhs[k3][j3][i3][third][2][0] - lhs[k1][j1][i1][first][0][0]*lhs[k2][j2][i2][second][2][0]
                              - lhs[k1][j1][i1][first][1][0]*lhs[k2][j2][i2][second][2][1]
                              - lhs[k1][j1][i1][first][2][0]*lhs[k2][j2][i2][second][2][2]
                              - lhs[k1][j1][i1][first][3][0]*lhs[k2][j2][i2][second][2][3]
                              - lhs[k1][j1][i1][first][4][0]*lhs[k2][j2][i2][second][2][4];
  lhs[k3][j3][i3][third][2][1] = lhs[k3][j3][i3][third][2][1] - lhs[k1][j1][i1][first][0][1]*lhs[k2][j2][i2][second][2][0]
                              - lhs[k1][j1][i1][first][1][1]*lhs[k2][j2][i2][second][2][1]
                              - lhs[k1][j1][i1][first][2][1]*lhs[k2][j2][i2][second][2][2]
                              - lhs[k1][j1][i1][first][3][1]*lhs[k2][j2][i2][second][2][3]
                              - lhs[k1][j1][i1][first][4][1]*lhs[k2][j2][i2][second][2][4];
  lhs[k3][j3][i3][third][2][2] = lhs[k3][j3][i3][third][2][2] - lhs[k1][j1][i1][first][0][2]*lhs[k2][j2][i2][second][2][0]
                              - lhs[k1][j1][i1][first][1][2]*lhs[k2][j2][i2][second][2][1]
                              - lhs[k1][j1][i1][first][2][2]*lhs[k2][j2][i2][second][2][2]
                              - lhs[k1][j1][i1][first][3][2]*lhs[k2][j2][i2][second][2][3]
                              - lhs[k1][j1][i1][first][4][2]*lhs[k2][j2][i2][second][2][4];
  lhs[k3][j3][i3][third][2][3] = lhs[k3][j3][i3][third][2][3] - lhs[k1][j1][i1][first][0][3]*lhs[k2][j2][i2][second][2][0]
                              - lhs[k1][j1][i1][first][1][3]*lhs[k2][j2][i2][second][2][1]
                              - lhs[k1][j1][i1][first][2][3]*lhs[k2][j2][i2][second][2][2]
                              - lhs[k1][j1][i1][first][3][3]*lhs[k2][j2][i2][second][2][3]
                              - lhs[k1][j1][i1][first][4][3]*lhs[k2][j2][i2][second][2][4];
  lhs[k3][j3][i3][third][2][4] = lhs[k3][j3][i3][third][2][4] - lhs[k1][j1][i1][first][0][4]*lhs[k2][j2][i2][second][2][0]
                              - lhs[k1][j1][i1][first][1][4]*lhs[k2][j2][i2][second][2][1]
                              - lhs[k1][j1][i1][first][2][4]*lhs[k2][j2][i2][second][2][2]
                              - lhs[k1][j1][i1][first][3][4]*lhs[k2][j2][i2][second][2][3]
                              - lhs[k1][j1][i1][first][4][4]*lhs[k2][j2][i2][second][2][4];
  lhs[k3][j3][i3][third][3][0] = lhs[k3][j3][i3][third][3][0] - lhs[k1][j1][i1][first][0][0]*lhs[k2][j2][i2][second][3][0]
                              - lhs[k1][j1][i1][first][1][0]*lhs[k2][j2][i2][second][3][1]
                              - lhs[k1][j1][i1][first][2][0]*lhs[k2][j2][i2][second][3][2]
                              - lhs[k1][j1][i1][first][3][0]*lhs[k2][j2][i2][second][3][3]
                              - lhs[k1][j1][i1][first][4][0]*lhs[k2][j2][i2][second][3][4];
  lhs[k3][j3][i3][third][3][1] = lhs[k3][j3][i3][third][3][1] - lhs[k1][j1][i1][first][0][1]*lhs[k2][j2][i2][second][3][0]
                              - lhs[k1][j1][i1][first][1][1]*lhs[k2][j2][i2][second][3][1]
                              - lhs[k1][j1][i1][first][2][1]*lhs[k2][j2][i2][second][3][2]
                              - lhs[k1][j1][i1][first][3][1]*lhs[k2][j2][i2][second][3][3]
                              - lhs[k1][j1][i1][first][4][1]*lhs[k2][j2][i2][second][3][4];
  lhs[k3][j3][i3][third][3][2] = lhs[k3][j3][i3][third][3][2] - lhs[k1][j1][i1][first][0][2]*lhs[k2][j2][i2][second][3][0]
                              - lhs[k1][j1][i1][first][1][2]*lhs[k2][j2][i2][second][3][1]
                              - lhs[k1][j1][i1][first][2][2]*lhs[k2][j2][i2][second][3][2]
                              - lhs[k1][j1][i1][first][3][2]*lhs[k2][j2][i2][second][3][3]
                              - lhs[k1][j1][i1][first][4][2]*lhs[k2][j2][i2][second][3][4];
  lhs[k3][j3][i3][third][3][3] = lhs[k3][j3][i3][third][3][3] - lhs[k1][j1][i1][first][0][3]*lhs[k2][j2][i2][second][3][0]
                              - lhs[k1][j1][i1][first][1][3]*lhs[k2][j2][i2][second][3][1]
                              - lhs[k1][j1][i1][first][2][3]*lhs[k2][j2][i2][second][3][2]
                              - lhs[k1][j1][i1][first][3][3]*lhs[k2][j2][i2][second][3][3]
                              - lhs[k1][j1][i1][first][4][3]*lhs[k2][j2][i2][second][3][4];
  lhs[k3][j3][i3][third][3][4] = lhs[k3][j3][i3][third][3][4] - lhs[k1][j1][i1][first][0][4]*lhs[k2][j2][i2][second][3][0]
                              - lhs[k1][j1][i1][first][1][4]*lhs[k2][j2][i2][second][3][1]
                              - lhs[k1][j1][i1][first][2][4]*lhs[k2][j2][i2][second][3][2]
                              - lhs[k1][j1][i1][first][3][4]*lhs[k2][j2][i2][second][3][3]
                              - lhs[k1][j1][i1][first][4][4]*lhs[k2][j2][i2][second][3][4];
  lhs[k3][j3][i3][third][4][0] = lhs[k3][j3][i3][third][4][0] - lhs[k1][j1][i1][first][0][0]*lhs[k2][j2][i2][second][4][0]
                              - lhs[k1][j1][i1][first][1][0]*lhs[k2][j2][i2][second][4][1]
                              - lhs[k1][j1][i1][first][2][0]*lhs[k2][j2][i2][second][4][2]
                              - lhs[k1][j1][i1][first][3][0]*lhs[k2][j2][i2][second][4][3]
                              - lhs[k1][j1][i1][first][4][0]*lhs[k2][j2][i2][second][4][4];
  lhs[k3][j3][i3][third][4][1] = lhs[k3][j3][i3][third][4][1] - lhs[k1][j1][i1][first][0][1]*lhs[k2][j2][i2][second][4][0]
                              - lhs[k1][j1][i1][first][1][1]*lhs[k2][j2][i2][second][4][1]
                              - lhs[k1][j1][i1][first][2][1]*lhs[k2][j2][i2][second][4][2]
                              - lhs[k1][j1][i1][first][3][1]*lhs[k2][j2][i2][second][4][3]
                              - lhs[k1][j1][i1][first][4][1]*lhs[k2][j2][i2][second][4][4];
  lhs[k3][j3][i3][third][4][2] = lhs[k3][j3][i3][third][4][2] - lhs[k1][j1][i1][first][0][2]*lhs[k2][j2][i2][second][4][0]
                              - lhs[k1][j1][i1][first][1][2]*lhs[k2][j2][i2][second][4][1]
                              - lhs[k1][j1][i1][first][2][2]*lhs[k2][j2][i2][second][4][2]
                              - lhs[k1][j1][i1][first][3][2]*lhs[k2][j2][i2][second][4][3]
                              - lhs[k1][j1][i1][first][4][2]*lhs[k2][j2][i2][second][4][4];
  lhs[k3][j3][i3][third][4][3] = lhs[k3][j3][i3][third][4][3] - lhs[k1][j1][i1][first][0][3]*lhs[k2][j2][i2][second][4][0]
                              - lhs[k1][j1][i1][first][1][3]*lhs[k2][j2][i2][second][4][1]
                              - lhs[k1][j1][i1][first][2][3]*lhs[k2][j2][i2][second][4][2]
                              - lhs[k1][j1][i1][first][3][3]*lhs[k2][j2][i2][second][4][3]
                              - lhs[k1][j1][i1][first][4][3]*lhs[k2][j2][i2][second][4][4];
  lhs[k3][j3][i3][third][4][4] = lhs[k3][j3][i3][third][4][4] - lhs[k1][j1][i1][first][0][4]*lhs[k2][j2][i2][second][4][0]
                              - lhs[k1][j1][i1][first][1][4]*lhs[k2][j2][i2][second][4][1]
                              - lhs[k1][j1][i1][first][2][4]*lhs[k2][j2][i2][second][4][2]
                              - lhs[k1][j1][i1][first][3][4]*lhs[k2][j2][i2][second][4][3]
                              - lhs[k1][j1][i1][first][4][4]*lhs[k2][j2][i2][second][4][4];
}


//#pragma acc routine
//x_binvcrhs( k, j, 0, BB, k, j, 0, CC, k, j, 0, lhs, rhs );
//x_binvcrhs( lhs[k][j][i][first], lhs[k][j][i][second], rhs[k][j][i] ); lhs[k1][j1][i1][first] lhs[k2][j2][i2][second] rhs[k3][j3][i3]
void x_binvcrhs(int k1, int j1, int i1, int first, int k2, int j2, int i2, int second, int k3, int j3, int i3, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5], double rhs[KMAX][JMAXP+1][IMAXP+1][5])//(double lhs[5][5], double c[5][5], double r[5]) //x_binvcrhs( lhs[k][j][i][BB], lhs[k][j][i][CC][, rhs[k][j][i][ );
{
  double pivot, coeff;

  pivot = 1.00/lhs[k1][j1][i1][first][0][0];
  lhs[k1][j1][i1][first][1][0] = lhs[k1][j1][i1][first][1][0]*pivot;
  lhs[k1][j1][i1][first][2][0] = lhs[k1][j1][i1][first][2][0]*pivot;
  lhs[k1][j1][i1][first][3][0] = lhs[k1][j1][i1][first][3][0]*pivot;
  lhs[k1][j1][i1][first][4][0] = lhs[k1][j1][i1][first][4][0]*pivot;
  lhs[k2][j2][i2][second][0][0] = lhs[k2][j2][i2][second][0][0]*pivot;
  lhs[k2][j2][i2][second][1][0] = lhs[k2][j2][i2][second][1][0]*pivot;
  lhs[k2][j2][i2][second][2][0] = lhs[k2][j2][i2][second][2][0]*pivot;
  lhs[k2][j2][i2][second][3][0] = lhs[k2][j2][i2][second][3][0]*pivot;
  lhs[k2][j2][i2][second][4][0] = lhs[k2][j2][i2][second][4][0]*pivot;
  rhs[k3][j3][i3][0]   = rhs[k3][j3][i3][0]  *pivot;

  coeff = lhs[k1][j1][i1][first][0][1];
  lhs[k1][j1][i1][first][1][1]= lhs[k1][j1][i1][first][1][1] - coeff*lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][1]= lhs[k1][j1][i1][first][2][1] - coeff*lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][1]= lhs[k1][j1][i1][first][3][1] - coeff*lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][1]= lhs[k1][j1][i1][first][4][1] - coeff*lhs[k1][j1][i1][first][4][0];
  lhs[k2][j2][i2][second][0][1] = lhs[k2][j2][i2][second][0][1] - coeff*lhs[k2][j2][i2][second][0][0];
  lhs[k2][j2][i2][second][1][1] = lhs[k2][j2][i2][second][1][1] - coeff*lhs[k2][j2][i2][second][1][0];
  lhs[k2][j2][i2][second][2][1] = lhs[k2][j2][i2][second][2][1] - coeff*lhs[k2][j2][i2][second][2][0];
  lhs[k2][j2][i2][second][3][1] = lhs[k2][j2][i2][second][3][1] - coeff*lhs[k2][j2][i2][second][3][0];
  lhs[k2][j2][i2][second][4][1] = lhs[k2][j2][i2][second][4][1] - coeff*lhs[k2][j2][i2][second][4][0];
  rhs[k3][j3][i3][1]   = rhs[k3][j3][i3][1]   - coeff*rhs[k3][j3][i3][0];

  coeff = lhs[k1][j1][i1][first][0][2];
  lhs[k1][j1][i1][first][1][2]= lhs[k1][j1][i1][first][1][2] - coeff*lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][2]= lhs[k1][j1][i1][first][2][2] - coeff*lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][2]= lhs[k1][j1][i1][first][3][2] - coeff*lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][2]= lhs[k1][j1][i1][first][4][2] - coeff*lhs[k1][j1][i1][first][4][0];
  lhs[k2][j2][i2][second][0][2] = lhs[k2][j2][i2][second][0][2] - coeff*lhs[k2][j2][i2][second][0][0];
  lhs[k2][j2][i2][second][1][2] = lhs[k2][j2][i2][second][1][2] - coeff*lhs[k2][j2][i2][second][1][0];
  lhs[k2][j2][i2][second][2][2] = lhs[k2][j2][i2][second][2][2] - coeff*lhs[k2][j2][i2][second][2][0];
  lhs[k2][j2][i2][second][3][2] = lhs[k2][j2][i2][second][3][2] - coeff*lhs[k2][j2][i2][second][3][0];
  lhs[k2][j2][i2][second][4][2] = lhs[k2][j2][i2][second][4][2] - coeff*lhs[k2][j2][i2][second][4][0];
  rhs[k3][j3][i3][2]   = rhs[k3][j3][i3][2]   - coeff*rhs[k3][j3][i3][0];

  coeff = lhs[k1][j1][i1][first][0][3];
  lhs[k1][j1][i1][first][1][3]= lhs[k1][j1][i1][first][1][3] - coeff*lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][3]= lhs[k1][j1][i1][first][2][3] - coeff*lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][3]= lhs[k1][j1][i1][first][3][3] - coeff*lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][3]= lhs[k1][j1][i1][first][4][3] - coeff*lhs[k1][j1][i1][first][4][0];
  lhs[k2][j2][i2][second][0][3] = lhs[k2][j2][i2][second][0][3] - coeff*lhs[k2][j2][i2][second][0][0];
  lhs[k2][j2][i2][second][1][3] = lhs[k2][j2][i2][second][1][3] - coeff*lhs[k2][j2][i2][second][1][0];
  lhs[k2][j2][i2][second][2][3] = lhs[k2][j2][i2][second][2][3] - coeff*lhs[k2][j2][i2][second][2][0];
  lhs[k2][j2][i2][second][3][3] = lhs[k2][j2][i2][second][3][3] - coeff*lhs[k2][j2][i2][second][3][0];
  lhs[k2][j2][i2][second][4][3] = lhs[k2][j2][i2][second][4][3] - coeff*lhs[k2][j2][i2][second][4][0];
  rhs[k3][j3][i3][3]   = rhs[k3][j3][i3][3]   - coeff*rhs[k3][j3][i3][0];

  coeff = lhs[k1][j1][i1][first][0][4];
  lhs[k1][j1][i1][first][1][4]= lhs[k1][j1][i1][first][1][4] - coeff*lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][4]= lhs[k1][j1][i1][first][2][4] - coeff*lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][4]= lhs[k1][j1][i1][first][3][4] - coeff*lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][4]= lhs[k1][j1][i1][first][4][4] - coeff*lhs[k1][j1][i1][first][4][0];
  lhs[k2][j2][i2][second][0][4] = lhs[k2][j2][i2][second][0][4] - coeff*lhs[k2][j2][i2][second][0][0];
  lhs[k2][j2][i2][second][1][4] = lhs[k2][j2][i2][second][1][4] - coeff*lhs[k2][j2][i2][second][1][0];
  lhs[k2][j2][i2][second][2][4] = lhs[k2][j2][i2][second][2][4] - coeff*lhs[k2][j2][i2][second][2][0];
  lhs[k2][j2][i2][second][3][4] = lhs[k2][j2][i2][second][3][4] - coeff*lhs[k2][j2][i2][second][3][0];
  lhs[k2][j2][i2][second][4][4] = lhs[k2][j2][i2][second][4][4] - coeff*lhs[k2][j2][i2][second][4][0];
  rhs[k3][j3][i3][4]   = rhs[k3][j3][i3][4]   - coeff*rhs[k3][j3][i3][0];


  pivot = 1.00/lhs[k1][j1][i1][first][1][1];
  lhs[k1][j1][i1][first][2][1] = lhs[k1][j1][i1][first][2][1]*pivot;
  lhs[k1][j1][i1][first][3][1] = lhs[k1][j1][i1][first][3][1]*pivot;
  lhs[k1][j1][i1][first][4][1] = lhs[k1][j1][i1][first][4][1]*pivot;
  lhs[k2][j2][i2][second][0][1] = lhs[k2][j2][i2][second][0][1]*pivot;
  lhs[k2][j2][i2][second][1][1] = lhs[k2][j2][i2][second][1][1]*pivot;
  lhs[k2][j2][i2][second][2][1] = lhs[k2][j2][i2][second][2][1]*pivot;
  lhs[k2][j2][i2][second][3][1] = lhs[k2][j2][i2][second][3][1]*pivot;
  lhs[k2][j2][i2][second][4][1] = lhs[k2][j2][i2][second][4][1]*pivot;
  rhs[k3][j3][i3][1]   = rhs[k3][j3][i3][1]  *pivot;

  coeff = lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][0]= lhs[k1][j1][i1][first][2][0] - coeff*lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][0]= lhs[k1][j1][i1][first][3][0] - coeff*lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][0]= lhs[k1][j1][i1][first][4][0] - coeff*lhs[k1][j1][i1][first][4][1];
  lhs[k2][j2][i2][second][0][0] = lhs[k2][j2][i2][second][0][0] - coeff*lhs[k2][j2][i2][second][0][1];
  lhs[k2][j2][i2][second][1][0] = lhs[k2][j2][i2][second][1][0] - coeff*lhs[k2][j2][i2][second][1][1];
  lhs[k2][j2][i2][second][2][0] = lhs[k2][j2][i2][second][2][0] - coeff*lhs[k2][j2][i2][second][2][1];
  lhs[k2][j2][i2][second][3][0] = lhs[k2][j2][i2][second][3][0] - coeff*lhs[k2][j2][i2][second][3][1];
  lhs[k2][j2][i2][second][4][0] = lhs[k2][j2][i2][second][4][0] - coeff*lhs[k2][j2][i2][second][4][1];
  rhs[k3][j3][i3][0]   = rhs[k3][j3][i3][0]   - coeff*rhs[k3][j3][i3][1];

  coeff = lhs[k1][j1][i1][first][1][2];
  lhs[k1][j1][i1][first][2][2]= lhs[k1][j1][i1][first][2][2] - coeff*lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][2]= lhs[k1][j1][i1][first][3][2] - coeff*lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][2]= lhs[k1][j1][i1][first][4][2] - coeff*lhs[k1][j1][i1][first][4][1];
  lhs[k2][j2][i2][second][0][2] = lhs[k2][j2][i2][second][0][2] - coeff*lhs[k2][j2][i2][second][0][1];
  lhs[k2][j2][i2][second][1][2] = lhs[k2][j2][i2][second][1][2] - coeff*lhs[k2][j2][i2][second][1][1];
  lhs[k2][j2][i2][second][2][2] = lhs[k2][j2][i2][second][2][2] - coeff*lhs[k2][j2][i2][second][2][1];
  lhs[k2][j2][i2][second][3][2] = lhs[k2][j2][i2][second][3][2] - coeff*lhs[k2][j2][i2][second][3][1];
  lhs[k2][j2][i2][second][4][2] = lhs[k2][j2][i2][second][4][2] - coeff*lhs[k2][j2][i2][second][4][1];
  rhs[k3][j3][i3][2]   = rhs[k3][j3][i3][2]   - coeff*rhs[k3][j3][i3][1];

  coeff = lhs[k1][j1][i1][first][1][3];
  lhs[k1][j1][i1][first][2][3]= lhs[k1][j1][i1][first][2][3] - coeff*lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][3]= lhs[k1][j1][i1][first][3][3] - coeff*lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][3]= lhs[k1][j1][i1][first][4][3] - coeff*lhs[k1][j1][i1][first][4][1];
  lhs[k2][j2][i2][second][0][3] = lhs[k2][j2][i2][second][0][3] - coeff*lhs[k2][j2][i2][second][0][1];
  lhs[k2][j2][i2][second][1][3] = lhs[k2][j2][i2][second][1][3] - coeff*lhs[k2][j2][i2][second][1][1];
  lhs[k2][j2][i2][second][2][3] = lhs[k2][j2][i2][second][2][3] - coeff*lhs[k2][j2][i2][second][2][1];
  lhs[k2][j2][i2][second][3][3] = lhs[k2][j2][i2][second][3][3] - coeff*lhs[k2][j2][i2][second][3][1];
  lhs[k2][j2][i2][second][4][3] = lhs[k2][j2][i2][second][4][3] - coeff*lhs[k2][j2][i2][second][4][1];
  rhs[k3][j3][i3][3]   = rhs[k3][j3][i3][3]   - coeff*rhs[k3][j3][i3][1];

  coeff = lhs[k1][j1][i1][first][1][4];
  lhs[k1][j1][i1][first][2][4]= lhs[k1][j1][i1][first][2][4] - coeff*lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][4]= lhs[k1][j1][i1][first][3][4] - coeff*lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][4]= lhs[k1][j1][i1][first][4][4] - coeff*lhs[k1][j1][i1][first][4][1];
  lhs[k2][j2][i2][second][0][4] = lhs[k2][j2][i2][second][0][4] - coeff*lhs[k2][j2][i2][second][0][1];
  lhs[k2][j2][i2][second][1][4] = lhs[k2][j2][i2][second][1][4] - coeff*lhs[k2][j2][i2][second][1][1];
  lhs[k2][j2][i2][second][2][4] = lhs[k2][j2][i2][second][2][4] - coeff*lhs[k2][j2][i2][second][2][1];
  lhs[k2][j2][i2][second][3][4] = lhs[k2][j2][i2][second][3][4] - coeff*lhs[k2][j2][i2][second][3][1];
  lhs[k2][j2][i2][second][4][4] = lhs[k2][j2][i2][second][4][4] - coeff*lhs[k2][j2][i2][second][4][1];
  rhs[k3][j3][i3][4]   = rhs[k3][j3][i3][4]   - coeff*rhs[k3][j3][i3][1];


  pivot = 1.00/lhs[k1][j1][i1][first][2][2];
  lhs[k1][j1][i1][first][3][2] = lhs[k1][j1][i1][first][3][2]*pivot;
  lhs[k1][j1][i1][first][4][2] = lhs[k1][j1][i1][first][4][2]*pivot;
  lhs[k2][j2][i2][second][0][2] = lhs[k2][j2][i2][second][0][2]*pivot;
  lhs[k2][j2][i2][second][1][2] = lhs[k2][j2][i2][second][1][2]*pivot;
  lhs[k2][j2][i2][second][2][2] = lhs[k2][j2][i2][second][2][2]*pivot;
  lhs[k2][j2][i2][second][3][2] = lhs[k2][j2][i2][second][3][2]*pivot;
  lhs[k2][j2][i2][second][4][2] = lhs[k2][j2][i2][second][4][2]*pivot;
  rhs[k3][j3][i3][2]   = rhs[k3][j3][i3][2]  *pivot;

  coeff = lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][0]= lhs[k1][j1][i1][first][3][0] - coeff*lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][0]= lhs[k1][j1][i1][first][4][0] - coeff*lhs[k1][j1][i1][first][4][2];
  lhs[k2][j2][i2][second][0][0] = lhs[k2][j2][i2][second][0][0] - coeff*lhs[k2][j2][i2][second][0][2];
  lhs[k2][j2][i2][second][1][0] = lhs[k2][j2][i2][second][1][0] - coeff*lhs[k2][j2][i2][second][1][2];
  lhs[k2][j2][i2][second][2][0] = lhs[k2][j2][i2][second][2][0] - coeff*lhs[k2][j2][i2][second][2][2];
  lhs[k2][j2][i2][second][3][0] = lhs[k2][j2][i2][second][3][0] - coeff*lhs[k2][j2][i2][second][3][2];
  lhs[k2][j2][i2][second][4][0] = lhs[k2][j2][i2][second][4][0] - coeff*lhs[k2][j2][i2][second][4][2];
  rhs[k3][j3][i3][0]   = rhs[k3][j3][i3][0]   - coeff*rhs[k3][j3][i3][2];

  coeff = lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][1]= lhs[k1][j1][i1][first][3][1] - coeff*lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][1]= lhs[k1][j1][i1][first][4][1] - coeff*lhs[k1][j1][i1][first][4][2];
  lhs[k2][j2][i2][second][0][1] = lhs[k2][j2][i2][second][0][1] - coeff*lhs[k2][j2][i2][second][0][2];
  lhs[k2][j2][i2][second][1][1] = lhs[k2][j2][i2][second][1][1] - coeff*lhs[k2][j2][i2][second][1][2];
  lhs[k2][j2][i2][second][2][1] = lhs[k2][j2][i2][second][2][1] - coeff*lhs[k2][j2][i2][second][2][2];
  lhs[k2][j2][i2][second][3][1] = lhs[k2][j2][i2][second][3][1] - coeff*lhs[k2][j2][i2][second][3][2];
  lhs[k2][j2][i2][second][4][1] = lhs[k2][j2][i2][second][4][1] - coeff*lhs[k2][j2][i2][second][4][2];
  rhs[k3][j3][i3][1]   = rhs[k3][j3][i3][1]   - coeff*rhs[k3][j3][i3][2];

  coeff = lhs[k1][j1][i1][first][2][3];
  lhs[k1][j1][i1][first][3][3]= lhs[k1][j1][i1][first][3][3] - coeff*lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][3]= lhs[k1][j1][i1][first][4][3] - coeff*lhs[k1][j1][i1][first][4][2];
  lhs[k2][j2][i2][second][0][3] = lhs[k2][j2][i2][second][0][3] - coeff*lhs[k2][j2][i2][second][0][2];
  lhs[k2][j2][i2][second][1][3] = lhs[k2][j2][i2][second][1][3] - coeff*lhs[k2][j2][i2][second][1][2];
  lhs[k2][j2][i2][second][2][3] = lhs[k2][j2][i2][second][2][3] - coeff*lhs[k2][j2][i2][second][2][2];
  lhs[k2][j2][i2][second][3][3] = lhs[k2][j2][i2][second][3][3] - coeff*lhs[k2][j2][i2][second][3][2];
  lhs[k2][j2][i2][second][4][3] = lhs[k2][j2][i2][second][4][3] - coeff*lhs[k2][j2][i2][second][4][2];
  rhs[k3][j3][i3][3]   = rhs[k3][j3][i3][3]   - coeff*rhs[k3][j3][i3][2];

  coeff = lhs[k1][j1][i1][first][2][4];
  lhs[k1][j1][i1][first][3][4]= lhs[k1][j1][i1][first][3][4] - coeff*lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][4]= lhs[k1][j1][i1][first][4][4] - coeff*lhs[k1][j1][i1][first][4][2];
  lhs[k2][j2][i2][second][0][4] = lhs[k2][j2][i2][second][0][4] - coeff*lhs[k2][j2][i2][second][0][2];
  lhs[k2][j2][i2][second][1][4] = lhs[k2][j2][i2][second][1][4] - coeff*lhs[k2][j2][i2][second][1][2];
  lhs[k2][j2][i2][second][2][4] = lhs[k2][j2][i2][second][2][4] - coeff*lhs[k2][j2][i2][second][2][2];
  lhs[k2][j2][i2][second][3][4] = lhs[k2][j2][i2][second][3][4] - coeff*lhs[k2][j2][i2][second][3][2];
  lhs[k2][j2][i2][second][4][4] = lhs[k2][j2][i2][second][4][4] - coeff*lhs[k2][j2][i2][second][4][2];
  rhs[k3][j3][i3][4]   = rhs[k3][j3][i3][4]   - coeff*rhs[k3][j3][i3][2];


  pivot = 1.00/lhs[k1][j1][i1][first][3][3];
  lhs[k1][j1][i1][first][4][3] = lhs[k1][j1][i1][first][4][3]*pivot;
  lhs[k2][j2][i2][second][0][3] = lhs[k2][j2][i2][second][0][3]*pivot;
  lhs[k2][j2][i2][second][1][3] = lhs[k2][j2][i2][second][1][3]*pivot;
  lhs[k2][j2][i2][second][2][3] = lhs[k2][j2][i2][second][2][3]*pivot;
  lhs[k2][j2][i2][second][3][3] = lhs[k2][j2][i2][second][3][3]*pivot;
  lhs[k2][j2][i2][second][4][3] = lhs[k2][j2][i2][second][4][3]*pivot;
  rhs[k3][j3][i3][3]   = rhs[k3][j3][i3][3]  *pivot;

  coeff = lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][0]= lhs[k1][j1][i1][first][4][0] - coeff*lhs[k1][j1][i1][first][4][3];
  lhs[k2][j2][i2][second][0][0] = lhs[k2][j2][i2][second][0][0] - coeff*lhs[k2][j2][i2][second][0][3];
  lhs[k2][j2][i2][second][1][0] = lhs[k2][j2][i2][second][1][0] - coeff*lhs[k2][j2][i2][second][1][3];
  lhs[k2][j2][i2][second][2][0] = lhs[k2][j2][i2][second][2][0] - coeff*lhs[k2][j2][i2][second][2][3];
  lhs[k2][j2][i2][second][3][0] = lhs[k2][j2][i2][second][3][0] - coeff*lhs[k2][j2][i2][second][3][3];
  lhs[k2][j2][i2][second][4][0] = lhs[k2][j2][i2][second][4][0] - coeff*lhs[k2][j2][i2][second][4][3];
  rhs[k3][j3][i3][0]   = rhs[k3][j3][i3][0]   - coeff*rhs[k3][j3][i3][3];

  coeff = lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][1]= lhs[k1][j1][i1][first][4][1] - coeff*lhs[k1][j1][i1][first][4][3];
  lhs[k2][j2][i2][second][0][1] = lhs[k2][j2][i2][second][0][1] - coeff*lhs[k2][j2][i2][second][0][3];
  lhs[k2][j2][i2][second][1][1] = lhs[k2][j2][i2][second][1][1] - coeff*lhs[k2][j2][i2][second][1][3];
  lhs[k2][j2][i2][second][2][1] = lhs[k2][j2][i2][second][2][1] - coeff*lhs[k2][j2][i2][second][2][3];
  lhs[k2][j2][i2][second][3][1] = lhs[k2][j2][i2][second][3][1] - coeff*lhs[k2][j2][i2][second][3][3];
  lhs[k2][j2][i2][second][4][1] = lhs[k2][j2][i2][second][4][1] - coeff*lhs[k2][j2][i2][second][4][3];
  rhs[k3][j3][i3][1]   = rhs[k3][j3][i3][1]   - coeff*rhs[k3][j3][i3][3];

  coeff = lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][2]= lhs[k1][j1][i1][first][4][2] - coeff*lhs[k1][j1][i1][first][4][3];
  lhs[k2][j2][i2][second][0][2] = lhs[k2][j2][i2][second][0][2] - coeff*lhs[k2][j2][i2][second][0][3];
  lhs[k2][j2][i2][second][1][2] = lhs[k2][j2][i2][second][1][2] - coeff*lhs[k2][j2][i2][second][1][3];
  lhs[k2][j2][i2][second][2][2] = lhs[k2][j2][i2][second][2][2] - coeff*lhs[k2][j2][i2][second][2][3];
  lhs[k2][j2][i2][second][3][2] = lhs[k2][j2][i2][second][3][2] - coeff*lhs[k2][j2][i2][second][3][3];
  lhs[k2][j2][i2][second][4][2] = lhs[k2][j2][i2][second][4][2] - coeff*lhs[k2][j2][i2][second][4][3];
  rhs[k3][j3][i3][2]   = rhs[k3][j3][i3][2]   - coeff*rhs[k3][j3][i3][3];

  coeff = lhs[k1][j1][i1][first][3][4];
  lhs[k1][j1][i1][first][4][4]= lhs[k1][j1][i1][first][4][4] - coeff*lhs[k1][j1][i1][first][4][3];
  lhs[k2][j2][i2][second][0][4] = lhs[k2][j2][i2][second][0][4] - coeff*lhs[k2][j2][i2][second][0][3];
  lhs[k2][j2][i2][second][1][4] = lhs[k2][j2][i2][second][1][4] - coeff*lhs[k2][j2][i2][second][1][3];
  lhs[k2][j2][i2][second][2][4] = lhs[k2][j2][i2][second][2][4] - coeff*lhs[k2][j2][i2][second][2][3];
  lhs[k2][j2][i2][second][3][4] = lhs[k2][j2][i2][second][3][4] - coeff*lhs[k2][j2][i2][second][3][3];
  lhs[k2][j2][i2][second][4][4] = lhs[k2][j2][i2][second][4][4] - coeff*lhs[k2][j2][i2][second][4][3];
  rhs[k3][j3][i3][4]   = rhs[k3][j3][i3][4]   - coeff*rhs[k3][j3][i3][3];


  pivot = 1.00/lhs[k1][j1][i1][first][4][4];
  lhs[k2][j2][i2][second][0][4] = lhs[k2][j2][i2][second][0][4]*pivot;
  lhs[k2][j2][i2][second][1][4] = lhs[k2][j2][i2][second][1][4]*pivot;
  lhs[k2][j2][i2][second][2][4] = lhs[k2][j2][i2][second][2][4]*pivot;
  lhs[k2][j2][i2][second][3][4] = lhs[k2][j2][i2][second][3][4]*pivot;
  lhs[k2][j2][i2][second][4][4] = lhs[k2][j2][i2][second][4][4]*pivot;
  rhs[k3][j3][i3][4]   = rhs[k3][j3][i3][4]  *pivot;

  coeff = lhs[k1][j1][i1][first][4][0];
  lhs[k2][j2][i2][second][0][0] = lhs[k2][j2][i2][second][0][0] - coeff*lhs[k2][j2][i2][second][0][4];
  lhs[k2][j2][i2][second][1][0] = lhs[k2][j2][i2][second][1][0] - coeff*lhs[k2][j2][i2][second][1][4];
  lhs[k2][j2][i2][second][2][0] = lhs[k2][j2][i2][second][2][0] - coeff*lhs[k2][j2][i2][second][2][4];
  lhs[k2][j2][i2][second][3][0] = lhs[k2][j2][i2][second][3][0] - coeff*lhs[k2][j2][i2][second][3][4];
  lhs[k2][j2][i2][second][4][0] = lhs[k2][j2][i2][second][4][0] - coeff*lhs[k2][j2][i2][second][4][4];
  rhs[k3][j3][i3][0]   = rhs[k3][j3][i3][0]   - coeff*rhs[k3][j3][i3][4];

  coeff = lhs[k1][j1][i1][first][4][1];
  lhs[k2][j2][i2][second][0][1] = lhs[k2][j2][i2][second][0][1] - coeff*lhs[k2][j2][i2][second][0][4];
  lhs[k2][j2][i2][second][1][1] = lhs[k2][j2][i2][second][1][1] - coeff*lhs[k2][j2][i2][second][1][4];
  lhs[k2][j2][i2][second][2][1] = lhs[k2][j2][i2][second][2][1] - coeff*lhs[k2][j2][i2][second][2][4];
  lhs[k2][j2][i2][second][3][1] = lhs[k2][j2][i2][second][3][1] - coeff*lhs[k2][j2][i2][second][3][4];
  lhs[k2][j2][i2][second][4][1] = lhs[k2][j2][i2][second][4][1] - coeff*lhs[k2][j2][i2][second][4][4];
  rhs[k3][j3][i3][1]   = rhs[k3][j3][i3][1]   - coeff*rhs[k3][j3][i3][4];

  coeff = lhs[k1][j1][i1][first][4][2];
  lhs[k2][j2][i2][second][0][2] = lhs[k2][j2][i2][second][0][2] - coeff*lhs[k2][j2][i2][second][0][4];
  lhs[k2][j2][i2][second][1][2] = lhs[k2][j2][i2][second][1][2] - coeff*lhs[k2][j2][i2][second][1][4];
  lhs[k2][j2][i2][second][2][2] = lhs[k2][j2][i2][second][2][2] - coeff*lhs[k2][j2][i2][second][2][4];
  lhs[k2][j2][i2][second][3][2] = lhs[k2][j2][i2][second][3][2] - coeff*lhs[k2][j2][i2][second][3][4];
  lhs[k2][j2][i2][second][4][2] = lhs[k2][j2][i2][second][4][2] - coeff*lhs[k2][j2][i2][second][4][4];
  rhs[k3][j3][i3][2]   = rhs[k3][j3][i3][2]   - coeff*rhs[k3][j3][i3][4];

  coeff = lhs[k1][j1][i1][first][4][3];
  lhs[k2][j2][i2][second][0][3] = lhs[k2][j2][i2][second][0][3] - coeff*lhs[k2][j2][i2][second][0][4];
  lhs[k2][j2][i2][second][1][3] = lhs[k2][j2][i2][second][1][3] - coeff*lhs[k2][j2][i2][second][1][4];
  lhs[k2][j2][i2][second][2][3] = lhs[k2][j2][i2][second][2][3] - coeff*lhs[k2][j2][i2][second][2][4];
  lhs[k2][j2][i2][second][3][3] = lhs[k2][j2][i2][second][3][3] - coeff*lhs[k2][j2][i2][second][3][4];
  lhs[k2][j2][i2][second][4][3] = lhs[k2][j2][i2][second][4][3] - coeff*lhs[k2][j2][i2][second][4][4];
  rhs[k3][j3][i3][3]   = rhs[k3][j3][i3][3]   - coeff*rhs[k3][j3][i3][4];
}


//#pragma acc routine
//x_binvrhs( k, j, isize, BB );//x_binvrhs( lhs[k1][j1][i1][first], rhs[k][j][i][ );
//x_binvrhs( k, j, isize, BB, k, j, isize, lhs, rhs );
//x_binvrhs( lhs[k][j][isize][BB], rhs[k][j][isize] );  lhs[k1][j1][i1][first]   rhs[k2][j2][i2]
void x_binvrhs(int k1, int j1, int i1, int first, int k2, int j2, int i2, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5], double rhs[KMAX][JMAXP+1][IMAXP+1][5])//(double lhs[5][5], double r[5])
{
  double pivot, coeff;

  pivot = 1.00/lhs[k1][j1][i1][first][0][0];
  lhs[k1][j1][i1][first][1][0] = lhs[k1][j1][i1][first][1][0]*pivot;
  lhs[k1][j1][i1][first][2][0] = lhs[k1][j1][i1][first][2][0]*pivot;
  lhs[k1][j1][i1][first][3][0] = lhs[k1][j1][i1][first][3][0]*pivot;
  lhs[k1][j1][i1][first][4][0] = lhs[k1][j1][i1][first][4][0]*pivot;
  rhs[k2][j2][i2][0]   = rhs[k2][j2][i2][0]  *pivot;

  coeff = lhs[k1][j1][i1][first][0][1];
  lhs[k1][j1][i1][first][1][1]= lhs[k1][j1][i1][first][1][1] - coeff*lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][1]= lhs[k1][j1][i1][first][2][1] - coeff*lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][1]= lhs[k1][j1][i1][first][3][1] - coeff*lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][1]= lhs[k1][j1][i1][first][4][1] - coeff*lhs[k1][j1][i1][first][4][0];
  rhs[k2][j2][i2][1]   = rhs[k2][j2][i2][1]   - coeff*rhs[k2][j2][i2][0];

  coeff = lhs[k1][j1][i1][first][0][2];
  lhs[k1][j1][i1][first][1][2]= lhs[k1][j1][i1][first][1][2] - coeff*lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][2]= lhs[k1][j1][i1][first][2][2] - coeff*lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][2]= lhs[k1][j1][i1][first][3][2] - coeff*lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][2]= lhs[k1][j1][i1][first][4][2] - coeff*lhs[k1][j1][i1][first][4][0];
  rhs[k2][j2][i2][2]   = rhs[k2][j2][i2][2]   - coeff*rhs[k2][j2][i2][0];

  coeff = lhs[k1][j1][i1][first][0][3];
  lhs[k1][j1][i1][first][1][3]= lhs[k1][j1][i1][first][1][3] - coeff*lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][3]= lhs[k1][j1][i1][first][2][3] - coeff*lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][3]= lhs[k1][j1][i1][first][3][3] - coeff*lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][3]= lhs[k1][j1][i1][first][4][3] - coeff*lhs[k1][j1][i1][first][4][0];
  rhs[k2][j2][i2][3]   = rhs[k2][j2][i2][3]   - coeff*rhs[k2][j2][i2][0];

  coeff = lhs[k1][j1][i1][first][0][4];
  lhs[k1][j1][i1][first][1][4]= lhs[k1][j1][i1][first][1][4] - coeff*lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][4]= lhs[k1][j1][i1][first][2][4] - coeff*lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][4]= lhs[k1][j1][i1][first][3][4] - coeff*lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][4]= lhs[k1][j1][i1][first][4][4] - coeff*lhs[k1][j1][i1][first][4][0];
  rhs[k2][j2][i2][4]   = rhs[k2][j2][i2][4]   - coeff*rhs[k2][j2][i2][0];


  pivot = 1.00/lhs[k1][j1][i1][first][1][1];
  lhs[k1][j1][i1][first][2][1] = lhs[k1][j1][i1][first][2][1]*pivot;
  lhs[k1][j1][i1][first][3][1] = lhs[k1][j1][i1][first][3][1]*pivot;
  lhs[k1][j1][i1][first][4][1] = lhs[k1][j1][i1][first][4][1]*pivot;
  rhs[k2][j2][i2][1]   = rhs[k2][j2][i2][1]  *pivot;

  coeff = lhs[k1][j1][i1][first][1][0];
  lhs[k1][j1][i1][first][2][0]= lhs[k1][j1][i1][first][2][0] - coeff*lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][0]= lhs[k1][j1][i1][first][3][0] - coeff*lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][0]= lhs[k1][j1][i1][first][4][0] - coeff*lhs[k1][j1][i1][first][4][1];
  rhs[k2][j2][i2][0]   = rhs[k2][j2][i2][0]   - coeff*rhs[k2][j2][i2][1];

  coeff = lhs[k1][j1][i1][first][1][2];
  lhs[k1][j1][i1][first][2][2]= lhs[k1][j1][i1][first][2][2] - coeff*lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][2]= lhs[k1][j1][i1][first][3][2] - coeff*lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][2]= lhs[k1][j1][i1][first][4][2] - coeff*lhs[k1][j1][i1][first][4][1];
  rhs[k2][j2][i2][2]   = rhs[k2][j2][i2][2]   - coeff*rhs[k2][j2][i2][1];

  coeff = lhs[k1][j1][i1][first][1][3];
  lhs[k1][j1][i1][first][2][3]= lhs[k1][j1][i1][first][2][3] - coeff*lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][3]= lhs[k1][j1][i1][first][3][3] - coeff*lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][3]= lhs[k1][j1][i1][first][4][3] - coeff*lhs[k1][j1][i1][first][4][1];
  rhs[k2][j2][i2][3]   = rhs[k2][j2][i2][3]   - coeff*rhs[k2][j2][i2][1];

  coeff = lhs[k1][j1][i1][first][1][4];
  lhs[k1][j1][i1][first][2][4]= lhs[k1][j1][i1][first][2][4] - coeff*lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][4]= lhs[k1][j1][i1][first][3][4] - coeff*lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][4]= lhs[k1][j1][i1][first][4][4] - coeff*lhs[k1][j1][i1][first][4][1];
  rhs[k2][j2][i2][4]   = rhs[k2][j2][i2][4]   - coeff*rhs[k2][j2][i2][1];


  pivot = 1.00/lhs[k1][j1][i1][first][2][2];
  lhs[k1][j1][i1][first][3][2] = lhs[k1][j1][i1][first][3][2]*pivot;
  lhs[k1][j1][i1][first][4][2] = lhs[k1][j1][i1][first][4][2]*pivot;
  rhs[k2][j2][i2][2]   = rhs[k2][j2][i2][2]  *pivot;

  coeff = lhs[k1][j1][i1][first][2][0];
  lhs[k1][j1][i1][first][3][0]= lhs[k1][j1][i1][first][3][0] - coeff*lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][0]= lhs[k1][j1][i1][first][4][0] - coeff*lhs[k1][j1][i1][first][4][2];
  rhs[k2][j2][i2][0]   = rhs[k2][j2][i2][0]   - coeff*rhs[k2][j2][i2][2];

  coeff = lhs[k1][j1][i1][first][2][1];
  lhs[k1][j1][i1][first][3][1]= lhs[k1][j1][i1][first][3][1] - coeff*lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][1]= lhs[k1][j1][i1][first][4][1] - coeff*lhs[k1][j1][i1][first][4][2];
  rhs[k2][j2][i2][1]   = rhs[k2][j2][i2][1]   - coeff*rhs[k2][j2][i2][2];

  coeff = lhs[k1][j1][i1][first][2][3];
  lhs[k1][j1][i1][first][3][3]= lhs[k1][j1][i1][first][3][3] - coeff*lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][3]= lhs[k1][j1][i1][first][4][3] - coeff*lhs[k1][j1][i1][first][4][2];
  rhs[k2][j2][i2][3]   = rhs[k2][j2][i2][3]   - coeff*rhs[k2][j2][i2][2];

  coeff = lhs[k1][j1][i1][first][2][4];
  lhs[k1][j1][i1][first][3][4]= lhs[k1][j1][i1][first][3][4] - coeff*lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][4]= lhs[k1][j1][i1][first][4][4] - coeff*lhs[k1][j1][i1][first][4][2];
  rhs[k2][j2][i2][4]   = rhs[k2][j2][i2][4]   - coeff*rhs[k2][j2][i2][2];


  pivot = 1.00/lhs[k1][j1][i1][first][3][3];
  lhs[k1][j1][i1][first][4][3] = lhs[k1][j1][i1][first][4][3]*pivot;
  rhs[k2][j2][i2][3]   = rhs[k2][j2][i2][3]  *pivot;

  coeff = lhs[k1][j1][i1][first][3][0];
  lhs[k1][j1][i1][first][4][0]= lhs[k1][j1][i1][first][4][0] - coeff*lhs[k1][j1][i1][first][4][3];
  rhs[k2][j2][i2][0]   = rhs[k2][j2][i2][0]   - coeff*rhs[k2][j2][i2][3];

  coeff = lhs[k1][j1][i1][first][3][1];
  lhs[k1][j1][i1][first][4][1]= lhs[k1][j1][i1][first][4][1] - coeff*lhs[k1][j1][i1][first][4][3];
  rhs[k2][j2][i2][1]   = rhs[k2][j2][i2][1]   - coeff*rhs[k2][j2][i2][3];

  coeff = lhs[k1][j1][i1][first][3][2];
  lhs[k1][j1][i1][first][4][2]= lhs[k1][j1][i1][first][4][2] - coeff*lhs[k1][j1][i1][first][4][3];
  rhs[k2][j2][i2][2]   = rhs[k2][j2][i2][2]   - coeff*rhs[k2][j2][i2][3];

  coeff = lhs[k1][j1][i1][first][3][4];
  lhs[k1][j1][i1][first][4][4]= lhs[k1][j1][i1][first][4][4] - coeff*lhs[k1][j1][i1][first][4][3];
  rhs[k2][j2][i2][4]   = rhs[k2][j2][i2][4]   - coeff*rhs[k2][j2][i2][3];


  pivot = 1.00/lhs[k1][j1][i1][first][4][4];
  rhs[k2][j2][i2][4]   = rhs[k2][j2][i2][4]  *pivot;

  coeff = lhs[k1][j1][i1][first][4][0];
  rhs[k2][j2][i2][0]   = rhs[k2][j2][i2][0]   - coeff*rhs[k2][j2][i2][4];

  coeff = lhs[k1][j1][i1][first][4][1];
  rhs[k2][j2][i2][1]   = rhs[k2][j2][i2][1]   - coeff*rhs[k2][j2][i2][4];

  coeff = lhs[k1][j1][i1][first][4][2];
  rhs[k2][j2][i2][2]   = rhs[k2][j2][i2][2]   - coeff*rhs[k2][j2][i2][4];

  coeff = lhs[k1][j1][i1][first][4][3];
  rhs[k2][j2][i2][3]   = rhs[k2][j2][i2][3]   - coeff*rhs[k2][j2][i2][4];
}*/

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
  #pragma acc parallel loop private(i,j,k,m,n,tmp1,tmp2,tmp3)
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
        fjac[k][j][i][0][0] = 0.0;
        fjac[k][j][i][1][0] = 1.0;
        fjac[k][j][i][2][0] = 0.0;
        fjac[k][j][i][3][0] = 0.0;
        fjac[k][j][i][4][0] = 0.0;

        fjac[k][j][i][0][1] = -(u[k][j][i][1] * tmp2 * u[k][j][i][1])
          + c2 * qs[k][j][i];
        fjac[k][j][i][1][1] = ( 2.0 - c2 ) * ( u[k][j][i][1] / u[k][j][i][0] );
        fjac[k][j][i][2][1] = - c2 * ( u[k][j][i][2] * tmp1 );
        fjac[k][j][i][3][1] = - c2 * ( u[k][j][i][3] * tmp1 );
        fjac[k][j][i][4][1] = c2;

        fjac[k][j][i][0][2] = - ( u[k][j][i][1]*u[k][j][i][2] ) * tmp2;
        fjac[k][j][i][1][2] = u[k][j][i][2] * tmp1;
        fjac[k][j][i][2][2] = u[k][j][i][1] * tmp1;
        fjac[k][j][i][3][2] = 0.0;
        fjac[k][j][i][4][2] = 0.0;

        fjac[k][j][i][0][3] = - ( u[k][j][i][1]*u[k][j][i][3] ) * tmp2;
        fjac[k][j][i][1][3] = u[k][j][i][3] * tmp1;
        fjac[k][j][i][2][3] = 0.0;
        fjac[k][j][i][3][3] = u[k][j][i][1] * tmp1;
        fjac[k][j][i][4][3] = 0.0;

        fjac[k][j][i][0][4] = ( c2 * 2.0 * square[k][j][i] - c1 * u[k][j][i][4] )
          * ( u[k][j][i][1] * tmp2 );
        fjac[k][j][i][1][4] = c1 *  u[k][j][i][4] * tmp1 
          - c2 * ( u[k][j][i][1]*u[k][j][i][1] * tmp2 + qs[k][j][i] );
        fjac[k][j][i][2][4] = - c2 * ( u[k][j][i][2]*u[k][j][i][1] ) * tmp2;
        fjac[k][j][i][3][4] = - c2 * ( u[k][j][i][3]*u[k][j][i][1] ) * tmp2;
        fjac[k][j][i][4][4] = c1 * ( u[k][j][i][1] * tmp1 );

        njac[k][j][i][0][0] = 0.0;
        njac[k][j][i][1][0] = 0.0;
        njac[k][j][i][2][0] = 0.0;
        njac[k][j][i][3][0] = 0.0;
        njac[k][j][i][4][0] = 0.0;

        njac[k][j][i][0][1] = - con43 * c3c4 * tmp2 * u[k][j][i][1];
        njac[k][j][i][1][1] =   con43 * c3c4 * tmp1;
        njac[k][j][i][2][1] =   0.0;
        njac[k][j][i][3][1] =   0.0;
        njac[k][j][i][4][1] =   0.0;

        njac[k][j][i][0][2] = - c3c4 * tmp2 * u[k][j][i][2];
        njac[k][j][i][1][2] =   0.0;
        njac[k][j][i][2][2] =   c3c4 * tmp1;
        njac[k][j][i][3][2] =   0.0;
        njac[k][j][i][4][2] =   0.0;

        njac[k][j][i][0][3] = - c3c4 * tmp2 * u[k][j][i][3];
        njac[k][j][i][1][3] =   0.0;
        njac[k][j][i][2][3] =   0.0;
        njac[k][j][i][3][3] =   c3c4 * tmp1;
        njac[k][j][i][4][3] =   0.0;

        njac[k][j][i][0][4] = - ( con43 * c3c4
            - c1345 ) * tmp3 * (u[k][j][i][1]*u[k][j][i][1])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][2]*u[k][j][i][2])
          - ( c3c4 - c1345 ) * tmp3 * (u[k][j][i][3]*u[k][j][i][3])
          - c1345 * tmp2 * u[k][j][i][4];

        njac[k][j][i][1][4] = ( con43 * c3c4
            - c1345 ) * tmp2 * u[k][j][i][1];
        njac[k][j][i][2][4] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][2];
        njac[k][j][i][3][4] = ( c3c4 - c1345 ) * tmp2 * u[k][j][i][3];
        njac[k][j][i][4][4] = ( c1345 ) * tmp1;
      }
      //---------------------------------------------------------------------
      // now jacobians set, so form left hand side in x direction
      //---------------------------------------------------------------------
      #pragma acc routine (lhsinit) worker
      lhsinit(k, j, isize, lhs);//x_lhsinit(lhs[k][j], isize);
      for (i = 1; i <= isize-1; i++) {
        tmp1 = dt * tx1;
        tmp2 = dt * tx2;

        lhs[k][j][i][AA][0][0] = - tmp2 * fjac[k][j][i-1][0][0]
          - tmp1 * njac[k][j][i-1][0][0]
          - tmp1 * dx1; 
        lhs[k][j][i][AA][1][0] = - tmp2 * fjac[k][j][i-1][1][0]
          - tmp1 * njac[k][j][i-1][1][0];
        lhs[k][j][i][AA][2][0] = - tmp2 * fjac[k][j][i-1][2][0]
          - tmp1 * njac[k][j][i-1][2][0];
        lhs[k][j][i][AA][3][0] = - tmp2 * fjac[k][j][i-1][3][0]
          - tmp1 * njac[k][j][i-1][3][0];
        lhs[k][j][i][AA][4][0] = - tmp2 * fjac[k][j][i-1][4][0]
          - tmp1 * njac[k][j][i-1][4][0];

        lhs[k][j][i][AA][0][1] = - tmp2 * fjac[k][j][i-1][0][1]
          - tmp1 * njac[k][j][i-1][0][1];
        lhs[k][j][i][AA][1][1] = - tmp2 * fjac[k][j][i-1][1][1]
          - tmp1 * njac[k][j][i-1][1][1]
          - tmp1 * dx2;
        lhs[k][j][i][AA][2][1] = - tmp2 * fjac[k][j][i-1][2][1]
          - tmp1 * njac[k][j][i-1][2][1];
        lhs[k][j][i][AA][3][1] = - tmp2 * fjac[k][j][i-1][3][1]
          - tmp1 * njac[k][j][i-1][3][1];
        lhs[k][j][i][AA][4][1] = - tmp2 * fjac[k][j][i-1][4][1]
          - tmp1 * njac[k][j][i-1][4][1];

        lhs[k][j][i][AA][0][2] = - tmp2 * fjac[k][j][i-1][0][2]
          - tmp1 * njac[k][j][i-1][0][2];
        lhs[k][j][i][AA][1][2] = - tmp2 * fjac[k][j][i-1][1][2]
          - tmp1 * njac[k][j][i-1][1][2];
        lhs[k][j][i][AA][2][2] = - tmp2 * fjac[k][j][i-1][2][2]
          - tmp1 * njac[k][j][i-1][2][2]
          - tmp1 * dx3;
        lhs[k][j][i][AA][3][2] = - tmp2 * fjac[k][j][i-1][3][2]
          - tmp1 * njac[k][j][i-1][3][2];
        lhs[k][j][i][AA][4][2] = - tmp2 * fjac[k][j][i-1][4][2]
          - tmp1 * njac[k][j][i-1][4][2];

        lhs[k][j][i][AA][0][3] = - tmp2 * fjac[k][j][i-1][0][3]
          - tmp1 * njac[k][j][i-1][0][3];
        lhs[k][j][i][AA][1][3] = - tmp2 * fjac[k][j][i-1][1][3]
          - tmp1 * njac[k][j][i-1][1][3];
        lhs[k][j][i][AA][2][3] = - tmp2 * fjac[k][j][i-1][2][3]
          - tmp1 * njac[k][j][i-1][2][3];
        lhs[k][j][i][AA][3][3] = - tmp2 * fjac[k][j][i-1][3][3]
          - tmp1 * njac[k][j][i-1][3][3]
          - tmp1 * dx4;
        lhs[k][j][i][AA][4][3] = - tmp2 * fjac[k][j][i-1][4][3]
          - tmp1 * njac[k][j][i-1][4][3];

        lhs[k][j][i][AA][0][4] = - tmp2 * fjac[k][j][i-1][0][4]
          - tmp1 * njac[k][j][i-1][0][4];
        lhs[k][j][i][AA][1][4] = - tmp2 * fjac[k][j][i-1][1][4]
          - tmp1 * njac[k][j][i-1][1][4];
        lhs[k][j][i][AA][2][4] = - tmp2 * fjac[k][j][i-1][2][4]
          - tmp1 * njac[k][j][i-1][2][4];
        lhs[k][j][i][AA][3][4] = - tmp2 * fjac[k][j][i-1][3][4]
          - tmp1 * njac[k][j][i-1][3][4];
        lhs[k][j][i][AA][4][4] = - tmp2 * fjac[k][j][i-1][4][4]
          - tmp1 * njac[k][j][i-1][4][4]
          - tmp1 * dx5;

        lhs[k][j][i][BB][0][0] = 1.0
          + tmp1 * 2.0 * njac[k][j][i][0][0]
          + tmp1 * 2.0 * dx1;
        lhs[k][j][i][BB][1][0] = tmp1 * 2.0 * njac[k][j][i][1][0];
        lhs[k][j][i][BB][2][0] = tmp1 * 2.0 * njac[k][j][i][2][0];
        lhs[k][j][i][BB][3][0] = tmp1 * 2.0 * njac[k][j][i][3][0];
        lhs[k][j][i][BB][4][0] = tmp1 * 2.0 * njac[k][j][i][4][0];

        lhs[k][j][i][BB][0][1] = tmp1 * 2.0 * njac[k][j][i][0][1];
        lhs[k][j][i][BB][1][1] = 1.0
          + tmp1 * 2.0 * njac[k][j][i][1][1]
          + tmp1 * 2.0 * dx2;
        lhs[k][j][i][BB][2][1] = tmp1 * 2.0 * njac[k][j][i][2][1];
        lhs[k][j][i][BB][3][1] = tmp1 * 2.0 * njac[k][j][i][3][1];
        lhs[k][j][i][BB][4][1] = tmp1 * 2.0 * njac[k][j][i][4][1];

        lhs[k][j][i][BB][0][2] = tmp1 * 2.0 * njac[k][j][i][0][2];
        lhs[k][j][i][BB][1][2] = tmp1 * 2.0 * njac[k][j][i][1][2];
        lhs[k][j][i][BB][2][2] = 1.0
          + tmp1 * 2.0 * njac[k][j][i][2][2]
          + tmp1 * 2.0 * dx3;
        lhs[k][j][i][BB][3][2] = tmp1 * 2.0 * njac[k][j][i][3][2];
        lhs[k][j][i][BB][4][2] = tmp1 * 2.0 * njac[k][j][i][4][2];

        lhs[k][j][i][BB][0][3] = tmp1 * 2.0 * njac[k][j][i][0][3];
        lhs[k][j][i][BB][1][3] = tmp1 * 2.0 * njac[k][j][i][1][3];
        lhs[k][j][i][BB][2][3] = tmp1 * 2.0 * njac[k][j][i][2][3];
        lhs[k][j][i][BB][3][3] = 1.0
          + tmp1 * 2.0 * njac[k][j][i][3][3]
          + tmp1 * 2.0 * dx4;
        lhs[k][j][i][BB][4][3] = tmp1 * 2.0 * njac[k][j][i][4][3];

        lhs[k][j][i][BB][0][4] = tmp1 * 2.0 * njac[k][j][i][0][4];
        lhs[k][j][i][BB][1][4] = tmp1 * 2.0 * njac[k][j][i][1][4];
        lhs[k][j][i][BB][2][4] = tmp1 * 2.0 * njac[k][j][i][2][4];
        lhs[k][j][i][BB][3][4] = tmp1 * 2.0 * njac[k][j][i][3][4];
        lhs[k][j][i][BB][4][4] = 1.0
          + tmp1 * 2.0 * njac[k][j][i][4][4]
          + tmp1 * 2.0 * dx5;

        lhs[k][j][i][CC][0][0] =  tmp2 * fjac[k][j][i+1][0][0]
          - tmp1 * njac[k][j][i+1][0][0]
          - tmp1 * dx1;
        lhs[k][j][i][CC][1][0] =  tmp2 * fjac[k][j][i+1][1][0]
          - tmp1 * njac[k][j][i+1][1][0];
        lhs[k][j][i][CC][2][0] =  tmp2 * fjac[k][j][i+1][2][0]
          - tmp1 * njac[k][j][i+1][2][0];
        lhs[k][j][i][CC][3][0] =  tmp2 * fjac[k][j][i+1][3][0]
          - tmp1 * njac[k][j][i+1][3][0];
        lhs[k][j][i][CC][4][0] =  tmp2 * fjac[k][j][i+1][4][0]
          - tmp1 * njac[k][j][i+1][4][0];

        lhs[k][j][i][CC][0][1] =  tmp2 * fjac[k][j][i+1][0][1]
          - tmp1 * njac[k][j][i+1][0][1];
        lhs[k][j][i][CC][1][1] =  tmp2 * fjac[k][j][i+1][1][1]
          - tmp1 * njac[k][j][i+1][1][1]
          - tmp1 * dx2;
        lhs[k][j][i][CC][2][1] =  tmp2 * fjac[k][j][i+1][2][1]
          - tmp1 * njac[k][j][i+1][2][1];
        lhs[k][j][i][CC][3][1] =  tmp2 * fjac[k][j][i+1][3][1]
          - tmp1 * njac[k][j][i+1][3][1];
        lhs[k][j][i][CC][4][1] =  tmp2 * fjac[k][j][i+1][4][1]
          - tmp1 * njac[k][j][i+1][4][1];

        lhs[k][j][i][CC][0][2] =  tmp2 * fjac[k][j][i+1][0][2]
          - tmp1 * njac[k][j][i+1][0][2];
        lhs[k][j][i][CC][1][2] =  tmp2 * fjac[k][j][i+1][1][2]
          - tmp1 * njac[k][j][i+1][1][2];
        lhs[k][j][i][CC][2][2] =  tmp2 * fjac[k][j][i+1][2][2]
          - tmp1 * njac[k][j][i+1][2][2]
          - tmp1 * dx3;
        lhs[k][j][i][CC][3][2] =  tmp2 * fjac[k][j][i+1][3][2]
          - tmp1 * njac[k][j][i+1][3][2];
        lhs[k][j][i][CC][4][2] =  tmp2 * fjac[k][j][i+1][4][2]
          - tmp1 * njac[k][j][i+1][4][2];

        lhs[k][j][i][CC][0][3] =  tmp2 * fjac[k][j][i+1][0][3]
          - tmp1 * njac[k][j][i+1][0][3];
        lhs[k][j][i][CC][1][3] =  tmp2 * fjac[k][j][i+1][1][3]
          - tmp1 * njac[k][j][i+1][1][3];
        lhs[k][j][i][CC][2][3] =  tmp2 * fjac[k][j][i+1][2][3]
          - tmp1 * njac[k][j][i+1][2][3];
        lhs[k][j][i][CC][3][3] =  tmp2 * fjac[k][j][i+1][3][3]
          - tmp1 * njac[k][j][i+1][3][3]
          - tmp1 * dx4;
        lhs[k][j][i][CC][4][3] =  tmp2 * fjac[k][j][i+1][4][3]
          - tmp1 * njac[k][j][i+1][4][3];

        lhs[k][j][i][CC][0][4] =  tmp2 * fjac[k][j][i+1][0][4]
          - tmp1 * njac[k][j][i+1][0][4];
        lhs[k][j][i][CC][1][4] =  tmp2 * fjac[k][j][i+1][1][4]
          - tmp1 * njac[k][j][i+1][1][4];
        lhs[k][j][i][CC][2][4] =  tmp2 * fjac[k][j][i+1][2][4]
          - tmp1 * njac[k][j][i+1][2][4];
        lhs[k][j][i][CC][3][4] =  tmp2 * fjac[k][j][i+1][3][4]
          - tmp1 * njac[k][j][i+1][3][4];
        lhs[k][j][i][CC][4][4] =  tmp2 * fjac[k][j][i+1][4][4]
          - tmp1 * njac[k][j][i+1][4][4]
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
      #pragma acc routine (binvcrhs) worker
      binvcrhs( k, j, 0, BB, k, j, 0, CC, k, j, 0, lhs, rhs );//x_binvcrhs( lhs[k][j][0][BB], lhs[k][j][0][CC], rhs[k][j][0] );

      //---------------------------------------------------------------------
      // begin inner most do loop
      // do all the elements of the cell unless last 
      //---------------------------------------------------------------------
      for (i = 1; i <= isize-1; i++) {
        //-------------------------------------------------------------------
        // rhs(i) = rhs(i) - A*rhs(i-1)
        //-------------------------------------------------------------------
        #pragma acc routine (matvec_sub) worker
        matvec_sub(k, j, i, AA, k, j, i-1, k, j, i, lhs, rhs);//x_matvec_sub(lhs[k][j][i][AA], rhs[k][j][i-1], rhs[k][j][i]);

        //-------------------------------------------------------------------
        // B(i) = B(i) - C(i-1)*A(i)
        //-------------------------------------------------------------------
        #pragma acc routine (matmul_sub) worker
        matmul_sub(k, j, i, AA, k, j, i-1, CC, k, j, i, BB, lhs);//x_matmul_sub(lhs[k][j][i][AA], lhs[k][j][i-1][CC], lhs[k][j][i][BB]);


        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[k][j][0] by b_inverse[k][j][0] and copy to rhs
        //-------------------------------------------------------------------
        #pragma acc routine (binvcrhs) worker
        binvcrhs( k, j, i, BB, k, j, i, CC, k, j, i, lhs, rhs );//x_binvcrhs( lhs[k][j][i][BB], lhs[k][j][i][CC], rhs[k][j][i] );
      }

      //---------------------------------------------------------------------
      // rhs(isize) = rhs(isize) - A*rhs(isize-1)
      //---------------------------------------------------------------------
      #pragma acc routine (matvec_sub) worker
      matvec_sub(k, j, isize, AA, k, j, isize-1, k, j, isize, lhs, rhs);//x_matvec_sub(lhs[k][j][isize][AA], rhs[k][j][isize-1], rhs[k][j][isize]);

      //---------------------------------------------------------------------
      // B(isize) = B(isize) - C(isize-1)*A(isize)
      //---------------------------------------------------------------------
      #pragma acc routine (matmul_sub) worker
      matmul_sub(k, j, isize, AA, k, j, isize-1, CC, k, j, isize, BB, lhs);//x_matmul_sub(lhs[k][j][isize][AA], lhs[k][j][isize-1][CC], lhs[k][j][isize][BB]);

      //---------------------------------------------------------------------
      // multiply rhs() by b_inverse() and copy to rhs
      //---------------------------------------------------------------------
      #pragma acc routine (binvrhs) worker
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
              - lhs[k][j][i][CC][n][m]*rhs[k][j][i+1][n];
          }
        }
      }
    }
  }
  if (timeron) timer_stop(t_xsolve);
}
