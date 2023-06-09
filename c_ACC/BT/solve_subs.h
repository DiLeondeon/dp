#include "header.h"

void lhsinit(int k, int j, int ni, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5])//(double lhs[][3][5][5], int ni)
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

void matvec_sub(int k1, int j1, int i1, int first, int k2, int j2, int i2, int k3, int j3, int i3, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5], double rhs[KMAX][JMAXP+1][IMAXP+1][5])//(double ablock[5][5], double avec[5], double bvec[5])//x_matvec_sub(k, j, i, AA, i-1);//x_matvec_sub(lhs[k][j][i][first], rhs[k][j][second], rhs[k][j][i]);
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
void matmul_sub(int k1, int j1, int i1, int first, int k2, int j2, int i2, int second, int k3, int j3, int i3, int third, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5])//(double ablock[5][5], double bblock[5][5], double cblock[5][5])
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
void binvcrhs(int k1, int j1, int i1, int first, int k2, int j2, int i2, int second, int k3, int j3, int i3, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5], double rhs[KMAX][JMAXP+1][IMAXP+1][5])//(double lhs[5][5], double c[5][5], double r[5]) //x_binvcrhs( lhs[k][j][i][BB], lhs[k][j][i][CC][, rhs[k][j][i][ );
{
  double pivot, coeff;
  int n,m,m1;

  for (m=0;m<5;m++){
    pivot = 1.00/lhs[k1][j1][i1][first][m][m];
    for (n=m+1;n<5;n++){
      lhs[k1][j1][i1][first][n][m] = lhs[k1][j1][i1][first][n][m]*pivot;
    }
    for (n=0;n<5;n++){
      lhs[k2][j2][i2][second][n][m] = lhs[k2][j2][i2][second][n][m]*pivot;
    }
    rhs[k3][j3][i3][m]   = rhs[k3][j3][i3][m]  *pivot;
    for (m1=0;m1<5;m1++){
      if (m1 != m){
        coeff = lhs[k1][j1][i1][first][m][m1];
        for (n=m+1;n<5;n++){
          lhs[k1][j1][i1][first][n][m1]= lhs[k1][j1][i1][first][n][m1] - coeff*lhs[k1][j1][i1][first][n][m];
        }
        for (n=0;n<5;n++){
          lhs[k2][j2][i2][second][n][m1] = lhs[k2][j2][i2][second][n][m1] - coeff*lhs[k2][j2][i2][second][n][m];
        }
        rhs[k3][j3][i3][m1]   = rhs[k3][j3][i3][m1]   - coeff*rhs[k3][j3][i3][m];
      }
    }
  }
  /*pivot = 1.00/lhs[k1][j1][i1][first][0][0];
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
  rhs[k3][j3][i3][3]   = rhs[k3][j3][i3][3]   - coeff*rhs[k3][j3][i3][4];*/
}


//#pragma acc routine
//x_binvrhs( k, j, isize, BB );//x_binvrhs( lhs[k1][j1][i1][first], rhs[k][j][i][ );
//x_binvrhs( k, j, isize, BB, k, j, isize, lhs, rhs );
//x_binvrhs( lhs[k][j][isize][BB], rhs[k][j][isize] );  lhs[k1][j1][i1][first]   rhs[k2][j2][i2]
void binvrhs(int k1, int j1, int i1, int first, int k2, int j2, int i2, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5], double rhs[KMAX][JMAXP+1][IMAXP+1][5])//(double lhs[5][5], double r[5])
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
}



/*for (m=0;m<5:m++){
  pivot = 1.00/lhs[k1][j1][i1][first][m][m];
  for (n=m+1;n<5;n++){
    lhs[k1][j1][i1][first][n][m] = lhs[k1][j1][i1][first][n][m]*pivot;
  }
  for (n=0;n<5;n++){
    lhs[k2][j2][i2][second][n][m] = lhs[k2][j2][i2][second][n][m]*pivot;
  }
  rhs[k3][j3][i3][m]   = rhs[k3][j3][i3][m]  *pivot;
  for (m1=0;m1<5;m1++){
    if (m1 != m){
      coeff = lhs[k1][j1][i1][first][m][m1];
      for (n=m+1;n<5;n++){
        lhs[k1][j1][i1][first][n][m1]= lhs[k1][j1][i1][first][n][m1] - coeff*lhs[k1][j1][i1][first][n][m];
      }
      for (n=0;n<5;n++){
        lhs[k2][j2][i2][second][n][m1] = lhs[k2][j2][i2][second][n][m1] - coeff*lhs[k2][j2][i2][second][n][m];
      }
      rhs[k3][j3][i3][m1]   = rhs[k3][j3][i3][m1]   - coeff*rhs[k3][j3][i3][m];
    }
  }
}*/