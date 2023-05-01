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
void y_lhsinit(double lhs[][3][5][5], int ni)
{
  int i, m, n;

  //---------------------------------------------------------------------
  // zero the whole left hand side for starters
  // set all diagonal values to 1. This is overkill, but convenient
  //---------------------------------------------------------------------
  i = 0;
  for (n = 0; n < 5; n++) {
    for (m = 0; m < 5; m++) {
      lhs[i][0][n][m] = 0.0;
      lhs[i][1][n][m] = 0.0;
      lhs[i][2][n][m] = 0.0;
    }
    lhs[i][1][n][n] = 1.0;
  }
  i = ni;
  for (n = 0; n < 5; n++) {
    for (m = 0; m < 5; m++) {
      lhs[i][0][n][m] = 0.0;
      lhs[i][1][n][m] = 0.0;
      lhs[i][2][n][m] = 0.0;
    }
    lhs[i][1][n][n] = 1.0;
  }
}

//#pragma acc routine
void y_matvec_sub(double ablock[5][5], double avec[5], double bvec[5])
{
  //---------------------------------------------------------------------
  // rhs[kc][jc][ic][i] = rhs[kc][jc][ic][i] 
  // $                  - lhs[ia][ablock][0][i]*
  //---------------------------------------------------------------------
  bvec[0] = bvec[0] - ablock[0][0]*avec[0]
                    - ablock[1][0]*avec[1]
                    - ablock[2][0]*avec[2]
                    - ablock[3][0]*avec[3]
                    - ablock[4][0]*avec[4];
  bvec[1] = bvec[1] - ablock[0][1]*avec[0]
                    - ablock[1][1]*avec[1]
                    - ablock[2][1]*avec[2]
                    - ablock[3][1]*avec[3]
                    - ablock[4][1]*avec[4];
  bvec[2] = bvec[2] - ablock[0][2]*avec[0]
                    - ablock[1][2]*avec[1]
                    - ablock[2][2]*avec[2]
                    - ablock[3][2]*avec[3]
                    - ablock[4][2]*avec[4];
  bvec[3] = bvec[3] - ablock[0][3]*avec[0]
                    - ablock[1][3]*avec[1]
                    - ablock[2][3]*avec[2]
                    - ablock[3][3]*avec[3]
                    - ablock[4][3]*avec[4];
  bvec[4] = bvec[4] - ablock[0][4]*avec[0]
                    - ablock[1][4]*avec[1]
                    - ablock[2][4]*avec[2]
                    - ablock[3][4]*avec[3]
                    - ablock[4][4]*avec[4];
}


//---------------------------------------------------------------------
// subtracts a(i,j,k) X b(i,j,k) from c(i,j,k)
//---------------------------------------------------------------------
//#pragma acc routine
void y_matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5])
{
  cblock[0][0] = cblock[0][0] - ablock[0][0]*bblock[0][0]
                              - ablock[1][0]*bblock[0][1]
                              - ablock[2][0]*bblock[0][2]
                              - ablock[3][0]*bblock[0][3]
                              - ablock[4][0]*bblock[0][4];
  cblock[0][1] = cblock[0][1] - ablock[0][1]*bblock[0][0]
                              - ablock[1][1]*bblock[0][1]
                              - ablock[2][1]*bblock[0][2]
                              - ablock[3][1]*bblock[0][3]
                              - ablock[4][1]*bblock[0][4];
  cblock[0][2] = cblock[0][2] - ablock[0][2]*bblock[0][0]
                              - ablock[1][2]*bblock[0][1]
                              - ablock[2][2]*bblock[0][2]
                              - ablock[3][2]*bblock[0][3]
                              - ablock[4][2]*bblock[0][4];
  cblock[0][3] = cblock[0][3] - ablock[0][3]*bblock[0][0]
                              - ablock[1][3]*bblock[0][1]
                              - ablock[2][3]*bblock[0][2]
                              - ablock[3][3]*bblock[0][3]
                              - ablock[4][3]*bblock[0][4];
  cblock[0][4] = cblock[0][4] - ablock[0][4]*bblock[0][0]
                              - ablock[1][4]*bblock[0][1]
                              - ablock[2][4]*bblock[0][2]
                              - ablock[3][4]*bblock[0][3]
                              - ablock[4][4]*bblock[0][4];
  cblock[1][0] = cblock[1][0] - ablock[0][0]*bblock[1][0]
                              - ablock[1][0]*bblock[1][1]
                              - ablock[2][0]*bblock[1][2]
                              - ablock[3][0]*bblock[1][3]
                              - ablock[4][0]*bblock[1][4];
  cblock[1][1] = cblock[1][1] - ablock[0][1]*bblock[1][0]
                              - ablock[1][1]*bblock[1][1]
                              - ablock[2][1]*bblock[1][2]
                              - ablock[3][1]*bblock[1][3]
                              - ablock[4][1]*bblock[1][4];
  cblock[1][2] = cblock[1][2] - ablock[0][2]*bblock[1][0]
                              - ablock[1][2]*bblock[1][1]
                              - ablock[2][2]*bblock[1][2]
                              - ablock[3][2]*bblock[1][3]
                              - ablock[4][2]*bblock[1][4];
  cblock[1][3] = cblock[1][3] - ablock[0][3]*bblock[1][0]
                              - ablock[1][3]*bblock[1][1]
                              - ablock[2][3]*bblock[1][2]
                              - ablock[3][3]*bblock[1][3]
                              - ablock[4][3]*bblock[1][4];
  cblock[1][4] = cblock[1][4] - ablock[0][4]*bblock[1][0]
                              - ablock[1][4]*bblock[1][1]
                              - ablock[2][4]*bblock[1][2]
                              - ablock[3][4]*bblock[1][3]
                              - ablock[4][4]*bblock[1][4];
  cblock[2][0] = cblock[2][0] - ablock[0][0]*bblock[2][0]
                              - ablock[1][0]*bblock[2][1]
                              - ablock[2][0]*bblock[2][2]
                              - ablock[3][0]*bblock[2][3]
                              - ablock[4][0]*bblock[2][4];
  cblock[2][1] = cblock[2][1] - ablock[0][1]*bblock[2][0]
                              - ablock[1][1]*bblock[2][1]
                              - ablock[2][1]*bblock[2][2]
                              - ablock[3][1]*bblock[2][3]
                              - ablock[4][1]*bblock[2][4];
  cblock[2][2] = cblock[2][2] - ablock[0][2]*bblock[2][0]
                              - ablock[1][2]*bblock[2][1]
                              - ablock[2][2]*bblock[2][2]
                              - ablock[3][2]*bblock[2][3]
                              - ablock[4][2]*bblock[2][4];
  cblock[2][3] = cblock[2][3] - ablock[0][3]*bblock[2][0]
                              - ablock[1][3]*bblock[2][1]
                              - ablock[2][3]*bblock[2][2]
                              - ablock[3][3]*bblock[2][3]
                              - ablock[4][3]*bblock[2][4];
  cblock[2][4] = cblock[2][4] - ablock[0][4]*bblock[2][0]
                              - ablock[1][4]*bblock[2][1]
                              - ablock[2][4]*bblock[2][2]
                              - ablock[3][4]*bblock[2][3]
                              - ablock[4][4]*bblock[2][4];
  cblock[3][0] = cblock[3][0] - ablock[0][0]*bblock[3][0]
                              - ablock[1][0]*bblock[3][1]
                              - ablock[2][0]*bblock[3][2]
                              - ablock[3][0]*bblock[3][3]
                              - ablock[4][0]*bblock[3][4];
  cblock[3][1] = cblock[3][1] - ablock[0][1]*bblock[3][0]
                              - ablock[1][1]*bblock[3][1]
                              - ablock[2][1]*bblock[3][2]
                              - ablock[3][1]*bblock[3][3]
                              - ablock[4][1]*bblock[3][4];
  cblock[3][2] = cblock[3][2] - ablock[0][2]*bblock[3][0]
                              - ablock[1][2]*bblock[3][1]
                              - ablock[2][2]*bblock[3][2]
                              - ablock[3][2]*bblock[3][3]
                              - ablock[4][2]*bblock[3][4];
  cblock[3][3] = cblock[3][3] - ablock[0][3]*bblock[3][0]
                              - ablock[1][3]*bblock[3][1]
                              - ablock[2][3]*bblock[3][2]
                              - ablock[3][3]*bblock[3][3]
                              - ablock[4][3]*bblock[3][4];
  cblock[3][4] = cblock[3][4] - ablock[0][4]*bblock[3][0]
                              - ablock[1][4]*bblock[3][1]
                              - ablock[2][4]*bblock[3][2]
                              - ablock[3][4]*bblock[3][3]
                              - ablock[4][4]*bblock[3][4];
  cblock[4][0] = cblock[4][0] - ablock[0][0]*bblock[4][0]
                              - ablock[1][0]*bblock[4][1]
                              - ablock[2][0]*bblock[4][2]
                              - ablock[3][0]*bblock[4][3]
                              - ablock[4][0]*bblock[4][4];
  cblock[4][1] = cblock[4][1] - ablock[0][1]*bblock[4][0]
                              - ablock[1][1]*bblock[4][1]
                              - ablock[2][1]*bblock[4][2]
                              - ablock[3][1]*bblock[4][3]
                              - ablock[4][1]*bblock[4][4];
  cblock[4][2] = cblock[4][2] - ablock[0][2]*bblock[4][0]
                              - ablock[1][2]*bblock[4][1]
                              - ablock[2][2]*bblock[4][2]
                              - ablock[3][2]*bblock[4][3]
                              - ablock[4][2]*bblock[4][4];
  cblock[4][3] = cblock[4][3] - ablock[0][3]*bblock[4][0]
                              - ablock[1][3]*bblock[4][1]
                              - ablock[2][3]*bblock[4][2]
                              - ablock[3][3]*bblock[4][3]
                              - ablock[4][3]*bblock[4][4];
  cblock[4][4] = cblock[4][4] - ablock[0][4]*bblock[4][0]
                              - ablock[1][4]*bblock[4][1]
                              - ablock[2][4]*bblock[4][2]
                              - ablock[3][4]*bblock[4][3]
                              - ablock[4][4]*bblock[4][4];
}


//#pragma acc routine
void y_binvcrhs(double lhs[5][5], double c[5][5], double r[5])
{
  double pivot, coeff;

  pivot = 1.00/lhs[0][0];
  lhs[1][0] = lhs[1][0]*pivot;
  lhs[2][0] = lhs[2][0]*pivot;
  lhs[3][0] = lhs[3][0]*pivot;
  lhs[4][0] = lhs[4][0]*pivot;
  c[0][0] = c[0][0]*pivot;
  c[1][0] = c[1][0]*pivot;
  c[2][0] = c[2][0]*pivot;
  c[3][0] = c[3][0]*pivot;
  c[4][0] = c[4][0]*pivot;
  r[0]   = r[0]  *pivot;

  coeff = lhs[0][1];
  lhs[1][1]= lhs[1][1] - coeff*lhs[1][0];
  lhs[2][1]= lhs[2][1] - coeff*lhs[2][0];
  lhs[3][1]= lhs[3][1] - coeff*lhs[3][0];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][0];
  c[0][1] = c[0][1] - coeff*c[0][0];
  c[1][1] = c[1][1] - coeff*c[1][0];
  c[2][1] = c[2][1] - coeff*c[2][0];
  c[3][1] = c[3][1] - coeff*c[3][0];
  c[4][1] = c[4][1] - coeff*c[4][0];
  r[1]   = r[1]   - coeff*r[0];

  coeff = lhs[0][2];
  lhs[1][2]= lhs[1][2] - coeff*lhs[1][0];
  lhs[2][2]= lhs[2][2] - coeff*lhs[2][0];
  lhs[3][2]= lhs[3][2] - coeff*lhs[3][0];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][0];
  c[0][2] = c[0][2] - coeff*c[0][0];
  c[1][2] = c[1][2] - coeff*c[1][0];
  c[2][2] = c[2][2] - coeff*c[2][0];
  c[3][2] = c[3][2] - coeff*c[3][0];
  c[4][2] = c[4][2] - coeff*c[4][0];
  r[2]   = r[2]   - coeff*r[0];

  coeff = lhs[0][3];
  lhs[1][3]= lhs[1][3] - coeff*lhs[1][0];
  lhs[2][3]= lhs[2][3] - coeff*lhs[2][0];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][0];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][0];
  c[0][3] = c[0][3] - coeff*c[0][0];
  c[1][3] = c[1][3] - coeff*c[1][0];
  c[2][3] = c[2][3] - coeff*c[2][0];
  c[3][3] = c[3][3] - coeff*c[3][0];
  c[4][3] = c[4][3] - coeff*c[4][0];
  r[3]   = r[3]   - coeff*r[0];

  coeff = lhs[0][4];
  lhs[1][4]= lhs[1][4] - coeff*lhs[1][0];
  lhs[2][4]= lhs[2][4] - coeff*lhs[2][0];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][0];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][0];
  c[0][4] = c[0][4] - coeff*c[0][0];
  c[1][4] = c[1][4] - coeff*c[1][0];
  c[2][4] = c[2][4] - coeff*c[2][0];
  c[3][4] = c[3][4] - coeff*c[3][0];
  c[4][4] = c[4][4] - coeff*c[4][0];
  r[4]   = r[4]   - coeff*r[0];


  pivot = 1.00/lhs[1][1];
  lhs[2][1] = lhs[2][1]*pivot;
  lhs[3][1] = lhs[3][1]*pivot;
  lhs[4][1] = lhs[4][1]*pivot;
  c[0][1] = c[0][1]*pivot;
  c[1][1] = c[1][1]*pivot;
  c[2][1] = c[2][1]*pivot;
  c[3][1] = c[3][1]*pivot;
  c[4][1] = c[4][1]*pivot;
  r[1]   = r[1]  *pivot;

  coeff = lhs[1][0];
  lhs[2][0]= lhs[2][0] - coeff*lhs[2][1];
  lhs[3][0]= lhs[3][0] - coeff*lhs[3][1];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][1];
  c[0][0] = c[0][0] - coeff*c[0][1];
  c[1][0] = c[1][0] - coeff*c[1][1];
  c[2][0] = c[2][0] - coeff*c[2][1];
  c[3][0] = c[3][0] - coeff*c[3][1];
  c[4][0] = c[4][0] - coeff*c[4][1];
  r[0]   = r[0]   - coeff*r[1];

  coeff = lhs[1][2];
  lhs[2][2]= lhs[2][2] - coeff*lhs[2][1];
  lhs[3][2]= lhs[3][2] - coeff*lhs[3][1];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][1];
  c[0][2] = c[0][2] - coeff*c[0][1];
  c[1][2] = c[1][2] - coeff*c[1][1];
  c[2][2] = c[2][2] - coeff*c[2][1];
  c[3][2] = c[3][2] - coeff*c[3][1];
  c[4][2] = c[4][2] - coeff*c[4][1];
  r[2]   = r[2]   - coeff*r[1];

  coeff = lhs[1][3];
  lhs[2][3]= lhs[2][3] - coeff*lhs[2][1];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][1];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][1];
  c[0][3] = c[0][3] - coeff*c[0][1];
  c[1][3] = c[1][3] - coeff*c[1][1];
  c[2][3] = c[2][3] - coeff*c[2][1];
  c[3][3] = c[3][3] - coeff*c[3][1];
  c[4][3] = c[4][3] - coeff*c[4][1];
  r[3]   = r[3]   - coeff*r[1];

  coeff = lhs[1][4];
  lhs[2][4]= lhs[2][4] - coeff*lhs[2][1];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][1];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][1];
  c[0][4] = c[0][4] - coeff*c[0][1];
  c[1][4] = c[1][4] - coeff*c[1][1];
  c[2][4] = c[2][4] - coeff*c[2][1];
  c[3][4] = c[3][4] - coeff*c[3][1];
  c[4][4] = c[4][4] - coeff*c[4][1];
  r[4]   = r[4]   - coeff*r[1];


  pivot = 1.00/lhs[2][2];
  lhs[3][2] = lhs[3][2]*pivot;
  lhs[4][2] = lhs[4][2]*pivot;
  c[0][2] = c[0][2]*pivot;
  c[1][2] = c[1][2]*pivot;
  c[2][2] = c[2][2]*pivot;
  c[3][2] = c[3][2]*pivot;
  c[4][2] = c[4][2]*pivot;
  r[2]   = r[2]  *pivot;

  coeff = lhs[2][0];
  lhs[3][0]= lhs[3][0] - coeff*lhs[3][2];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][2];
  c[0][0] = c[0][0] - coeff*c[0][2];
  c[1][0] = c[1][0] - coeff*c[1][2];
  c[2][0] = c[2][0] - coeff*c[2][2];
  c[3][0] = c[3][0] - coeff*c[3][2];
  c[4][0] = c[4][0] - coeff*c[4][2];
  r[0]   = r[0]   - coeff*r[2];

  coeff = lhs[2][1];
  lhs[3][1]= lhs[3][1] - coeff*lhs[3][2];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][2];
  c[0][1] = c[0][1] - coeff*c[0][2];
  c[1][1] = c[1][1] - coeff*c[1][2];
  c[2][1] = c[2][1] - coeff*c[2][2];
  c[3][1] = c[3][1] - coeff*c[3][2];
  c[4][1] = c[4][1] - coeff*c[4][2];
  r[1]   = r[1]   - coeff*r[2];

  coeff = lhs[2][3];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][2];
  c[0][3] = c[0][3] - coeff*c[0][2];
  c[1][3] = c[1][3] - coeff*c[1][2];
  c[2][3] = c[2][3] - coeff*c[2][2];
  c[3][3] = c[3][3] - coeff*c[3][2];
  c[4][3] = c[4][3] - coeff*c[4][2];
  r[3]   = r[3]   - coeff*r[2];

  coeff = lhs[2][4];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][2];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][2];
  c[0][4] = c[0][4] - coeff*c[0][2];
  c[1][4] = c[1][4] - coeff*c[1][2];
  c[2][4] = c[2][4] - coeff*c[2][2];
  c[3][4] = c[3][4] - coeff*c[3][2];
  c[4][4] = c[4][4] - coeff*c[4][2];
  r[4]   = r[4]   - coeff*r[2];


  pivot = 1.00/lhs[3][3];
  lhs[4][3] = lhs[4][3]*pivot;
  c[0][3] = c[0][3]*pivot;
  c[1][3] = c[1][3]*pivot;
  c[2][3] = c[2][3]*pivot;
  c[3][3] = c[3][3]*pivot;
  c[4][3] = c[4][3]*pivot;
  r[3]   = r[3]  *pivot;

  coeff = lhs[3][0];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][3];
  c[0][0] = c[0][0] - coeff*c[0][3];
  c[1][0] = c[1][0] - coeff*c[1][3];
  c[2][0] = c[2][0] - coeff*c[2][3];
  c[3][0] = c[3][0] - coeff*c[3][3];
  c[4][0] = c[4][0] - coeff*c[4][3];
  r[0]   = r[0]   - coeff*r[3];

  coeff = lhs[3][1];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][3];
  c[0][1] = c[0][1] - coeff*c[0][3];
  c[1][1] = c[1][1] - coeff*c[1][3];
  c[2][1] = c[2][1] - coeff*c[2][3];
  c[3][1] = c[3][1] - coeff*c[3][3];
  c[4][1] = c[4][1] - coeff*c[4][3];
  r[1]   = r[1]   - coeff*r[3];

  coeff = lhs[3][2];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][3];
  c[0][2] = c[0][2] - coeff*c[0][3];
  c[1][2] = c[1][2] - coeff*c[1][3];
  c[2][2] = c[2][2] - coeff*c[2][3];
  c[3][2] = c[3][2] - coeff*c[3][3];
  c[4][2] = c[4][2] - coeff*c[4][3];
  r[2]   = r[2]   - coeff*r[3];

  coeff = lhs[3][4];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][3];
  c[0][4] = c[0][4] - coeff*c[0][3];
  c[1][4] = c[1][4] - coeff*c[1][3];
  c[2][4] = c[2][4] - coeff*c[2][3];
  c[3][4] = c[3][4] - coeff*c[3][3];
  c[4][4] = c[4][4] - coeff*c[4][3];
  r[4]   = r[4]   - coeff*r[3];


  pivot = 1.00/lhs[4][4];
  c[0][4] = c[0][4]*pivot;
  c[1][4] = c[1][4]*pivot;
  c[2][4] = c[2][4]*pivot;
  c[3][4] = c[3][4]*pivot;
  c[4][4] = c[4][4]*pivot;
  r[4]   = r[4]  *pivot;

  coeff = lhs[4][0];
  c[0][0] = c[0][0] - coeff*c[0][4];
  c[1][0] = c[1][0] - coeff*c[1][4];
  c[2][0] = c[2][0] - coeff*c[2][4];
  c[3][0] = c[3][0] - coeff*c[3][4];
  c[4][0] = c[4][0] - coeff*c[4][4];
  r[0]   = r[0]   - coeff*r[4];

  coeff = lhs[4][1];
  c[0][1] = c[0][1] - coeff*c[0][4];
  c[1][1] = c[1][1] - coeff*c[1][4];
  c[2][1] = c[2][1] - coeff*c[2][4];
  c[3][1] = c[3][1] - coeff*c[3][4];
  c[4][1] = c[4][1] - coeff*c[4][4];
  r[1]   = r[1]   - coeff*r[4];

  coeff = lhs[4][2];
  c[0][2] = c[0][2] - coeff*c[0][4];
  c[1][2] = c[1][2] - coeff*c[1][4];
  c[2][2] = c[2][2] - coeff*c[2][4];
  c[3][2] = c[3][2] - coeff*c[3][4];
  c[4][2] = c[4][2] - coeff*c[4][4];
  r[2]   = r[2]   - coeff*r[4];

  coeff = lhs[4][3];
  c[0][3] = c[0][3] - coeff*c[0][4];
  c[1][3] = c[1][3] - coeff*c[1][4];
  c[2][3] = c[2][3] - coeff*c[2][4];
  c[3][3] = c[3][3] - coeff*c[3][4];
  c[4][3] = c[4][3] - coeff*c[4][4];
  r[3]   = r[3]   - coeff*r[4];
}


//#pragma acc routine
void y_binvrhs(double lhs[5][5], double r[5])
{
  double pivot, coeff;

  pivot = 1.00/lhs[0][0];
  lhs[1][0] = lhs[1][0]*pivot;
  lhs[2][0] = lhs[2][0]*pivot;
  lhs[3][0] = lhs[3][0]*pivot;
  lhs[4][0] = lhs[4][0]*pivot;
  r[0]   = r[0]  *pivot;

  coeff = lhs[0][1];
  lhs[1][1]= lhs[1][1] - coeff*lhs[1][0];
  lhs[2][1]= lhs[2][1] - coeff*lhs[2][0];
  lhs[3][1]= lhs[3][1] - coeff*lhs[3][0];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][0];
  r[1]   = r[1]   - coeff*r[0];

  coeff = lhs[0][2];
  lhs[1][2]= lhs[1][2] - coeff*lhs[1][0];
  lhs[2][2]= lhs[2][2] - coeff*lhs[2][0];
  lhs[3][2]= lhs[3][2] - coeff*lhs[3][0];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][0];
  r[2]   = r[2]   - coeff*r[0];

  coeff = lhs[0][3];
  lhs[1][3]= lhs[1][3] - coeff*lhs[1][0];
  lhs[2][3]= lhs[2][3] - coeff*lhs[2][0];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][0];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][0];
  r[3]   = r[3]   - coeff*r[0];

  coeff = lhs[0][4];
  lhs[1][4]= lhs[1][4] - coeff*lhs[1][0];
  lhs[2][4]= lhs[2][4] - coeff*lhs[2][0];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][0];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][0];
  r[4]   = r[4]   - coeff*r[0];


  pivot = 1.00/lhs[1][1];
  lhs[2][1] = lhs[2][1]*pivot;
  lhs[3][1] = lhs[3][1]*pivot;
  lhs[4][1] = lhs[4][1]*pivot;
  r[1]   = r[1]  *pivot;

  coeff = lhs[1][0];
  lhs[2][0]= lhs[2][0] - coeff*lhs[2][1];
  lhs[3][0]= lhs[3][0] - coeff*lhs[3][1];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][1];
  r[0]   = r[0]   - coeff*r[1];

  coeff = lhs[1][2];
  lhs[2][2]= lhs[2][2] - coeff*lhs[2][1];
  lhs[3][2]= lhs[3][2] - coeff*lhs[3][1];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][1];
  r[2]   = r[2]   - coeff*r[1];

  coeff = lhs[1][3];
  lhs[2][3]= lhs[2][3] - coeff*lhs[2][1];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][1];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][1];
  r[3]   = r[3]   - coeff*r[1];

  coeff = lhs[1][4];
  lhs[2][4]= lhs[2][4] - coeff*lhs[2][1];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][1];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][1];
  r[4]   = r[4]   - coeff*r[1];


  pivot = 1.00/lhs[2][2];
  lhs[3][2] = lhs[3][2]*pivot;
  lhs[4][2] = lhs[4][2]*pivot;
  r[2]   = r[2]  *pivot;

  coeff = lhs[2][0];
  lhs[3][0]= lhs[3][0] - coeff*lhs[3][2];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][2];
  r[0]   = r[0]   - coeff*r[2];

  coeff = lhs[2][1];
  lhs[3][1]= lhs[3][1] - coeff*lhs[3][2];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][2];
  r[1]   = r[1]   - coeff*r[2];

  coeff = lhs[2][3];
  lhs[3][3]= lhs[3][3] - coeff*lhs[3][2];
  lhs[4][3]= lhs[4][3] - coeff*lhs[4][2];
  r[3]   = r[3]   - coeff*r[2];

  coeff = lhs[2][4];
  lhs[3][4]= lhs[3][4] - coeff*lhs[3][2];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][2];
  r[4]   = r[4]   - coeff*r[2];


  pivot = 1.00/lhs[3][3];
  lhs[4][3] = lhs[4][3]*pivot;
  r[3]   = r[3]  *pivot;

  coeff = lhs[3][0];
  lhs[4][0]= lhs[4][0] - coeff*lhs[4][3];
  r[0]   = r[0]   - coeff*r[3];

  coeff = lhs[3][1];
  lhs[4][1]= lhs[4][1] - coeff*lhs[4][3];
  r[1]   = r[1]   - coeff*r[3];

  coeff = lhs[3][2];
  lhs[4][2]= lhs[4][2] - coeff*lhs[4][3];
  r[2]   = r[2]   - coeff*r[3];

  coeff = lhs[3][4];
  lhs[4][4]= lhs[4][4] - coeff*lhs[4][3];
  r[4]   = r[4]   - coeff*r[3];


  pivot = 1.00/lhs[4][4];
  r[4]   = r[4]  *pivot;

  coeff = lhs[4][0];
  r[0]   = r[0]   - coeff*r[4];

  coeff = lhs[4][1];
  r[1]   = r[1]   - coeff*r[4];

  coeff = lhs[4][2];
  r[2]   = r[2]   - coeff*r[4];

  coeff = lhs[4][3];
  r[3]   = r[3]   - coeff*r[4];
}

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
  #pragma acc parallel loop private(i,j,k,m,n,tmp1,tmp2,tmp3)
  for (k = 1; k <= grid_points[2]-2; k++) {
    //#pragma acc loop
    for (i = 1; i <= grid_points[0]-2; i++) {
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
      y_lhsinit(lhs[k][i], jsize);
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
      y_binvcrhs( lhs[k][i][0][BB], lhs[k][i][0][CC], rhs[k][0][i] );

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
        y_matvec_sub(lhs[k][i][j][AA], rhs[k][j-1][i], rhs[k][j][i]);

        //-------------------------------------------------------------------
        // B(j) = B(j) - C(j-1)*A(j)
        //-------------------------------------------------------------------
        //#pragma acc routine (matmul_sub) worker
        y_matmul_sub(lhs[k][i][j][AA], lhs[k][i][j-1][CC], lhs[k][i][j][BB]);

        //-------------------------------------------------------------------
        // multiply c[k][j][i] by b_inverse and copy back to c
        // multiply rhs[k][0][i] by b_inverse[k][0][i] and copy to rhs
        //-------------------------------------------------------------------
        //#pragma acc routine (binvcrhs) worker
        y_binvcrhs( lhs[k][i][j][BB], lhs[k][i][j][CC], rhs[k][j][i] );
      }

      //---------------------------------------------------------------------
      // rhs(jsize) = rhs(jsize) - A*rhs(jsize-1)
      //---------------------------------------------------------------------
      //#pragma acc routine (matvec_sub) worker
      y_matvec_sub(lhs[k][i][jsize][AA], rhs[k][jsize-1][i], rhs[k][jsize][i]);

      //---------------------------------------------------------------------
      // B(jsize) = B(jsize) - C(jsize-1)*A(jsize)
      // matmul_sub(AA,i,jsize,k,c,
      // $              CC,i,jsize-1,k,c,BB,i,jsize,k)
      //---------------------------------------------------------------------
      //#pragma acc routine (matmul_sub) worker
      y_matmul_sub(lhs[k][i][jsize][AA], lhs[k][i][jsize-1][CC], lhs[k][i][jsize][BB]);

      //---------------------------------------------------------------------
      // multiply rhs(jsize) by b_inverse(jsize) and copy to rhs
      //---------------------------------------------------------------------
      //#pragma acc routine (binvrhs) worker
      y_binvrhs( lhs[k][i][jsize][BB], rhs[k][jsize][i] );

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
