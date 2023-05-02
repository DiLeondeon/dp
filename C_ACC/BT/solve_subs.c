//#pragma acc routine
/*void matvec_sub(double ablock[5][5], double avec[5], double bvec[5])
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
void matmul_sub(double ablock[5][5], double bblock[5][5], double cblock[5][5])
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
void binvcrhs(double lhs[5][5], double c[5][5], double r[5])
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
void binvrhs(double lhs[5][5], double r[5])
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
*/
#pragma acc routine
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
#pragma acc routine
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
#pragma acc routine
void binvcrhs(int k1, int j1, int i1, int first, int k2, int j2, int i2, int second, int k3, int j3, int i3, double lhs[PROBLEM_SIZE+1][PROBLEM_SIZE+1][PROBLEM_SIZE+1][3][5][5], double rhs[KMAX][JMAXP+1][IMAXP+1][5])//(double lhs[5][5], double c[5][5], double r[5]) //x_binvcrhs( lhs[k][j][i][BB], lhs[k][j][i][CC][, rhs[k][j][i][ );
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
#pragma acc routine
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