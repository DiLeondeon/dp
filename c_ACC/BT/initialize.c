#include "header.h"

//---------------------------------------------------------------------
// This subroutine initializes the field variable u using 
// tri-linear transfinite interpolation of the boundary values     
//---------------------------------------------------------------------
/*void exact_solution_1(double xi, double eta, double zeta, double dtemp[5], double ce[5][13])
{
  int m;

  for (m = 0; m < 5; m++) {
    dtemp[m] =  ce[m][0] +
      xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
      eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
      zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
      zeta*ce[m][12])));
  }
}*/

void initialize()
{
  int i, j, k, m, ix, iy, iz;
  double xi, eta, zeta, Pface[2][3][5], Pxi, Peta, Pzeta, temp[5];

  //---------------------------------------------------------------------
  // Later (in compute_rhs) we compute 1/u for every element. A few of 
  // the corner elements are not used, but it convenient (and faster) 
  // to compute the whole thing with a simple loop. Make sure those 
  // values are nonzero by initializing the whole thing here. 
  //---------------------------------------------------------------------
  #pragma acc parallel loop private(i,j,k,m) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      //#pragma acc loop
      for (i = 0; i <= grid_points[0]-1; i++) {
        for (m = 0; m < 5; m++) {
          u[k][j][i][m] = 1.0;
        }
      }
    }
  }

  //---------------------------------------------------------------------
  // first store the "interpolated" values everywhere on the grid    
  //---------------------------------------------------------------------
  #pragma acc parallel loop private(i,j,k,zeta,eta,xi,ix,iy,iz,m,Pface,Pxi,Peta,Pzeta,m) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      //#pragma acc loop
      for (i = 0; i <= grid_points[0]-1; i++) {
        zeta = (double)(k) * dnzm1;
        eta = (double)(j) * dnym1;
        xi = (double)(i) * dnxm1;

        for (ix = 0; ix < 2; ix++) {
          //#pragma acc routine (exact_solution) worker
          //exact_solution_1((double)ix, eta, zeta, &Pface[ix][0][0], ce);
          //void exact_solution_1(double xi, double eta, double zeta, double dtemp[5], double ce[5][13])
          int m;

          for (m = 0; m < 5; m++) {
            Pface[ix][0][m] =  ce[m][0] +
            (double)ix*(ce[m][1] + (double)ix*(ce[m][4] + (double)ix*(ce[m][7] + (double)ix*ce[m][10]))) +
            eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
            zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
            zeta*ce[m][12])));
          }
        }

        for (iy = 0; iy < 2; iy++) {
          //#pragma acc routine (exact_solution) worker
          //exact_solution_1(xi, (double)iy , zeta, &Pface[iy][1][0], ce);
          int m;

          for (m = 0; m < 5; m++) {
            Pface[iy][1][m] =  ce[m][0] +
            xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
            (double)iy*(ce[m][2] + (double)iy*(ce[m][5] + (double)iy*(ce[m][8] + (double)iy*ce[m][11])))+
            zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
            zeta*ce[m][12])));
          }
        }

        for (iz = 0; iz < 2; iz++) {
          //#pragma acc routine (exact_solution) worker
          //exact_solution_1(xi, eta, (double)iz, &Pface[iz][2][0], ce);
          int m;

        for (m = 0; m < 5; m++) {
            Pface[iz][2][m] =  ce[m][0] +
            xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
            eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
            (double)iz*(ce[m][3] + (double)iz*(ce[m][6] + (double)iz*(ce[m][9] + 
            (double)iz*ce[m][12])));
          }
        }

        for (m = 0; m < 5; m++) {
          Pxi   = xi   * Pface[1][0][m] + (1.0-xi)   * Pface[0][0][m];
          Peta  = eta  * Pface[1][1][m] + (1.0-eta)  * Pface[0][1][m];
          Pzeta = zeta * Pface[1][2][m] + (1.0-zeta) * Pface[0][2][m];

          u[k][j][i][m] = Pxi + Peta + Pzeta - 
                          Pxi*Peta - Pxi*Pzeta - Peta*Pzeta + 
                          Pxi*Peta*Pzeta;
        }
      }
    }
  }

  //---------------------------------------------------------------------
  // now store the exact values on the boundaries        
  //---------------------------------------------------------------------

  //---------------------------------------------------------------------
  // west face                                                  
  //---------------------------------------------------------------------
  i = 0;
  xi = 0.0;
  #pragma acc parallel loop private(j,k,m,zeta,eta,temp) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      zeta = (double)(k) * dnzm1;
      eta = (double)(j) * dnym1;
      //#pragma acc routine (exact_solution) worker
      //exact_solution_1(xi, eta, zeta, temp, ce);
      //void exact_solution_1(double xi, double eta, double zeta, double dtemp[5], double ce[5][13])
      int m;

      for (m = 0; m < 5; m++) {
        temp[m] =  ce[m][0] +
        xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
        eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
        zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
        zeta*ce[m][12])));
      }
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // east face                                                      
  //---------------------------------------------------------------------
  i = grid_points[0]-1;
  xi = 1.0;
  #pragma acc parallel loop private(j,k,m,zeta,eta,temp) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      zeta = (double)(k) * dnzm1;
      eta = (double)(j) * dnym1;
      //#pragma acc routine (exact_solution) worker
      //exact_solution_1(xi, eta, zeta, temp, ce);
      int m;

      for (m = 0; m < 5; m++) {
        temp[m] =  ce[m][0] +
        xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
        eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
        zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
        zeta*ce[m][12])));
      }
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // south face                                                 
  //---------------------------------------------------------------------
  j = 0;
  eta = 0.0;
  #pragma acc parallel loop private(i,k,m,zeta,xi,temp) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (i = 0; i <= grid_points[0]-1; i++) {
      zeta = (double)(k) * dnzm1;
      xi = (double)(i) * dnxm1;
      //#pragma acc routine (exact_solution) worker
      //exact_solution_1(xi, eta, zeta, temp, ce);
      int m;

      for (m = 0; m < 5; m++) {
        temp[m] =  ce[m][0] +
        xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
        eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
        zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
        zeta*ce[m][12])));
      }
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // north face                                    
  //---------------------------------------------------------------------
  j = grid_points[1]-1;
  eta = 1.0;
  #pragma acc parallel loop private(i,k,m,zeta,xi,temp) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (i = 0; i <= grid_points[0]-1; i++) {
      zeta = (double)(k) * dnzm1;
      xi = (double)(i) * dnxm1;
      //#pragma acc routine (exact_solution) worker
      //exact_solution_1(xi, eta, zeta, temp, ce);
      int m;

      for (m = 0; m < 5; m++) {
        temp[m] =  ce[m][0] +
        xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
        eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
        zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
        zeta*ce[m][12])));
      }
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // bottom face                                       
  //---------------------------------------------------------------------
  k = 0;
  zeta = 0.0;
  #pragma acc parallel loop private(i,j,m,eta,xi,temp) 
  for (j = 0; j <= grid_points[1]-1; j++) {
    //#pragma acc loop
    for (i =0; i <= grid_points[0]-1; i++) {
      eta = (double)(j) * dnym1;
      xi = (double)(i) *dnxm1;
      //#pragma acc routine (exact_solution) worker
      //exact_solution_1(xi, eta, zeta, temp, ce);
      int m;

      for (m = 0; m < 5; m++) {
        temp[m] =  ce[m][0] +
        xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
        eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
        zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
        zeta*ce[m][12])));
      }
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }

  //---------------------------------------------------------------------
  // top face     
  //---------------------------------------------------------------------
  k = grid_points[2]-1;
  zeta = 1.0;
  #pragma acc parallel loop private(i,j,m,eta,xi,temp) 
  for (j = 0; j <= grid_points[1]-1; j++) {
    //#pragma acc loop
    for (i = 0; i <= grid_points[0]-1; i++) {
      eta = (double)(j) * dnym1;
      xi = (double)(i) * dnxm1;
      //#pragma acc routine (exact_solution) worker
      //exact_solution_1(xi, eta, zeta, temp, ce);
      int m;

      for (m = 0; m < 5; m++) {
        temp[m] =  ce[m][0] +
        xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
        eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
        zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
        zeta*ce[m][12])));
      }
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
}


/*void lhsinit(double lhs[][3][5][5], int ni)
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
}*/
/*#pragma acc routine
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
}*/
