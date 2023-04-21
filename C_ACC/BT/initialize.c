#include "header.h"

//---------------------------------------------------------------------
// This subroutine initializes the field variable u using 
// tri-linear transfinite interpolation of the boundary values     
//---------------------------------------------------------------------
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
  #pragma acc kernels//#pragma acc parallel loop private(i,j,k,m) 
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
  #pragma acc kernels//#pragma acc parallel loop private(i,j,k,zeta,eta,xi,ix,iy,iz,m,Pface,Pxi,Peta,Pzeta) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      //#pragma acc loop
      for (i = 0; i <= grid_points[0]-1; i++) {
        zeta = (double)(k) * dnzm1;
        eta = (double)(j) * dnym1;
        xi = (double)(i) * dnxm1;

        for (ix = 0; ix < 2; ix++) {
          //#pragma acc routine (exact_solution) (exact_solution) worker
          exact_solution((double)ix, eta, zeta, &Pface[ix][0][0]);
        }

        for (iy = 0; iy < 2; iy++) {
          //#pragma acc routine (exact_solution) (exact_solution) worker
          exact_solution(xi, (double)iy , zeta, &Pface[iy][1][0]);
        }

        for (iz = 0; iz < 2; iz++) {
          //#pragma acc routine (exact_solution) worker
          exact_solution(xi, eta, (double)iz, &Pface[iz][2][0]);
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
  #pragma acc kernels//#pragma acc parallel loop private(j,k,m,zeta,eta,temp) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      zeta = (double)(k) * dnzm1;
      eta = (double)(j) * dnym1;
      //#pragma acc routine (exact_solution) worker
      exact_solution(xi, eta, zeta, temp);
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
  #pragma acc kernels//#pragma acc parallel loop private(j,k,m,zeta,eta,temp) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      zeta = (double)(k) * dnzm1;
      eta = (double)(j) * dnym1;
      //#pragma acc routine (exact_solution) worker
      exact_solution(xi, eta, zeta, temp);
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
  #pragma acc kernels//#pragma acc parallel loop private(i,k,m,zeta,xi,temp) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (i = 0; i <= grid_points[0]-1; i++) {
      zeta = (double)(k) * dnzm1;
      xi = (double)(i) * dnxm1;
      //#pragma acc routine (exact_solution) worker
      exact_solution(xi, eta, zeta, temp);
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
  #pragma acc kernels//#pragma acc parallel loop private(i,k,m,zeta,xi,temp) 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (i = 0; i <= grid_points[0]-1; i++) {
      zeta = (double)(k) * dnzm1;
      xi = (double)(i) * dnxm1;
      //#pragma acc routine (exact_solution) worker
      exact_solution(xi, eta, zeta, temp);
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
  #pragma acc kernels//#pragma acc parallel loop private(i,j,m,eta,xi,temp) 
  for (j = 0; j <= grid_points[1]-1; j++) {
    //#pragma acc loop
    for (i =0; i <= grid_points[0]-1; i++) {
      eta = (double)(j) * dnym1;
      xi = (double)(i) *dnxm1;
      //#pragma acc routine (exact_solution) worker
      exact_solution(xi, eta, zeta, temp);
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
  #pragma acc kernels//#pragma acc parallel loop private(i,j,m,eta,xi,temp) 
  for (j = 0; j <= grid_points[1]-1; j++) {
    //#pragma acc loop
    for (i = 0; i <= grid_points[0]-1; i++) {
      eta = (double)(j) * dnym1;
      xi = (double)(i) * dnxm1;
      //#pragma acc routine (exact_solution) worker
      exact_solution(xi, eta, zeta, temp);
      for (m = 0; m < 5; m++) {
        u[k][j][i][m] = temp[m];
      }
    }
  }
}


void lhsinit(double lhs[][3][5][5], int ni)
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
