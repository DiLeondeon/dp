#include <math.h>
#include "header.h"


//---------------------------------------------------------------------
// this function computes the norm of the difference between the
// computed solution and the exact solution
//---------------------------------------------------------------------
void error_norm(double rms[5])
{
  int i, j, k, m, d;
  double xi, eta, zeta, u_exact[5], add;
  double rms_local[5];

  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }
  //#pragma acc enter data copyin(rms[0:5])
  #pragma acc kernels//#pragma acc parallel private(i,j,k,m,zeta,eta,xi,add,u_exact,rms_local)
  {
  for (m = 0; m < 5; m++) {
    rms_local[m] = 0.0;
  }
  //#pragma acc loop 
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      //#pragma acc loop
      for (i = 0; i <= grid_points[0]-1; i++) {
        zeta = (double)(k) * dnzm1;
        eta = (double)(j) * dnym1;
        xi = (double)(i) * dnxm1;
        #pragma acc routine (exact_solution) worker
        exact_solution(xi, eta, zeta, u_exact);

        for (m = 0; m < 5; m++) {
          add = u[k][j][i][m]-u_exact[m];
          rms_local[m] = rms_local[m] + add*add;
        }
      }
    }
  }
  for (m = 0; m < 5; m++) {
    //#pragma acc atomic
    rms[m] += rms_local[m];
  }
  } //end parallel

  for (m = 0; m < 5; m++) {
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    }
    rms[m] = sqrt(rms[m]);
  }
}


void rhs_norm(double rms[5])
{
  int i, j, k, d, m;
  double add;
  double rms_local[5];

  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  } 

  #pragma acc kernels//#pragma acc parallel private(i,j,k,m,add,rms_local)
  {
  for (m = 0; m < 5; m++) {
    rms_local[m] = 0.0;
  }
  //#pragma acc loop
  for (k = 1; k <= grid_points[2]-2; k++) {
    //#pragma acc loop
    for (j = 1; j <= grid_points[1]-2; j++) {
      //#pragma acc loop
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (m = 0; m < 5; m++) {
          add = rhs[k][j][i][m];
          rms_local[m] = rms_local[m] + add*add;
        } 
      } 
    } 
  } 
  for (m = 0; m < 5; m++) {
    //#pragma acc atomic
    rms[m] += rms_local[m];
  }
  } //end parallel

  for (m = 0; m < 5; m++) {
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    } 
    rms[m] = sqrt(rms[m]);
  } 
  //#pragma acc exit data copyout(rms[0:5])
}
