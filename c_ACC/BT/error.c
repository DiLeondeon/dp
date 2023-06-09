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
  double rms_local_0 = 0.0, rms_local_1 = 0.0, rms_local_2 = 0.0, rms_local_3 = 0.0, rms_local_4 = 0.0;


  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  }

  //#pragma acc parallel private(i,j,k,m,zeta,eta,xi,add,u_exact,m) 
  {
  //#pragma acc loop reduction(+:rms_local_0,rms_local_1,rms_local_2,rms_local_3,rms_local_4)
  for (k = 0; k <= grid_points[2]-1; k++) {
    //#pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      //#pragma acc loop
      for (i = 0; i <= grid_points[0]-1; i++) {
        zeta = (double)(k) * dnzm1;
        eta = (double)(j) * dnym1;
        xi = (double)(i) * dnxm1;
        //#pragma acc routine (exact_solution) worker
        //exact_solution(xi, eta, zeta, u_exact, ce);
        //void exact_solution(double xi, double eta, double zeta, double dtemp[5], double ce[5][13])
        int m;

        for (m = 0; m < 5; m++) {
          u_exact[m] =  ce[m][0] +
          xi*(ce[m][1] + xi*(ce[m][4] + xi*(ce[m][7] + xi*ce[m][10]))) +
          eta*(ce[m][2] + eta*(ce[m][5] + eta*(ce[m][8] + eta*ce[m][11])))+
          zeta*(ce[m][3] + zeta*(ce[m][6] + zeta*(ce[m][9] + 
          zeta*ce[m][12])));
        }

        add = u[k][j][i][0]-u_exact[0];
        rms_local_0 = rms_local_0 + add*add;
        add = u[k][j][i][1]-u_exact[1];
        rms_local_1 = rms_local_1 + add*add;
        add = u[k][j][i][2]-u_exact[2];
        rms_local_2 = rms_local_2 + add*add;
        add = u[k][j][i][3]-u_exact[3];
        rms_local_3 = rms_local_3 + add*add;
        add = u[k][j][i][4]-u_exact[4];
        rms_local_4 = rms_local_4 + add*add;
      }
    }
  }
  }
  rms[0] += rms_local_0;
  rms[1] += rms_local_1;
  rms[2] += rms_local_2;
  rms[3] += rms_local_3;
  rms[4] += rms_local_4;


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
  double rms_local_0 = 0.0, rms_local_1 = 0.0, rms_local_2 = 0.0, rms_local_3 = 0.0, rms_local_4 = 0.0;

  for (m = 0; m < 5; m++) {
    rms[m] = 0.0;
  } 

  //#pragma acc parallel private(i,j,k,add)
  {
  //#pragma acc loop reduction(+:rms_local_0,rms_local_1,rms_local_2,rms_local_3,rms_local_4)
  for (k = 1; k <= grid_points[2]-2; k++) {
    //#pragma acc loop
    for (j = 1; j <= grid_points[1]-2; j++) {
      //#pragma acc loop
      for (i = 1; i <= grid_points[0]-2; i++) {
          add = rhs[k][j][i][0];
          rms_local_0 = rms_local_0 + add*add;
          add = rhs[k][j][i][1];
          rms_local_1 = rms_local_1 + add*add;
          add = rhs[k][j][i][2];
          rms_local_2 = rms_local_2 + add*add;
          add = rhs[k][j][i][3];
          rms_local_3 = rms_local_3 + add*add;
          add = rhs[k][j][i][4];
          rms_local_4 = rms_local_4 + add*add;
      } 
    } 
  } 
  }
  rms[0] += rms_local_0;
  rms[1] += rms_local_1;
  rms[2] += rms_local_2;
  rms[3] += rms_local_3;
  rms[4] += rms_local_4;

  for (m = 0; m < 5; m++) {
    for (d = 0; d < 3; d++) {
      rms[m] = rms[m] / (double)(grid_points[d]-2);
    } 
    rms[m] = sqrt(rms[m]);
  } 
}
