#include "header.h"

//---------------------------------------------------------------------
// compute the right hand side based on exact solution
//---------------------------------------------------------------------
void exact_rhs()
{
  double dtemp[5], xi, eta, zeta, dtpp;
  int m, i, j, k, ip1, im1, jp1, jm1, km1, kp1;

  #pragma acc enter data create(forcing[0:KMAX][0:JMAXP+1][0:IMAXP+1][0:5])

  //---------------------------------------------------------------------
  // initialize                                  
  //---------------------------------------------------------------------
  #pragma acc parallel loop private(i,j,k,m)
  for (k = 0; k <= grid_points[2]-1; k++) {
    #pragma acc loop
    for (j = 0; j <= grid_points[1]-1; j++) {
      #pragma acc loop
      for (i = 0; i <= grid_points[0]-1; i++) {
        for (m = 0; m < 5; m++) {
          forcing[k][j][i][m] = 0.0;
        }
      }
    } 
  }
    
  #pragma acc enter data create(ue[0:PROBLEM_SIZE+1][0:PROBLEM_SIZE+1][0:PROBLEM_SIZE+1][0:5],buf[0:PROBLEM_SIZE+1][0:PROBLEM_SIZE+1][0:PROBLEM_SIZE+1][0:5],\
                                cuf[0:PROBLEM_SIZE+1][0:PROBLEM_SIZE+1][0:PROBLEM_SIZE+1],q[0:PROBLEM_SIZE+1][0:PROBLEM_SIZE+1][0:PROBLEM_SIZE+1])

  //---------------------------------------------------------------------
  // xi-direction flux differences                      
  //---------------------------------------------------------------------
  #pragma acc parallel loop private(i,j,k,m,zeta,eta,xi,dtpp,im1,ip1,dtemp)
  for (k = 1; k <= grid_points[2]-2; k++) {
    #pragma acc loop
    for (j = 1; j <= grid_points[1]-2; j++) {
      #pragma acc loop
      for (i = 0; i <= grid_points[0]-1; i++) {
        zeta = (double)(k) * dnzm1;
        eta = (double)(j) * dnym1;
        xi = (double)(i) * dnxm1;
        #pragma acc routine (exact_solution) worker
        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[k][j][i][m] = dtemp[m];
        }

        dtpp = 1.0 / dtemp[0];

        for (m = 1; m < 5; m++) {
          buf[k][j][i][m] = dtpp * dtemp[m];
        }

        cuf[k][j][i]    = buf[k][j][i][1] * buf[k][j][i][1];
        buf[k][j][i][0] = cuf[k][j][i] + buf[k][j][i][2] * buf[k][j][i][2] + buf[k][j][i][3] * buf[k][j][i][3];
        q[k][j][i] = 0.5*(buf[k][j][i][1]*ue[k][j][i][1] + buf[k][j][i][2]*ue[k][j][i][2] +
                    buf[k][j][i][3]*ue[k][j][i][3]);
      }

      for (i = 1; i <= grid_points[0]-2; i++) {
        im1 = i-1;
        ip1 = i+1;

        forcing[k][j][i][0] = forcing[k][j][i][0] -
          tx2*( ue[k][j][ip1][1]-ue[k][j][im1][1] )+
          dx1tx1*(ue[k][j][ip1][0]-2.0*ue[k][j][i][0]+ue[k][j][im1][0]);

        forcing[k][j][i][1] = forcing[k][j][i][1] - tx2 * (
            (ue[k][j][ip1][1]*buf[k][j][ip1][1]+c2*(ue[k][j][ip1][4]-q[k][j][ip1]))-
            (ue[k][j][im1][1]*buf[k][j][im1][1]+c2*(ue[k][j][im1][4]-q[k][j][im1])))+
          xxcon1*(buf[k][j][ip1][1]-2.0*buf[k][j][i][1]+buf[k][j][im1][1])+
          dx2tx1*( ue[k][j][ip1][1]-2.0* ue[k][j][i][1]+ue[k][j][im1][1]);

        forcing[k][j][i][2] = forcing[k][j][i][2] - tx2 * (
            ue[k][j][ip1][2]*buf[k][j][ip1][1]-ue[k][j][im1][2]*buf[k][j][im1][1])+
          xxcon2*(buf[k][j][ip1][2]-2.0*buf[k][j][i][2]+buf[k][j][im1][2])+
          dx3tx1*( ue[k][j][ip1][2]-2.0*ue[k][j][i][2] +ue[k][j][im1][2]);

        forcing[k][j][i][3] = forcing[k][j][i][3] - tx2*(
            ue[k][j][ip1][3]*buf[k][j][ip1][1]-ue[k][j][im1][3]*buf[k][j][im1][1])+
          xxcon2*(buf[k][j][ip1][3]-2.0*buf[k][j][i][3]+buf[k][j][im1][3])+
          dx4tx1*( ue[k][j][ip1][3]-2.0* ue[k][j][i][3]+ ue[k][j][im1][3]);

        forcing[k][j][i][4] = forcing[k][j][i][4] - tx2*(
            buf[k][j][ip1][1]*(c1*ue[k][j][ip1][4]-c2*q[k][j][ip1])-
            buf[k][j][im1][1]*(c1*ue[k][j][im1][4]-c2*q[k][j][im1]))+
          0.5*xxcon3*(buf[k][j][ip1][0]-2.0*buf[k][j][i][0]+
              buf[k][j][im1][0])+
          xxcon4*(cuf[k][j][ip1]-2.0*cuf[k][j][i]+cuf[k][j][im1])+
          xxcon5*(buf[k][j][ip1][4]-2.0*buf[k][j][i][4]+buf[k][j][im1][4])+
          dx5tx1*( ue[k][j][ip1][4]-2.0* ue[k][j][i][4]+ ue[k][j][im1][4]);
      }

      //---------------------------------------------------------------------
      // Fourth-order dissipation                         
      //---------------------------------------------------------------------
      for (m = 0; m < 5; m++) {
        i = 1;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (5.0*ue[k][j][i][m] - 4.0*ue[k][j][i+1][m] +ue[k][j][i+2][m]);
        i = 2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (-4.0*ue[k][j][i-1][m] + 6.0*ue[k][j][i][m] -
            4.0*ue[k][j][i+1][m] +     ue[k][j][i+2][m]);
      }

      for (i = 3; i <= grid_points[0]-4; i++) {
        for (m = 0; m < 5; m++) {
          forcing[k][j][i][m] = forcing[k][j][i][m] - dssp*
            (ue[k][j][i-2][m] - 4.0*ue[k][j][i-1][m] +
             6.0*ue[k][j][i][m] - 4.0*ue[k][j][i+1][m] + ue[k][j][i+2][m]);
        }
      }

      for (m = 0; m < 5; m++) {
        i = grid_points[0]-3;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[k][j][i-2][m] - 4.0*ue[k][j][i-1][m] +
           6.0*ue[k][j][i][m] - 4.0*ue[k][j][i+1][m]);
        i = grid_points[0]-2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[k][j][i-2][m] - 4.0*ue[k][j][i-1][m] + 5.0*ue[k][j][i][m]);
      }
    }
  }

  //---------------------------------------------------------------------
  // eta-direction flux differences             
  //---------------------------------------------------------------------
  #pragma acc parallel loop private(i,j,k,m,zeta,eta,xi,dtpp,jm1,jp1,dtemp)
  for (k = 1; k <= grid_points[2]-2; k++) {
    #pragma acc loop
    for (i = 1; i <= grid_points[0]-2; i++) {
      #pragma acc loop
      for (j = 0; j <= grid_points[1]-1; j++) {
        zeta = (double)(k) * dnzm1;
        xi = (double)(i) * dnxm1;
        eta = (double)(j) * dnym1;
        #pragma acc routine (exact_solution) worker
        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[k][i][j][m] = dtemp[m];
        }

        dtpp = 1.0/dtemp[0];

        for (m = 1; m < 5; m++) {
          buf[k][i][j][m] = dtpp * dtemp[m];
        }

        cuf[k][i][j]    = buf[k][i][j][2] * buf[k][i][j][2];
        buf[k][i][j][0] = cuf[k][i][j] + buf[k][i][j][1] * buf[k][i][j][1] + buf[k][i][j][3] * buf[k][i][j][3];
        q[k][i][j] = 0.5*(buf[k][i][j][1]*ue[k][i][j][1] + buf[k][i][j][2]*ue[k][i][j][2] +
                    buf[k][i][j][3]*ue[k][i][j][3]);
      }

      for (j = 1; j <= grid_points[1]-2; j++) {
        jm1 = j-1;
        jp1 = j+1;

        forcing[k][j][i][0] = forcing[k][j][i][0] -
          ty2*( ue[k][i][jp1][2]-ue[k][i][jm1][2] )+
          dy1ty1*(ue[k][i][jp1][0]-2.0*ue[k][i][j][0]+ue[k][i][jm1][0]);

        forcing[k][j][i][1] = forcing[k][j][i][1] - ty2*(
            ue[k][i][jp1][1]*buf[k][i][jp1][2]-ue[k][i][jm1][1]*buf[k][i][jm1][2])+
          yycon2*(buf[k][i][jp1][1]-2.0*buf[k][i][j][1]+buf[k][i][jm1][1])+
          dy2ty1*( ue[k][i][jp1][1]-2.0* ue[k][i][j][1]+ ue[k][i][jm1][1]);

        forcing[k][j][i][2] = forcing[k][j][i][2] - ty2*(
            (ue[k][i][jp1][2]*buf[k][i][jp1][2]+c2*(ue[k][i][jp1][4]-q[k][i][jp1]))-
            (ue[k][i][jm1][2]*buf[k][i][jm1][2]+c2*(ue[k][i][jm1][4]-q[k][i][jm1])))+
          yycon1*(buf[k][i][jp1][2]-2.0*buf[k][i][j][2]+buf[k][i][jm1][2])+
          dy3ty1*( ue[k][i][jp1][2]-2.0*ue[k][i][j][2] +ue[k][i][jm1][2]);

        forcing[k][j][i][3] = forcing[k][j][i][3] - ty2*(
            ue[k][i][jp1][3]*buf[k][i][jp1][2]-ue[k][i][jm1][3]*buf[k][i][jm1][2])+
          yycon2*(buf[k][i][jp1][3]-2.0*buf[k][i][j][3]+buf[k][i][jm1][3])+
          dy4ty1*( ue[k][i][jp1][3]-2.0*ue[k][i][j][3]+ ue[k][i][jm1][3]);

        forcing[k][j][i][4] = forcing[k][j][i][4] - ty2*(
            buf[k][i][jp1][2]*(c1*ue[k][i][jp1][4]-c2*q[k][i][jp1])-
            buf[k][i][jm1][2]*(c1*ue[k][i][jm1][4]-c2*q[k][i][jm1]))+
          0.5*yycon3*(buf[k][i][jp1][0]-2.0*buf[k][i][j][0]+
              buf[k][i][jm1][0])+
          yycon4*(cuf[k][i][jp1]-2.0*cuf[k][i][j]+cuf[k][i][jm1])+
          yycon5*(buf[k][i][jp1][4]-2.0*buf[k][i][j][4]+buf[k][i][jm1][4])+
          dy5ty1*(ue[k][i][jp1][4]-2.0*ue[k][i][j][4]+ue[k][i][jm1][4]);
      }

      //---------------------------------------------------------------------
      // Fourth-order dissipation                      
      //---------------------------------------------------------------------
      for (m = 0; m < 5; m++) {
        j = 1;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (5.0*ue[k][i][j][m] - 4.0*ue[k][i][j+1][m] +ue[k][i][j+2][m]);
        j = 2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (-4.0*ue[k][i][j-1][m] + 6.0*ue[k][i][j][m] -
           4.0*ue[k][i][j+1][m] +       ue[k][i][j+2][m]);
      }

      for (j = 3; j <= grid_points[1]-4; j++) {
        for (m = 0; m < 5; m++) {
          forcing[k][j][i][m] = forcing[k][j][i][m] - dssp*
            (ue[k][i][j-2][m] - 4.0*ue[k][i][j-1][m] +
             6.0*ue[k][i][j][m] - 4.0*ue[k][i][j+1][m] + ue[k][i][j+2][m]);
        }
      }

      for (m = 0; m < 5; m++) {
        j = grid_points[1]-3;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[k][i][j-2][m] - 4.0*ue[k][i][j-1][m] +
           6.0*ue[k][i][j][m] - 4.0*ue[k][i][j+1][m]);
        j = grid_points[1]-2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[k][i][j-2][m] - 4.0*ue[k][i][j-1][m] + 5.0*ue[k][i][j][m]);
      }
    }
  }

  //---------------------------------------------------------------------
  // zeta-direction flux differences                     
  //---------------------------------------------------------------------
  #pragma acc parallel loop private(i,j,k,m,zeta,eta,xi,dtpp,km1,kp1,dtemp)
  for (j = 1; j <= grid_points[1]-2; j++) {
    #pragma acc loop
    for (i = 1; i <= grid_points[0]-2; i++) {
      #pragma acc loop
      for (k = 0; k <= grid_points[2]-1; k++) {
        eta = (double)(j) * dnym1;
        xi = (double)(i) * dnxm1;
        zeta = (double)(k) * dnzm1;
        #pragma acc routine (exact_solution) worker
        exact_solution(xi, eta, zeta, dtemp);
        for (m = 0; m < 5; m++) {
          ue[j][i][k][m] = dtemp[m];
        }

        dtpp = 1.0/dtemp[0];

        for (m = 1; m < 5; m++) {
          buf[j][i][k][m] = dtpp * dtemp[m];
        }

        cuf[j][i][k]    = buf[j][i][k][3] * buf[j][i][k][3];
        buf[j][i][k][0] = cuf[j][i][k] + buf[j][i][k][1] * buf[j][i][k][1] + buf[j][i][k][2] * buf[j][i][k][2];
        q[j][i][k] = 0.5*(buf[j][i][k][1]*ue[j][i][k][1] + buf[j][i][k][2]*ue[j][i][k][2] +
                    buf[j][i][k][3]*ue[j][i][k][3]);
      }

      for (k = 1; k <= grid_points[2]-2; k++) {
        km1 = k-1;
        kp1 = k+1;

        forcing[k][j][i][0] = forcing[k][j][i][0] -
          tz2*( ue[j][i][kp1][3]-ue[j][i][km1][3] )+
          dz1tz1*(ue[j][i][kp1][0]-2.0*ue[j][i][k][0]+ue[j][i][km1][0]);

        forcing[k][j][i][1] = forcing[k][j][i][1] - tz2 * (
            ue[j][i][kp1][1]*buf[j][i][kp1][3]-ue[j][i][km1][1]*buf[j][i][km1][3])+
          zzcon2*(buf[j][i][kp1][1]-2.0*buf[j][i][k][1]+buf[j][i][km1][1])+
          dz2tz1*( ue[j][i][kp1][1]-2.0* ue[j][i][k][1]+ ue[j][i][km1][1]);

        forcing[k][j][i][2] = forcing[k][j][i][2] - tz2 * (
            ue[j][i][kp1][2]*buf[j][i][kp1][3]-ue[j][i][km1][2]*buf[j][i][km1][3])+
          zzcon2*(buf[j][i][kp1][2]-2.0*buf[j][i][k][2]+buf[j][i][km1][2])+
          dz3tz1*(ue[j][i][kp1][2]-2.0*ue[j][i][k][2]+ue[j][i][km1][2]);

        forcing[k][j][i][3] = forcing[k][j][i][3] - tz2 * (
            (ue[j][i][kp1][3]*buf[j][i][kp1][3]+c2*(ue[j][i][kp1][4]-q[j][i][kp1]))-
            (ue[j][i][km1][3]*buf[j][i][km1][3]+c2*(ue[j][i][km1][4]-q[j][i][km1])))+
          zzcon1*(buf[j][i][kp1][3]-2.0*buf[j][i][k][3]+buf[j][i][km1][3])+
          dz4tz1*( ue[j][i][kp1][3]-2.0*ue[j][i][k][3] +ue[j][i][km1][3]);

        forcing[k][j][i][4] = forcing[k][j][i][4] - tz2 * (
            buf[j][i][kp1][3]*(c1*ue[j][i][kp1][4]-c2*q[j][i][kp1])-
            buf[j][i][km1][3]*(c1*ue[j][i][km1][4]-c2*q[j][i][km1]))+
          0.5*zzcon3*(buf[j][i][kp1][0]-2.0*buf[j][i][k][0]
              +buf[j][i][km1][0])+
          zzcon4*(cuf[j][i][kp1]-2.0*cuf[j][i][k]+cuf[j][i][km1])+
          zzcon5*(buf[j][i][kp1][4]-2.0*buf[j][i][k][4]+buf[j][i][km1][4])+
          dz5tz1*( ue[j][i][kp1][4]-2.0*ue[j][i][k][4]+ ue[j][i][km1][4]);
      }

      //---------------------------------------------------------------------
      // Fourth-order dissipation                        
      //---------------------------------------------------------------------
      for (m = 0; m < 5; m++) {
        k = 1;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (5.0*ue[j][i][k][m] - 4.0*ue[j][i][k+1][m] +ue[j][i][k+2][m]);
        k = 2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (-4.0*ue[j][i][k-1][m] + 6.0*ue[j][i][k][m] -
           4.0*ue[j][i][k+1][m] +       ue[j][i][k+2][m]);
      }

      for (k = 3; k <= grid_points[2]-4; k++) {
        for (m = 0; m < 5; m++) {
          forcing[k][j][i][m] = forcing[k][j][i][m] - dssp*
            (ue[j][i][k-2][m] - 4.0*ue[j][i][k-1][m] +
             6.0*ue[j][i][k][m] - 4.0*ue[j][i][k+1][m] + ue[j][i][k+2][m]);
        }
      }

      for (m = 0; m < 5; m++) {
        k = grid_points[2]-3;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[j][i][k-2][m] - 4.0*ue[j][i][k-1][m] +
           6.0*ue[j][i][k][m] - 4.0*ue[j][i][k+1][m]);
        k = grid_points[2]-2;
        forcing[k][j][i][m] = forcing[k][j][i][m] - dssp *
          (ue[j][i][k-2][m] - 4.0*ue[j][i][k-1][m] + 5.0*ue[j][i][k][m]);
      }

    }
  }
  #pragma acc exit data delete(ue,buf,cuf,q)
  //---------------------------------------------------------------------
  // now change the sign of the forcing function, 
  //---------------------------------------------------------------------
  #pragma acc parallel loop private(i,j,k,m)
  for (k = 1; k <= grid_points[2]-2; k++) {
    #pragma acc loop
    for (j = 1; j <= grid_points[1]-2; j++) {
      #pragma acc loop
      for (i = 1; i <= grid_points[0]-2; i++) {
        for (m = 0; m < 5; m++) {
          forcing[k][j][i][m] = -1.0 * forcing[k][j][i][m];
        }
      }
    }
  }
}
