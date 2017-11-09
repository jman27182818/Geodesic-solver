#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>
#include "metricSchwar.h"

//The initial parameters which need to be specified are x_0^{\mu}, t_0^{\mu}, and e^{\mu}_{a}(0)

//define the dimensions
#define DIM 4

//tells whether the metric is diagonal or not 0 = True, 1= Fals
#define DIAGONAL 0

//initial point (in the global coordinates)
const double x0[4] = {0,1e-12,3.14159/2,0};


//initial velocity (make sure it is normalized)
const double v0[4] = {1.1547,0.5*1.1547,0,0};
//const double v0[4] = {1,0,0,0};



double Christoffel(double t, double x, double y, double z, int i,
		   int k, int l){
  double term1, term2, CH;
  
  term1=0;
  term2=0;
  int mm=0;
  int vv=0;
  int ss=0;
  for(mm=0; mm< 4; mm++){
    term1 += g(t,x,y,z,mm,k)*ginvderv(t,x,y,z,mm,i,l)
      + g(t,x,y,z,mm,l)*ginvderv(t,x,y,z,mm,i,k);}
  term1 =(-1)*term1;
  
  for(mm=0; mm< 4; mm++){
    for(vv=0; vv< 4; vv++){
      for(ss=0; ss< 4; ss++){
	term2 += ginv(t,x,y,z, i,mm)
	  *g(t,x,y,z,vv,k)*g(t,x,y,z,ss,l)
	  *ginvderv(t,x,y,z,ss,vv,mm);}}}
  CH = term1+term2;      

  return 0.5*CH;

}


int main(){

  double T = 0;
  double X = 800;
  double Y = (3.141592653)/2;
  double Z = 0;
  //metric
  int jj;
  for(jj=0; jj < 4;jj++){
    printf("%g \t %g \t %g \t %g \n", 
	   g(T,X,Y,Z,jj,0),
	   g(T,X,Y,Z,jj,1),
	   g(T,X,Y,Z,jj,2),
	   g(T,X,Y,Z,jj,3));
      }

  printf("\n \n");

  for(jj=0; jj < 4;jj++){
    printf("%g \t %g \t %g \t %g \n", 
	   ginv(T,X,Y,Z,jj,0),
	   ginv(T,X,Y,Z,jj,1),
	   ginv(T,X,Y,Z,jj,2),
	   ginv(T,X,Y,Z,jj,3));
      }


  printf("\n \n");

  int kk = 1;
  for(jj=0; jj < 4;jj++){
    printf("%g \t %g \t %g \t %g \n", 
	   ginvderv(T,X,Y,Z,jj,0,kk),
	   ginvderv(T,X,Y,Z,jj,1,kk),
	   ginvderv(T,X,Y,Z,jj,2,kk),
	   ginvderv(T,X,Y,Z,jj,3,kk));
      }


  printf("\n \n");

  int cc = 0;
  for(jj=0; jj < 4;jj++){
    printf("%g \t %g \t %g \t %g \n", 
	   Christoffel(T,X,Y,Z,cc,jj,0),
	   Christoffel(T,X,Y,Z,cc,jj,1),
	   Christoffel(T,X,Y,Z,cc,jj,2),
	   Christoffel(T,X,Y,Z,cc,jj,3));
      }

  //  printf("%g \n",Christoffel(3,1,2,3,3,3,0));
  return 0;
}



		
