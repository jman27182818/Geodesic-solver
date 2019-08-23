#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include "math_functions.h"
/* This is the metric file for cartesian desitter.  Use when looking
   at the origin.  If one is not looking at the origin one can use
   metricDSsphere.h as well */


/* To solve geodesics one needs to be specified is g_{ij} and
   dg_{ij}/dx^{k} and g^{ij}.*/




/*defines the  curvature K */
const double K = 0;


// the scale factor function a(t) 
double a(double t){
  double ao=1;
  return ao*sqrt(t);
}

// the derivative of the scale factor function, i.e  a'(t) 
double da(double t){
  double ao=1;
  return 0.5*ao*pow(t,-0.5);
}

//metric function g_{ij}
//TODO: define a struct in order for g to take in a vector or array
double g(double t, double x, double y, double z,int i, int j){
  double r = sqrt(x*x+y*y+z*z);
   double Xc[4]={t,x,y,z};

   if(i > 0 && j > 0){
     return a(t)*a(t)*(kron(i,j)+ K*Xc[i]*Xc[j]*pow(1-K*r*r,-1));
   }
   else if(i==0 && j==0){
     return -1;}
   else {return 0;}
  
}

//inverse metric
double ginv(double t, double x, double y, double z,int i, int j){
  double Xc[4]={t,x,y,z};

   if(i > 0 && j > 0){
     return pow(a(t),-2)*
       ((1-K*Xc[i]*Xc[j])*kron(i,j) -K*(1-kron(i,j))*Xc[i]*Xc[j]);
   }
   else if(i==0 && j==0){
     return -1;}
   else {return 0;}


}


//To solve goedesic equations the fundamental object is dg^{ij}/dx^k

double ginv0(double t, double x, double y, double z, int i, int j){
  double Xc[4]={t,x,y,z};

   if(i > 0 && j > 0){
     return -2*pow(a(t),-3)*da(t)*
       ((1-K*Xc[i]*Xc[j])*kron(i,j) -K*(1-kron(i,j))*Xc[i]*Xc[j]);
   }
   else if(i==0 && j==0){
     return 0;}
   else {return 0;}


}

double ginv1(double t, double x, double y, double z, int i, int j){
  double Xc[4]={t,x,y,z};
  double fudge = (i == j && j==1)? -2*K*pow(a(t),-2): -K*pow(a(t),-2);
  if(i==j && j==0){
    return 0;
  }
  else if(i==1 && j > 0){
    return fudge*Xc[j];
  }
  else if(j==1 && i > 0){
    return fudge*Xc[i];
  }
  else {return 0;}
}

double ginv2(double t, double x, double y, double z, int i, int j){
  double Xc[4]={t,x,y,z};
  double fudge = (i == j && j==2)? -2*K*pow(a(t),-2): -K*pow(a(t),-2);
  if(i==j && j==0){
    return 0;
  }
  else if(i==2 && j > 0){
    return fudge*Xc[j];
  }
  else if(j==2 && i > 0){
    return fudge*Xc[i];
  }
  else {return 0;}

}


double ginv3(double t, double x, double y, double z, int i, int j){

  double Xc[4]={t,x,y,z};
  double fudge = (i == j && j==3)? -2*K*pow(a(t),-2): -K*pow(a(t),-2);
  if(i==j && j==0){
    return 0;
  }
  else if(i==3 && j > 0){
    return fudge*Xc[j];
  }
  else if(j==3 && i > 0){
    return fudge*Xc[i];
  }
  else {return 0;}



}


// ginvderv(t,x,y,z,i,j,k)=dg_{ij}/dx^k
double ginvderv(double t, double x, double y, double z, 
		int i, int j,int k){
  switch(k) {
  case 0:
    return ginv0(t,x,y,z,i,j);
  case 1:
    return ginv1(t,x,y,z,i,j);
  case 2:
    return ginv2(t,x,y,z,i,j);
  case 3:
    return ginv3(t,x,y,z,i,j);
  default:
    printf("incorrect index \n");
    exit(1);
    break;
  }

}












//TODO: calculate acceleration, apparent horizon
