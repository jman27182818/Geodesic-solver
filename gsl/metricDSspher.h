#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include "math_functions.h"


/* To solve geodesics one needs to be specified is g_{ij} and
   dg_{ij}/dx^{k} and g^{ij}.*/




/*defines the horizon size (this is specific to this metric ie static
  Desitter) */
const double a = 100;


//metric function g_{ij}
//TODO: define a struct in order for g to take in a vector or array
double g(double t, double r, double th, double phi,int i, int j){
  double gtt;
  gtt = 1-pow(r/a,2);
  if(abs(i-j) > 0){
    return 0;}
  else if(i==0){
    return -gtt;
}
  else if(i==1){
    return 1/gtt;
}
  else if(i==2){
    return r*r;
}
  else if(i==3){
    return r*r*sin(th)*sin(th);
}
  else{
    printf("invalid index");
    return 0;
}
}

//inverse metric
double ginv(double t, double r, double th, double phi,int i, int j){
  double gtt;
  gtt = 1-pow(r/a,2);
  if(abs(i-j) > 0){
    return 0;}
  else if(i==0){
    return -(1/gtt);
}
  else if(i==1){
    return gtt;
}
  else if(i==2){
    return pow(r*r,-1);
}
  else if(i==3){
    return pow(r*r*sin(th)*sin(th),-1);
}
  else{
    printf("invalid index \n");
    return 0;
}
}


//To solve goedesic equations the fundamental object is dg^{ij}/dx^k

double ginv0(double t, double r, double th, double phi, int i, int j){
  return 0;
}

double ginv1(double t, double r, double th, double phi, int i, int j){
  
  if(j != i){return 0;}
  else if(j==0){return (-2)*r*a*a*pow(a*a-r*r,-2);}
  else if(j==1){return (-2)*r*pow(a,-2);}
  else if(j==2){return (-2)*pow(r,-3);}
  else if(j==3){return (-2)*pow(sin(th),-2)*pow(r,-3);}
  else{
    printf("invalid index\n");
    return 0;}
}

double ginv2(double t, double r, double th, double phi, int i, int j){
  if(i == 3 && j == 3){
    return (-2)*pow(tan(th),-1)*pow(sin(th),-2)*pow(r,-2);}
  else{return 0;}
}


double ginv3(double t, double r, double th, double phi, int i, int j){
  return 0;
}


// ginvderv(t,r,th,phi,i,j,k)=dg_{ij}/dx^k
double ginvderv(double t, double r, double th, double phi, 
		int i, int j,int k){
  switch(k) {
  case 0:
    return ginv0(t,r,th,phi,i,j);
  case 1:
    return ginv1(t,r,th,phi,i,j);
  case 2:
    return ginv2(t,r,th,phi,i,j);
  case 3:
    return ginv3(t,r,th,phi,i,j);
  default:
    printf("incorrect index \n");
    exit(1);
    break;
  }

}






//TODO: calculate acceleration, apparent horizon
