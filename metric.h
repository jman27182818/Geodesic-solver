#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>

/* This is the metric file for cartesian desitter.  Use when looking
   at the origin.  If one is not looking at the origin one can use
   metricDSsphere.h as well */
const double pi = 3.14159265358979323846;

/* To solve geodesics one needs to be specified is g_{ij} and
   dg_{ij}/dx^{k} and g^{ij}.*/




/*defines the horizon size (this is specific to this metric ie static
  Desitter) */
const double a = 100;


//metric function g_{ij}
//TODO: define a struct in order for g to take in a vector or array
double g(double t, double x, double y, double z,int i, int j){
  double r = sqrt(x*x+y*y+z*z);
  double rdif= r*r-a*a;
  double Xc[4]={t,x,y,z};
  if( i==j){
    switch(j) {
    case 0:
      return -1+pow(r/a,2);
    case 1:
      return 1- x*x*pow(rdif,-1);
    case 2:
      return 1- y*y*pow(rdif,-1);
    case 3:
      return 1- z*z*pow(rdif,-1);
    }
  }
  else if(i > 0 && j >0){
    return (-1)*Xc[i]*Xc[j]*pow(rdif,-1);
      }
  else {return 0;}
  
}

//inverse metric
double ginv(double t, double x, double y, double z,int i, int j){
  double r = sqrt(x*x+y*y+z*z);
  double Xc[4]={t,x,y,z};
  if( i==j){
    switch(j) {
    case 0:
      return pow(-1+pow(r/a,2),-1);
    case 1:
      return 1- x*x*pow(a,-2);
    case 2:
      return 1- y*y*pow(a,-2);
    case 3:
      return 1- z*z*pow(a,-2);
    }
  }
  else if(i > 0 && j >0){
    return (-1)*Xc[i]*Xc[j]*pow(a,-2);
      }
  else {return 0;}

}


//To solve goedesic equations the fundamental object is dg^{ij}/dx^k

double ginv0(double t, double x, double y, double z, int i, int j){
  return 0;
}

double ginv1(double t, double x, double y, double z, int i, int j){
  double r = sqrt(x*x+y*y+z*z);
  double Xc[4]={t,x,y,z};
  double fudge = (i == j && j==1)? 2 : 1;
  if(i==j && j==0){
    return -2*x*a*a*pow(r*r-a*a,-2);
  }
  else if(i==1 && j > 0){
    return (-1)*fudge*Xc[j]*pow(a,-2);
  }
  else if(j==1 && i > 0){
    return (-1)*fudge*Xc[i]*pow(a,-2);
  }
  else {return 0;}
}

double ginv2(double t, double x, double y, double z, int i, int j){
  double r = sqrt(x*x+y*y+z*z);
  double Xc[4]={t,x,y,z};
  double fudge = (i == j && j==2)? 2 : 1;
  if(i==j && j==0){
    return -2*y*a*a*pow(r*r-a*a,-2);
  }
  else if(i==2 && j > 0){
    return -fudge*Xc[j]*pow(a,-2);
  }
  else if(j==2 && i > 0){
    return -fudge*Xc[i]*pow(a,-2);
  }
  else {return 0;}
}


double ginv3(double t, double x, double y, double z, int i, int j){

  double r = sqrt(x*x+y*y+z*z);
  double Xc[4]={t,x,y,z};
  double fudge = (i == j && j==3)? 2 : 1;
  if(i==j && j==0){
    return -2*z*a*a*pow(r*r-a*a,-2);
  }
  else if(i==3 && j > 0){
    return -fudge*Xc[j]*pow(a,-2);
  }
  else if(j==3 && i > 0){
    return -fudge*Xc[i]*pow(a,-2);
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


// Christoffel \Gamma^{i}_{jk} (automatically calculated)

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





//TODO: calculate acceleration, apparent horizon
