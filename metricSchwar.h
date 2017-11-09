#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>

//defines the mass of the black hole (this is specific to this metric)
const double M = 100;


//define the dimensions
#define DIM 4

//tells whether the metric is diagonal or not 0 = True, 1= Fals
#define DIAGONAL 0



//metric function g_{ij}
double g(double t, double r, double th, double phi,int i, int j){
  double RR = 1-2*M*pow(r,-1);
  if(i != j){
    return 0;}
  else if(i==0){
    return -RR;}
  else if(i==1){
    return pow(RR,-1);}
  else if(i==2){
    return r*r;}
  else if(i==3){
    return r*r*sin(th)*sin(th);}
  else{
    printf("invalid index");
    return 0;}
}


double ginv(double t, double r, double th, double phi,int i, int j){
  if(i != j){
    return 0;}
  else if(i==0){
    return pow((-1)*(1-(2*M)/r),-1);}
  else if(i==1){
    return pow(1-(2*M)/r,1);}
  else if(i==2){
    return pow(r*r,-1);}
  else if(i==3){
    return pow(r*r*sin(th)*sin(th),-1);}
  else{
    printf("invalid index");
    return 0;}
}

//To solve goedesic equations the fundamental object is dg^{ij}/dx^k

//ginv0 = dg/dt
double ginv0(double t, double r, double th, double phi,int i, int j){
  return 0;
}


double ginv1(double t, double r, double th, double phi,int i, int j){
  if(i != j){
    return 0;}
  else if(i==0){
    return 2*M*pow(r-2*M,-2);}
  else if(i==1){
    return 2*M*pow(r,-2);}
  else if(i==2){
    return -2*pow(r,-3);}
  else if(i==3){
    return -2*pow(r,-3)*pow(sin(th),-2);}
  else{
    printf("invalid index");
    return 0;}
}



double ginv2(double t, double r, double th, double phi,int i, int j){
  if(i==j && i==3){
    return -2*pow(sin(th)*r,-3)*cos(th);}
  else {return 0;}
}


double ginv3(double t, double r, double th, double phi,int i, int j){
  return 0;
}



// ginvderv(t,x,y,z,i,j,k)=dg_{ij}/dx^k
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
