#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>

//defines the mass of the black hole (this is specific to this metric)
#define M = 100;


//define the dimensions
#define DIM 4

//tells whether the metric is diagonal or not 0 = True, 1= Fals
#define DIAGONAL 0


//initial point (in the global coordinates)
double x0[4]={0,800,pi/2,0};

//initial velocity
double v0[3]={0,0,0};

//metric function g_{ij}
double g(double t, double r, double th, double phi,int i, int j){
  if(i != j){
    return 0;}
  else if(i=0){
    return (-1)*(1-(2*M)/r);}
  else if(i=1){
    return pow(1-(2*M)/r,-1);}
  else if(i=2){
    return r*r;}
  else if(i=3){
    return r*r*sin(th)*sin(th);}
  else{
    printf("invalid index")}
}

// Chirstoffel \Gamma^{i}_{jk}

double Christoffel(double t, double r, double th, double phi, int i,
		   int j, int k){




}


//TODO: calculate acceleration, apparent horizon
