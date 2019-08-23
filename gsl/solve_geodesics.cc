#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_linalg.h>

//change metric file as desired.
#include "metricSchwar.h"


//The initial parameters which need to be specified are x_0^{\mu}, t_0^{\mu}, and e^{\mu}_{a}(0)

//define the dimensions
#define DIM 4

//tells whether the metric is diagonal or not 0 = True, 1= Fals
#define DIAGONAL 0

//initial point (in the global coordinates)
const double x0[4] = {0,400,3.14/2,0};
//note try not to use theta=0

//initial velocity v^{\mu}_0 (make sure it is approximately normalized)
const double rdt = sqrt(2*M*pow(x0[1],-1))+0.1;

const double v0[4] = {sqrt(1+rdt*rdt-2*M*pow(x0[1],-1))*pow(1-2*M*pow(x0[1],-1),-1),rdt,0,0};




int geodesics (double tau, const double y[], double f[],
	  void *params){
  double *kparams = (double*) params;
  double sign = kparams[0];
  int uu=0;
  for(uu=0; uu < 8; uu++){
    f[uu] = 0;
}

  //geodesic as hamiltonian flow setup
  int jj=0;
  int ii=0;
  int kk=0;
    for(ii=0; ii < 4; ii++){
      for(jj=0; jj < 4; jj++){
      f[ii] += sign*ginv(y[0],y[1],y[2],y[3],ii,jj)*y[jj+4];
      }}
    ii=4;
    int v;
  for(ii=4; ii < 8; ii++){
    for(jj=0; jj < 4; jj++){
      for(kk=0; kk < 4; kk++){
	v=ii-4;
	f[ii] += sign*(-0.5)*ginvderv(y[0],y[1],y[2],y[3],jj,kk,v)*y[jj+4]*y[kk+4];
	    }}}

  return GSL_SUCCESS;
}



int
main (void)
{
  double mu = -1.0;
  //  int *nothing=NULL;
  gsl_odeiv2_system sys = {geodesics, NULL, 8, &mu};
     
       gsl_odeiv2_driver * d = 
         gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
					1e-6, 1e-6, 0.0);
       int i,ss,qq;
       double tau = 0.0, tauf = 1000, tau0 = 0;

       double p[DIM]={0,0,0,0};
       for(ss = 0; ss < DIM; ss++){
	 for(qq = 0; qq < DIM; qq++){
	    p[ss] += g(x0[0],x0[1],x0[2],x0[3],ss,qq)*v0[qq];
	 }}

       
       double y[8] = { x0[0],x0[1],x0[2],x0[3],p[0],p[1],p[2],p[3]};

       //printf("%g %g %g %g \n", v0[0],v0[1],v0[2],v0[3]);
       int Nsample = 200;


       //prints (tau, t, r, theta,phi)
       // will stop if it gets too close to the horizon.
       for (i = 0; i < Nsample; i++)
         {
           double taui = tau0 + (i * (tauf - tau0) / Nsample);
           int status = gsl_odeiv2_driver_apply (d, &tau, taui, y);
	   
           if (status != GSL_SUCCESS  || (y[1] < (2+0.1) * M ))
	     {
	       printf ("too close to horizon or gsl status %i \n", status);
	       break;
	     }
	   printf ("%g %g %g %g %g\n",tau, y[0], y[1], y[2],y[3]);
         }
     
       gsl_odeiv2_driver_free (d);
       return 0;
}
