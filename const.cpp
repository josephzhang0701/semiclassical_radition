#include "const.h"
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>


double pi=M_PI;                        //Defininition of Pi  
double lambda0=5.2759e-05;//0.4e-4;//5.2759e-05;//0.8e-4;     //wavelength
double k=2.0*pi/lambda0;   // wave vector
double c=2.99792458e10;  // light velocity in vacuum
double ww0=k*c;          // laser frequency
double tp0=lambda0/c;  //2.6685fs  
double tp=8.0*tp0;//3.0*tp0;//10.0*tp0;// // laser pulse duration---1ps
double phase0=0.0/180.0*pi;  // initial phase of laser pulse
//global parameters

double electronm=9.1093826e-28;     //electron mass
double eq0=4.80320441e-10;          //electron electric quantity
double mev=1.60217653e-6;           // MeV in the unit of Joule

double me=1.0*electronm;    // mass of particle
double eq=-1.0*eq0;        // electric quantity of particle
double hb=1.054571726*1e-34*1e7;//Reduced planck constant
double af=pow(eq,2.0)/hb/c;   //fine structure constant
double lambdac=hb/me/c; //compton wavelength
double tauc=lambdac/c; //compton time
double w0=5.0*lambda0;//1.0*lambda0;   // laser radius
double zr=k*pow(w0,2.0)/2.0;
double eps0=w0/zr;   //lax fields variable
double s=ww0*tp/(2.0*sqrt(2.0*log(2.0)));//s=ww0*tp/sqrt(2.0*log(2.0));    //s= 8.004669384955495 if lambda0=1e-4;
double shift=3.0;//5.0;
double theta_m=shift*s;
double q=0.707;
double q1=1.0/(2.0*sqrt(2.0));//1.0/2.0*(sqrt(2.0));//1.0;//2.0;
double E0=q*me*c*ww0/sqrt(eq*eq);
//double E02=q1*me*c*ww0/sqrt(eq*eq);
double ey_ratio =1.0;
//monte carlo distribution of electrons

double k2=2.0*k;
double ww02=k2*c;
double phase02=90.0/180.0*pi;
double w02=w0*1.0; //laser radius
double zr2=k2*pow(w02,2.0)/2.0;
double eps02=w02/zr2;
double E02=(0.5)*E0;
   			    
double cylinderlength=10.0*lambda0; //length of the cylinder 1mm
double cylinder_start=0.0;   // the start position of cylinder in longitudinal direction
double cylinder_end=5.0*lambda0;   // the end position of cylinder in longitudinal direction 
double cylinderradius=1.0*lambda0;   // radius of the cylinder 
double radius_y=1.0*lambda0;                
double Ekinmean=20000.0*mev;//8000.0*mev;           // mean initial kinetic energy of the electrons 
double Wk=1022.0*mev/me/c/c;//8000.0*mev/me/c/c;            
double Ekinsigma=0.0;//Ekinmean*0.02;  // standard deviation of initial kin. energies 
double anglemean=180.0/180.0*pi;           // mean initial injection angle of electron bunch 
double anglesigma=0.0;//0.02;  // standard deviation of initial injection angle of electron bunch 


int NP=1;
int NS=15e3;		


















