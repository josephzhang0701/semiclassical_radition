#ifndef CONST_H
#define CONST_H

#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
//全高斯单位制 cm g s

constexpr double pi=M_PI;                        //Defininition of Pi
constexpr double lambda0=5.2759e-05;//0.4e-4;//5.2759e-05;//0.8e-4;     //wavelength
constexpr double hb=1.054571726*1e-34*1e7;//Reduced planck constant
constexpr double c=2.99792458e10;  // light velocity in vacuum
constexpr double electronm=9.1093826e-28;     //electron mass
constexpr double eq0=4.80320441e-10;          //electron electric quantity
constexpr double mev=1.60217653e-6;           // MeV in the unit of ergs, 1 Joule == 1e7 ergs

constexpr double me_c2 = 8.1871047868e-07;

const double k=2.0*pi/lambda0;   // wave vector
const double ww0=k*c;          // laser frequency
const double tp0=lambda0/c;  //2.6685fs
const double tp=8.0*tp0;//3.0*tp0;//10.0*tp0;// // laser pulse duration---1ps
const double phase0=0.0/180.0*pi;  // initial phase of laser pulse
//global parameters


const double me=1.0*electronm;    // mass of particle
const double eq=-1.0*eq0;        // electric quantity of particle

const double af=pow(eq,2.0)/hb/c;   //fine structure constant
const double lambdac=hb/me/c; //compton wavelength
const double tauc=lambdac/c; //compton time
const double w0=5.0*lambda0;//1.0*lambda0;   // laser radius
const double zr=k*pow(w0,2.0)/2.0;
const double eps0=w0/zr;   //lax fields variable
const double s=ww0*tp/(2.0*sqrt(2.0*log(2.0)));//s=ww0*tp/sqrt(2.0*log(2.0));    //s= 8.004669384955495 if lambda0=1e-4;
const double shift=3.0;//5.0;
const double theta_m=shift*s;
const double q=0.707;     //峰值强度, 可以改, 无量纲, 文章里的eta用的是平均强度
// rms == 0.707 * peak
const double q1=1.0/(2.0*sqrt(2.0));//1.0/2.0*(sqrt(2.0));//1.0;//2.0;
const double E0=q*me*c*ww0/sqrt(eq*eq);
//const double E02=q1*me*c*ww0/sqrt(eq*eq);
const double ey_ratio =1.0;
//monte carlo distribution of electrons

const double k2=2.0*k;
const double ww02=k2*c;
const double phase02=90.0/180.0*pi;
const double w02=w0*1.0; //laser radius
const double zr2=k2*pow(w02,2.0)/2.0;
const double eps02=w02/zr2;
const double E02=(0.5)*E0;

const double cylinderlength=10.0*lambda0; //length of the cylinder 1mm
const double cylinder_start=0.0;   // the start position of cylinder in longitudinal direction
const double cylinder_end=5.0*lambda0;   // the end position of cylinder in longitudinal direction
const double cylinderradius=1.0*lambda0;   // radius of the cylinder
const double radius_y=1.0*lambda0;
const double Ekinmean=20000.0*mev;        //8000.0*mev; // mean initial kinetic energy of the electrons
const double Wk=1022.0*mev/me/c/c;//8000.0*mev/me/c/c;
const double Ekinsigma=0.0;//Ekinmean*0.02; // standard deviation of initial kin. energies
const double anglemean=180.0/180.0*pi;           // mean initial injection angle of an electron bunch
const double anglesigma=0.0;//0.02;  // standard deviation of an initial injection angle of electron bunch


constexpr int NP=1;
constexpr int NS=15e3;


#endif

