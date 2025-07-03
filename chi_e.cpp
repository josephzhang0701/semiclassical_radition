#include "const.h"
#include "chi_e.h"
#include <math.h>
#include <cmath>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_integration.h>
#include <iostream>
#include "fresnel.h"


ret_em_prob Wr(const double* g,double t){

//loading the fields

double Ex,Ey,Ez,Bx,By,Bz;
ret_em_prob internal_result;
double theta=ww0*t+k*g[4];

        Ex=E0*cos(theta);
        Ey=0.0;
        Ez=0.0;
        Bx=0.0;
        By=-Ex;
        Bz=0.0;

double vb=sqrt(g[1]*g[1]+g[3]*g[3]+g[5]*g[5]);
double grg=sqrt(vb*vb+1.0);       // gamma

double j1=g[1]*Ex+g[3]*Ey+g[5]*Ez;//p*E
double j2=grg*Ex+g[3]*Bz-g[5]*By;//grg*E+p^B
double j3=grg*Ey+g[5]*Bx-g[1]*Bz;
double j4=grg*Ez+g[1]*By-g[3]*Bx;


  double chi_r=hb*sqrt(eq*eq)/me/me/c/c/c*sqrt(j2*j2+j3*j3+j4*j4-j1*j1);

double a1= eq/me/c/grg*(j2/grg-j1*g[1]/grg/grg);

double a2= eq/me/c/grg*(j3/grg-j1*g[3]/grg/grg);

double a3= eq/me/c/grg*(j4/grg-j1*g[5]/grg/grg);



internal_result.chi_r =  chi_r;
internal_result.ax =  a1;
internal_result.ay =  a2;
internal_result.az =  a3;
return internal_result;

}

