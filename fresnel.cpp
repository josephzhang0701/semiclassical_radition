#include <math.h>
#include "complex.h"
#define EPS 6.0e-8//6.0e-8
#define MAXIT 100//100
#define FPMIN 1.0e-30
#define XMIN 1.5
//#define PI 3.1415927
//#define PIBY2 (PI/2.0)
double PI=M_PI;
double PIBY2=(PI/2.0);
//Here EPS is the relative error; MAXIT is the maximum number of iterations allowed; FPMIN
//is a number near the smallest representable doubleing-point number; XMIN is the dividing line
//between using the series and continued fraction.
#define TRUE 1
#define ONE Complex(1.0,0.0)
#include <stdio.h>
#include <stdlib.h>
#include "fresnel.h"
void frenel(double x, double *c, double *s)
//Computes the Fresnel integrals S(x) and C(x) for all real x.
{
void nrerror(char error_text[]);
int k,n,odd;
double a,ax,fact,pix2,sign,sum,sumc,sums,term,test;
fcomplex b,cc,d,h,del,cs;
ax=fabs(x);
if (ax < sqrt(FPMIN)) { //Special case: avoid failure of convergence
*s=0.0; //test because of underflow.
*c=ax;
} else if (ax <= XMIN) { //Evaluate both series simultaneously.
sum=sums=0.0;
sumc=ax;
sign=1.0;
fact=PIBY2*ax*ax;
odd=TRUE;
term=ax;
n=3;
for (k=1;k<=MAXIT;k++) {
term *= fact/k;
sum += sign*term/n;
test=fabs(sum)*EPS;
if (odd) {
sign = -sign;
sums=sum;
sum=sumc;
} else {
sumc=sum;
sum=sums;
}
if (term < test) break;
odd=!odd;
n += 2;
}
if (k > MAXIT) nrerror("series failed in frenel");
*s=sums;
*c=sumc;
} else { //Evaluate continued fraction by modified
pix2=PI*ax*ax; //Lentz's method (§5.2).
b=Complex(1.0,-pix2);
cc=Complex(1.0/FPMIN,0.0);
d=h=Cdiv(ONE,b);
n = -1;
for (k=2;k<=MAXIT;k++) {
n += 2;
a = -n*(n+1);
b=Cadd(b,Complex(4.0,0.0));
d=Cdiv(ONE,Cadd(RCmul(a,d),b));// Denominators cannot be zero.
cc=Cadd(b,Cdiv(Complex(a,0.0),cc));
del=Cmul(cc,d);
h=Cmul(h,del);
if (fabs(del.r-1.0)+fabs(del.i) < EPS) break;
}
if (k > MAXIT) nrerror("cf failed in frenel");
h=Cmul(Complex(ax,-ax),h);
cs=Cmul(Complex(0.5,0.5),
Csub(ONE,Cmul(Complex(cos(0.5*pix2),sin(0.5*pix2)),h)));
*c=cs.r;
*s=cs.i;
}
if (x < 0.0) { //Use antisymmetry.
*c = -(*c);
*s = -(*s);
}
}



void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{fprintf(stderr,"Numerical Recipes run-time error...\n");
fprintf(stderr,"%s\n",error_text);
fprintf(stderr,"...now exiting to system...\n");
exit(1);
}

double sinc(double x)
{
    if (x == 0.0)
        return 1.0;
    return sin(x)/x;
}


double sign (double x)
  {
   double y;
   if (x>0.0) y=1.0;
   else if (x<0.0) y=-1.0;
   else y=0.0;
   return y;
  }

double heaviside (double x)
  { double y= 0.5 * (sign(x) + 1.0);
    return y;
  }

