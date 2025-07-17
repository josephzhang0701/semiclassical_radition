#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "const.h"
#include "rndinit.h"
//#include "fields.h"
#include <cmath>

void rndinit (double *y0array, double *z0array, double *x0array, double* Ekinarray, double *anglearray, double *anglearray2, int particle_nr2) {
  
  const gsl_rng_type * T1;
  gsl_rng * r1;
  gsl_rng_env_setup();
  
  T1=gsl_rng_default;
  r1=gsl_rng_alloc(T1);
  
  gsl_rng_set(r1, 314159265*(particle_nr2+1));
  
  double r_ro,r_ph;
  for (int i=0; i<NP; ++i) {
    do{
      do{
        r_ro = gsl_ran_gaussian (r1,5.0*lambda0);
        } while (r_ro==0);   
        r_ph = gsl_ran_flat (r1, 0.0,2.0*pi);
        y0array[i]=r_ro*cos(r_ph);
        z0array[i]=r_ro*sin(r_ph);
        x0array[i]=0.0;
        //x0array[i]=gsl_ran_flat (r1, cylinder_start,cylinder_end);     ///区间[cylinder_start, cylinder_end]中，均匀分布的一个随机数
    }   while (r_ro*r_ro>=cylinderradius*cylinderradius);    ///确保生成的 (x0array[i], y0array[i]) 点在圆内

    Ekinarray[i]=Ekinmean + gsl_ran_gaussian (r1,Ekinsigma);	
    anglearray[i]=anglemean + gsl_ran_gaussian (r1,anglesigma);	
   // anglearray[i]=anglemean + gsl_ran_flat (r1,-anglesigma,anglesigma);
    anglearray2[i]=gsl_ran_flat (r1, 0.0,2.0*pi);		
    
  }
  
  gsl_rng_free (r1);
}
