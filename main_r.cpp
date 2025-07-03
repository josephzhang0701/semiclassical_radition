#include <iostream>
#include <iomanip>
#include <fstream> 
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <string>
#include <sstream>

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_statistics.h>

#include "const.h"
#include "chi_e.h"
#include "fresnel.h"
#include "vector3d.h"

using namespace std;
using namespace SCATMECH;

int
func (double t, const double g[], double dg[],                             //Definition ODEs
      void *params)
{

 double Ex,Ey,Ez,Bx,By,Bz;
 double theta=ww0*t+k*g[4];

        Ex=E0*cos(theta);
        Ey=0.0;
        Ez=0.0;
        Bx=0.0;
        By=-Ex;
        Bz=0.0;


 ret_em_prob parameter_em_per_step;
 parameter_em_per_step=Wr(g,t);

double vb=sqrt(g[1]*g[1]+g[3]*g[3]+g[5]*g[5]);
double grg=sqrt(vb*vb+1.0);       // gamma

  dg[0] = g[1]*c/grg;                                              
  dg[2] = g[3]*c/grg;                               				
  dg[4] = g[5]*c/grg;   
  dg[1] = eq/me/c*(Ex+g[3]/grg*Bz-g[5]/grg*By);
  dg[3] = eq/me/c*(Ey+g[5]/grg*Bx-g[1]/grg*Bz);
  dg[5] = eq/me/c*(Ez+g[1]/grg*By-g[3]/grg*Bx);

  return GSL_SUCCESS;
}


vector<double> linespace(double start, double ed, int num) {
    // catch rarely, throw often
    if (num < 2) {
        throw new exception();
    }
    int partitions = num - 1;
    vector<double> pts;
    // length of each segment    
    double length = (ed - start) / partitions; 
    // first, not to change
    pts.push_back(start);
    for (int i = 1; i < num - 1; i ++) {
        pts.push_back(start + i * length);
    }
    // last, not to change
    pts.push_back(ed);
    return pts;
}

int  main(int argc, char *argv[])
{

  int particle_nr(std::atoi(argv[1]));

  std::stringstream filename1;
  filename1 << "data/pbw_" <<particle_nr<< ".dat";
  ofstream fout1;
  fout1.open (filename1.str().c_str());

  std::stringstream filename2;
  filename2 << "data/pbx_" <<particle_nr<< ".dat";
  ofstream fout2;
  fout2.open (filename2.str().c_str());

  std::stringstream filename3;
  filename3 << "data/trace_" <<particle_nr<< ".dat";
  ofstream fout3;
  fout3.open (filename3.str().c_str());

 double rr0=(Ekinmean/(me * c * c)) + 1.0;
 int i_m1=100;
 int j_m1=100;
 int k_m1=100;
vector<double> x_list=linespace(-(7.0)/rr0,(7.0)/rr0,i_m1);
vector<double> y_list=linespace(-(7.0)/rr0,(7.0)/rr0,j_m1);
vector<double> w_list=linespace(0.0000001*rr0*me*c*c/hb,0.9999999*rr0*me*c*c/hb,k_m1);

double dtx=x_list[1]-x_list[0];
double dty=y_list[1]-y_list[0];
double dw=w_list[1]-w_list[0];

  double fx,fy,fz,kkx,g1,n_a,nz;
  double TIrx[i_m1][j_m1][k_m1];
  double TIry[i_m1][j_m1][k_m1];
  double TIrz[i_m1][j_m1][k_m1];
  double TIix[i_m1][j_m1][k_m1];
  double TIiy[i_m1][j_m1][k_m1];
  double TIiz[i_m1][j_m1][k_m1];
  double TIr_J[i_m1][j_m1][k_m1];
  double TIi_J[i_m1][j_m1][k_m1];

for (int ki=0;ki<k_m1;ki++){
  for (int i=0;i<i_m1;i++){
       for (int j=0;j<j_m1;j++){
           TIrx[i][j][ki]=0.0;
           TIry[i][j][ki]=0.0;
           TIrz[i][j][ki]=0.0;
           TIix[i][j][ki]=0.0;
           TIiy[i][j][ki]=0.0;
           TIiz[i][j][ki]=0.0;
           TIr_J[i][j][ki]=0.0;
           TIi_J[i][j][ki]=0.0;
        }
  }
}


 int nnt=2;
 for (int tth=1; tth<nnt; tth++) 
 { 
   
    if (argc!=2) {
    std::cerr << "Wrong numer of arguments.\n";
    exit(EXIT_FAILURE);
    }
   
   int particle_nr2(std::atoi(argv[1]));
   particle_nr2=tth+(particle_nr-1)*(nnt-1);

 /////////////////////////////////////////////////////

  const gsl_odeiv_step_type * T 
   = gsl_odeiv_step_rkf45;
    
  gsl_odeiv_step * sdgl 
    = gsl_odeiv_step_alloc (T, 6);
  gsl_odeiv_control * cdgl 
    = gsl_odeiv_control_y_new (1e-06, 1e-8); ///////// Please check the precision carefully 
  gsl_odeiv_evolve * edgl 
    = gsl_odeiv_evolve_alloc (6);
    
      
  gsl_odeiv_system sys = {func, 0, 6, 0};

  double ww;
  double ti = 0.0;
  double tf = 7.5*tp0;//100.0*tp0;//tp0*120.0;
  double dt = 0.001*tp0;
  double t;
  double g[6];
  double x0array[NP];
  double y0array[NP];
  double z0array[NP];
  double Ekinarray[NP];
  double anglearray[NP];
  double anglearray2[NP];
  double aa=dt;   ////////////////// Be careful and please check the value can make some differences???
  double h=aa;
  double Pe;                ////// electron momentum
  double grg;
  double TT=1.0e-3;

  const gsl_rng_type * T1;
  gsl_rng * r1;
  gsl_rng_env_setup();
  
  T1=gsl_rng_default;
  r1=gsl_rng_alloc(T1);
  
  gsl_rng_set(r1, 314159265*(particle_nr2+1));


 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      int i=0;
      rr0 = (Ekinmean/(me * c * c)) + 1.0;
      grg=rr0;
      double p0 = sqrt(rr0*rr0-1.0);
      g[0] = 0.0;
      g[2] = 0.0;
      g[4] = 0.0;
      g[1] = 0.0;
      g[3] = 0.0;
      g[5] = p0;
    
      double Beta0,vx,vy,vz,de_z,de_vz;
      Beta0=1.0-1.0/2.0/grg/grg;
     
      double ax,ay,az,nx,ny;  
      double www,tx,ty;
	  
      double Ex,theta;
        
              t=ti;
          while (t <= tf)
	     {

                                      
        //******************electron trajectory***********************************
 
                 double tff = t+dt;
	         while (t < tff)
	             {
	               int status = gsl_odeiv_evolve_apply (edgl, cdgl, sdgl,&sys, &t, tff,&h, g);
	               if (status != GSL_SUCCESS)
		             break;
	             } 

	              gsl_odeiv_evolve_reset (edgl);
	              gsl_odeiv_step_reset (sdgl);
                h=dt;
                Pe = sqrt(g[1]*g[1]+g[3]*g[3]+g[5]*g[5]);
                grg=sqrt(Pe*Pe+1.0);       // gamma
//acceleration
                ret_em_prob parameter_em_per_step=Wr(g,t);
                ax=parameter_em_per_step.ax;
                ay=parameter_em_per_step.ay;
                az=parameter_em_per_step.az;
                vx=g[1]/grg;
                vy=g[3]/grg;
                vz=g[5]/grg;

                fout3<<t/tp0<<" "<<g[0]<<" "<<g[2]<<" "<<g[4]<<" "<<vx<<" "<<vy<<" "<<vz<<endl;

		de_z=g[4]/c-Beta0*t;
		de_vz=vz-Beta0;

      for (int i=0;i<i_m1;i++)
          {
            tx=x_list[i];//[i];//tx=nx; 

      for (int j=0;j<j_m1;j++)
          {
            ty=y_list[j];//phi=y00+dy*j;//ty=ny;

	    for (int ki=0;ki<k_m1;ki++)
          {     
		
          nz=1.0-(tx*tx+ty*ty)/2.0;
	        www=hb*w_list[ki]/me/c/c;
          ww=rr0/(rr0-www)*w_list[ki];//w'=ele/pos*w;
             
	        g1=1.0/pow(1.0/2.0/rr0/rr0+(tx*tx+ty*ty)/2.0-de_vz-tx*vx-ty*vy,2.0);
	        fx=g1*(ty*((tx-vx)*ay-(ty-vy)*ax)-(1.0/2.0/rr0/rr0-(tx*tx+ty*ty)/2.0-de_vz)*ax+(tx-vx)*az);
          fy=g1*((ty-vy)*az-(1.0/2.0/rr0/rr0-(tx*tx+ty*ty)/2.0-de_vz)*ay-tx*((tx-vx)*ay-(ty-vy)*ax));
          fz=g1*(tx*((1.0/2.0/rr0/rr0-(tx*tx+ty*ty)/2.0-de_vz)*ax-(tx-vx)*az)-ty*((ty-vy)*az-1.0/2.0/rr0/rr0-(tx*tx+ty*ty)/2.0-de_vz)*ay);
	        n_a=g1*(tx*ax+ty*ay+az);//dot(n,a)*g;
	        kkx=ww*((1.0/2.0/rr0/rr0+(tx*tx+ty*ty)/2.0)*t-de_z-tx*g[0]/c-ty*g[2]/c);
                 
           TIrx[i][j][ki]=TIrx[i][j][ki]+fx*cos(kkx)*dt;
           TIry[i][j][ki]=TIry[i][j][ki]+fy*cos(kkx)*dt;
           TIrz[i][j][ki]=TIrz[i][j][ki]+fz*cos(kkx)*dt;
           TIix[i][j][ki]=TIix[i][j][ki]+fx*sin(kkx)*dt;
           TIiy[i][j][ki]=TIiy[i][j][ki]+fy*sin(kkx)*dt;
           TIiz[i][j][ki]=TIiz[i][j][ki]+fz*sin(kkx)*dt;
           TIr_J[i][j][ki]=TIr_J[i][j][ki]+n_a*cos(kkx)*dt;
           TIi_J[i][j][ki]=TIi_J[i][j][ki]+n_a*sin(kkx)*dt;
		   
             }//ki
			   }//j
		   }//i
	  }//t

  gsl_odeiv_evolve_free (edgl);
  gsl_odeiv_control_free (cdgl);
  gsl_odeiv_step_free (sdgl);

 }

 Vector3D<double> Ireal,Iimag;
 double Breal,Bimag,Jreal,Jimag,PP,Bsqua,Areal,Aimag,Asqua;
 double rr02,C,D1,D2,www,tx,ty,domega;


 for (int k=0;k<k_m1;k++){
       PP=0;

        for (int i=0;i<i_m1;i++){
            tx=x_list[i];

              for (int j=0;j<j_m1;j++){
               ty=y_list[j];

               domega=dtx*dty;
               www=hb*w_list[k]/me/c/c;
               rr02=rr0-www;

               C=eq*eq/4.0/pi/pi/hb;
               D1=(pow(rr02,2)+pow(rr0,2))/2.0/pow(rr0,2);
               D2=pow(www,2)/2.0/pow(rr0,4);

               Ireal={TIrx[i][j][k],TIry[i][j][k],TIrz[i][j][k]};
               Iimag={TIix[i][j][k],TIiy[i][j][k],TIiz[i][j][k]};
               Jreal=TIr_J[i][j][k];
               Jimag=TIi_J[i][j][k];

               Areal=dot(Ireal,Ireal);
               Aimag=dot(Iimag,Iimag);
               Breal=Jreal*Jreal;
               Bimag=Jimag*Jimag;
               Asqua=Areal+Aimag;
               Bsqua=Breal+Bimag;

               PP=PP+C*(D1*Asqua+D2*Bsqua)*domega;
                }
          }
       fout1 << std::fixed << std::scientific << std::setprecision(15) <<hb*w_list[k]/me/c/c/rr0<<" "<<PP<<endl;
 }


for (int i=0;i<i_m1;i++){
    tx=x_list[i];

    for (int j=0;j<j_m1;j++){
        ty=y_list[j];
        PP=0.0;

		    for (int k=0;k<k_m1-1;k++){

               www=hb*w_list[k]/me/c/c/hb;
               rr02=rr0-www;

               C=eq*eq/4.0/pi/pi;
               D1=(pow(rr02,2)+pow(rr0,2))/2.0/pow(rr0,2);
               D2=pow(www,2)/2.0/pow(rr0,4);

               Ireal={TIrx[i][j][k],TIry[i][j][k],TIrz[i][j][k]};
               Iimag={TIix[i][j][k],TIiy[i][j][k],TIiz[i][j][k]};
               Jreal=TIr_J[i][j][k];
               Jimag=TIi_J[i][j][k];

               Areal=dot(Ireal,Ireal);
               Aimag=dot(Iimag,Iimag);
               Breal=Jreal*Jreal;
               Bimag=Jimag*Jimag;
               Asqua=Areal+Aimag;
               Bsqua=Breal+Bimag;

               PP=PP+C*(D1*Asqua+D2*Bsqua)*dw;

                }
    fout2<<x_list[i]<<" "<<y_list[j]<<" "<<PP<<endl;
          }
}

fout1.close();
fout2.close();
fout3.close();
return 0; 
}
