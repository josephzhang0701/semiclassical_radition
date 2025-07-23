#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <string>

#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_statistics.h>

#include "const.h"
#include "vector3d.h" //点乘叉乘并矢
#include "linespace.h"
//下面是祖传代码, 不一定用到
#include "chi_e.h"
#include "fresnel.h"

using namespace std;
using namespace SCATMECH;

inline void get_fields(double t, const double g[], double& Ex, double& Ey, double& Ez, double& Bx, double& By, double& Bz) {
    double theta = ww0*t + k*g[4];
    Ex = E0 * cos(theta);
    Ey = 0.0;
    Ez = 0.0;
    Bx = 0.0;
    By = -Ex;
    Bz = 0.0;
}

// ODE definition
int func(double t, const double g[], double dg[], void* params) {
    double Ex, Ey, Ez, Bx, By, Bz;
    get_fields(t, g, Ex, Ey, Ez, Bx, By, Bz);

    ret_em_prob parameter_em_per_step;
    parameter_em_per_step = Wr(g, t, Ex, Ey, Ez, Bx, By, Bz);

    double vb = sqrt(g[1]*g[1] + g[3]*g[3] + g[5]*g[5]);
    double grg = sqrt(vb*vb + 1.0);

    dg[0] = g[1]*c/grg;
    dg[2] = g[3]*c/grg;
    dg[4] = g[5]*c/grg;
    dg[1] = eq/me/c*(Ex + g[3]/grg*Bz - g[5]/grg*By);
    dg[3] = eq/me/c*(Ey + g[5]/grg*Bx - g[1]/grg*Bz);
    dg[5] = eq/me/c*(Ez + g[1]/grg*By - g[3]/grg*Bx);

    return GSL_SUCCESS;
}

//vector<double> linespace(double start, double ed, int num) {
//    // catch rarely, throw often
//    if (num < 2) {
//        throw new exception();
//    }
//    int partitions = num - 1;
//    vector<double> pts;
//    // length of each segment
//    double length = (ed - start) / partitions;
//    // first, not to change
//    pts.push_back(start);
//    for (int i = 1; i < num - 1; i ++) {
//        pts.push_back(start + i * length);
//    }
//    // last, not to change
//    pts.push_back(ed);
//    return pts;
//}

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Wrong number of arguments. Usage: ./irun.x <particle_nr>\n";
        return EXIT_FAILURE;
    }

    int particle_nr = atoi(argv[1]);

    // prepare output files
    ostringstream fn1, fn2, fn3;
    fn1 << "data/pbw_" << particle_nr << ".dat";
    fn2 << "data/pbx_" << particle_nr << ".dat";
    fn3 << "data/trace_" << particle_nr << ".dat";
    ofstream fout1(fn1.str()), fout2(fn2.str()), fout3(fn3.str());

    // grid sizes
    int i_m1 = 200;
    int j_m1 = 200;
    int k_m1 = 200;

    // linspace grids
    double rr0 = (Ekinmean/(me*c*c)) + 1.0;     //TODO: rr0好像是gamma_initial
    //TODO: 网格选择的讲究, Compton Scattering光斑大小
//    vector<double> x_list = linespace(-4.5/rr0, 4.5/rr0, i_m1);
//    vector<double> y_list = linespace(-4.5/rr0, 4.5/rr0, j_m1);
//    vector<double> w_list = linespace(1e-7*rr0*me*c*c/hb, 0.9999999*rr0*me*c*c/hb, k_m1);   //避免除以0


    auto x_list = linespace(-4.5/rr0, 4.5/rr0, i_m1).generate();
    auto y_list = linespace(-4.5/rr0, 4.5/rr0, j_m1).generate();
    auto w_list = linespace(
            1e-7*rr0*me*c*c/hb,
            0.9999999*rr0*me*c*c/hb,
            k_m1
    ).generate();


    double dtx = x_list[1] - x_list[0];
    double dty = y_list[1] - y_list[0];
    double dw  = w_list[1] - w_list[0];

    // allocate on heap
    size_t total = size_t(i_m1) * j_m1 * k_m1;
    vector<double> TIrx(total), TIry(total), TIrz(total);
    vector<double> TIix(total), TIiy(total), TIiz(total);
    vector<double> TIr_J(total), TIi_J(total);

    auto IDX = [j_m1, k_m1](int i, int j, int k) -> size_t {
        return (size_t(i) * j_m1 + j) * k_m1 + k;
    };

    // zero arrays
    for (int ki = 0; ki < k_m1; ++ki)
        for (int i = 0; i < i_m1; ++i)
            for (int j = 0; j < j_m1; ++j) {
                size_t idx = IDX(i,j,ki);
                TIrx[idx] = TIry[idx] = TIrz[idx] = 0.0;
                TIix[idx] = TIiy[idx] = TIiz[idx] = 0.0;
                TIr_J[idx] = TIi_J[idx] = 0.0;
            }

    int nnt = 1;    //hard code 粒子数, 电子片位置可以设计
    // TODO: 大量粒子模拟的时候需要改成外面传入
    for (int tth = 1; tth <= nnt; ++tth) {
        // tth 保证所有核心一起启动
        //particle_nr2 是核心数, 提交任务传来的
        int particle_nr2 = tth + (particle_nr - 1)*(nnt - 1);

        // GSL ODE setup
        const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
        gsl_odeiv_step    *sdgl = gsl_odeiv_step_alloc(T, 6);
        gsl_odeiv_control *cdgl = gsl_odeiv_control_y_new(1e-6, 1e-8);
        gsl_odeiv_evolve  *edgl = gsl_odeiv_evolve_alloc(6);
        gsl_odeiv_system    sys = {func, nullptr, 6, nullptr};

        double ww;
        double ti = 0.0;
        double tf = 7.5*tp0;        //tp0==lambda0/c
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
        double fx,fy,fz,kkx,g1,n_a;

        // 随机数撒种, 和电子片有关
        // TODO: 注意电子片
        const gsl_rng_type * T1;        //声明一个指向GSL随机数类型的指针。
        gsl_rng * r1;
        gsl_rng_env_setup();

        T1=gsl_rng_default;
        r1=gsl_rng_alloc(T1);

        gsl_rng_set(r1, 314159265*(particle_nr2+1));        //设定随机数生成器的“种子”，使得每个particle_nr2对应不同的随机数序列。


        // initial conditions
        int i=0;
        rr0 = (Ekinmean/(me * c * c)) + 1.0;
        grg=rr0;
        double p0 = sqrt(rr0*rr0-1.0);
        // 电子ini x
        g[0] = 0.0;
        g[2] = 0.0;
        g[4] = 0.0;
        // 电子ini p
        g[1] = 0.0;
        g[3] = 0.0;
        g[5] = p0;

        double Beta0,vx,vy,vz,de_z,de_vz;
        Beta0=1.0-1.0/2.0/grg/grg;

        double ax, ay, az, nx, ny, nz;
        double www,tx,ty;

        double Ex, Ey, Ez, Bx, By, Bz, theta;

        t=ti;

        while (t <= tf) {
            //******************electron trajectory***********************************
            double tff = t + dt;
            while (t < tff) {
                int status = gsl_odeiv_evolve_apply(edgl, cdgl, sdgl, &sys, &t, tff, &h, g);
                if (status != GSL_SUCCESS) break;
            }
            gsl_odeiv_evolve_reset(edgl);
            gsl_odeiv_step_reset(sdgl);
            h = dt;

            // trajectory output
            Pe = sqrt(g[1]*g[1] + g[3]*g[3] + g[5]*g[5]);
            grg = sqrt(g[1]*g[1] + g[3]*g[3] + g[5]*g[5] + 1.0); //gamma

            get_fields(t, g, Ex, Ey, Ez, Bx, By, Bz);
            ret_em_prob parameter_em_per_step = Wr(g, t, Ex, Ey, Ez, Bx, By, Bz);
            ax=parameter_em_per_step.ax;
            ay=parameter_em_per_step.ay;
            az=parameter_em_per_step.az;
            vx = g[1]/grg; vy = g[3]/grg; vz = g[5]/grg;
            fout3 << t/tp0 << " " << g[0] << " " << g[2] << " " << g[4]
                  << " " << vx << " " << vy << " " << vz << "\n";

            de_z=g[4]/c-Beta0*t;
            de_vz=vz-Beta0;

            // accumulate integrals 计算辐射的角分布呢
            for (int i = 0; i < i_m1; ++i) {
                tx = x_list[i];
                for (int j = 0; j < j_m1; ++j) {
                    ty = y_list[j];
                    for (int ki2 = 0; ki2 < k_m1; ++ki2) {
                        size_t idx = IDX(i,j,ki2);
                        nz = 1.0 - (tx*tx+ty*ty)/2.0;
                        www = hb*w_list[ki2]/me/c/c;
                        ww = rr0/(rr0-www)*w_list[ki2];
                        g1 = 1.0/pow(1.0/2.0/rr0/rr0 + (tx*tx+ty*ty)/2.0 - de_vz - tx*vx - ty*vy, 2.0);
                        fx = g1*(ty*((tx-vx)*ay - (ty-vy)*ax)
                                        - (1.0/2.0/rr0/rr0 - (tx*tx+ty*ty)/2.0 - de_vz)*ax
                                        + (tx-vx)*az);
                        fy = g1*((ty-vy)*az
                                        - (1.0/2.0/rr0/rr0 - (tx*tx+ty*ty)/2.0 - de_vz)*ay
                                        - tx*((tx-vx)*ay - (ty-vy)*ax));
                        fz = g1*(tx*((1.0/2.0/rr0/rr0 - (tx*tx+ty*ty)/2.0 - de_vz)*ax - (tx-vx)*az)
                                        - ty*((ty-vy)*az
                                              - (1.0/2.0/rr0/rr0 - (tx*tx+ty*ty)/2.0 - de_vz)*ay));
                        n_a = g1*(tx*ax + ty*ay + az);
                        kkx = ww*((1.0/2.0/rr0/rr0 + (tx*tx+ty*ty)/2.0)*t - de_z
                                         - tx*g[0]/c - ty*g[2]/c);
                        TIrx[idx] += fx * cos(kkx) * dt;
                        TIry[idx] += fy * cos(kkx) * dt;
                        TIrz[idx] += fz * cos(kkx) * dt;
                        TIix[idx] += fx * sin(kkx) * dt;
                        TIiy[idx] += fy * sin(kkx) * dt;
                        TIiz[idx] += fz * sin(kkx) * dt;
                        TIr_J[idx] += n_a * cos(kkx) * dt;
                        TIi_J[idx] += n_a * sin(kkx) * dt;
                    }  //ki2
                } //j
            } //i
        } //t
        gsl_odeiv_evolve_free(edgl);
        gsl_odeiv_control_free(cdgl);
        gsl_odeiv_step_free(sdgl);
    }

    Vector3D<double> Ireal,Iimag;
    double Breal,Bimag,Jreal,Jimag,PP,Areal,Aimag;
    double rr02,C,D1,D2,www,tx,ty,domega;

    // compute and write spectra 能谱
    for (int ki = 0; ki < k_m1; ++ki) {
        PP = 0.0;
        for (int i = 0; i < i_m1; ++i) {
            tx = x_list[i]; //空间的网格
            for (int j = 0; j < j_m1; ++j) {
                domega = dtx * dty;     //小角度近似
                size_t idx = IDX(i,j,ki);
                rr02 = rr0 - hb*w_list[ki]/me/c/c;
                C = eq*eq/4.0/pi/pi/hb;
                D1 = (pow(rr02,2) + pow(rr0,2))/(2.0*pow(rr0,2));
                D2 = pow(hb*w_list[ki]/me/c/c,2)/(2.0*pow(rr0,4));
                Ireal = Vector3D<double>{TIrx[idx], TIry[idx], TIrz[idx]};
                Iimag = Vector3D<double>{TIix[idx], TIiy[idx], TIiz[idx]};
                Jreal = TIr_J[idx], Jimag = TIi_J[idx];
                Areal = dot(Ireal,Ireal);
                Aimag = dot(Iimag,Iimag);
                Breal = Jreal*Jreal;
                Bimag = Jimag*Jimag;
                PP += C*(D1*(Areal+Aimag) + D2*(Breal+Bimag)) * domega;
            }
        }
        fout1 << fixed << scientific << setprecision(15)
              << hb*w_list[ki]/me/c/c/rr0 << " " << PP << endl;
    }

    // angular distribution
    for (int i = 0; i < i_m1; ++i) {
        tx = x_list[i];
        for (int j = 0; j < j_m1; ++j) {
            PP = 0.0;
            for (int ki = 0; ki < k_m1 - 1; ++ki) {
                size_t idx = IDX(i,j,ki);
                rr02 = rr0 - hb*w_list[ki]/me/c/c;
                C = eq*eq/4.0/pi/pi;
                D1 = (pow(rr02,2) + pow(rr0,2))/(2.0*pow(rr0,2));
                D2 = pow(hb*w_list[ki]/me/c/c,2)/(2.0*pow(rr0,4));
                Ireal = Vector3D<double>{TIrx[idx], TIry[idx], TIrz[idx]};
                Iimag = Vector3D<double>{TIix[idx], TIiy[idx], TIiz[idx]};
                Jreal = TIr_J[idx], Jimag = TIi_J[idx];
                Areal = dot(Ireal,Ireal);
                Aimag = dot(Iimag,Iimag);
                Breal = Jreal*Jreal;
                Bimag = Jimag*Jimag;
                PP += C*(D1*(Areal+Aimag) + D2*(Breal+Bimag)) * dw;
            }
            fout2 << tx << " " << y_list[j] << " " << PP << endl;
        }
    }

    fout1.close();
    fout2.close();
    fout3.close();
    return 0;
}