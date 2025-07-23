//
// Created by Joseph Zhang on 7/23/25.
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "linespace.h"
#include "const.h"
#include "vector3d.h"  // 点乘叉乘并矢
#include "chi_e.h"

using namespace std;
using namespace SCATMECH;

// Simulation time parameters (editable)
inline double get_t0() { return 0.0; }
inline double get_dt() { return 0.001 * tp0; }
inline double get_tf() { return 7.5   * tp0; }

// Grid parameters
static const int i_m1 = 300;
static const int j_m1 = 300;
static const int freq_m1 = 200;

// Initial gamma for grid (equals el_gamma_ini)
double el_gamma_ini = (Ekinmean / (me * c * c)) + 1.0;

// Linspace grids
static const auto x_list = linespace(-7.5/el_gamma_ini, 7.5/el_gamma_ini, i_m1).generate();
static const auto y_list = linespace(-7.5/el_gamma_ini, 7.5/el_gamma_ini, j_m1).generate();
static const auto freq_list = linespace(
        1e-7*el_gamma_ini*me*c*c/hb,
        (1-1e-7)*el_gamma_ini*me*c*c/hb,
        freq_m1
).generate();

// Grid spacings
static const double dtx = x_list[1] - x_list[0];
static const double dty = y_list[1] - y_list[0];
static const double dw  = freq_list[1] - freq_list[0];

// 3D index helper
inline size_t IDX(int i, int j, int k) {
    return (size_t(i) * j_m1 + j) * freq_m1 + k;
}

// -----------------------------------------------------------------------------
// 1) ParticleInit holds one particle's editable initial state.
struct ParticleInit {
    double x, y, z;    // initial position
    double px, py, pz; // initial momentum
};

// 1) Define all your particles here.  Directly edit these entries.
vector<ParticleInit> initializeParticles() {
    double el_gamma_ini = (Ekinmean/(me * c * c)) + 1.0;
    double p0 = sqrt(el_gamma_ini*el_gamma_ini-1.0);
    return {
            //   x     y     z      px    py    pz
            { 0.0,  0.0,  0.0,   0.0,  0.0,  p0 },   // Particle 0
            { 0.0,  0.0, -0.99*lambda0,   0.0,  0.0,  p0 },   // Particle 1
            //… add more as needed
    };
}

// inline field definition unchanged
inline void get_fields(double t, const double g[],
                       double& Ex, double& Ey, double& Ez,
                       double& Bx, double& By, double& Bz) {
    double theta = ww0*t + k*g[4];
    Ex = E0 * cos(theta);
    Ey = 0.0;
    Ez = 0.0;
    Bx = 0.0;
    By = -Ex;
    Bz = 0.0;
}

// unchanged ODE RHS
int func(double t, const double g[], double dg[], void* /*params*/) {
    double Ex, Ey, Ez, Bx, By, Bz;
    get_fields(t, g, Ex, Ey, Ez, Bx, By, Bz);
//    auto em = Wr(g, t, Ex, Ey, Ez, Bx, By, Bz);

    double vb = sqrt(g[1]*g[1] + g[3]*g[3] + g[5]*g[5]);
    double gr = sqrt(vb*vb + 1.0);

    // dr/dt = c·p/γ
    dg[0] = g[1]*c/gr;
    dg[2] = g[3]*c/gr;
    dg[4] = g[5]*c/gr;
    // dp/dt = Lorentz force
    dg[1] = eq/me/c*(Ex + (g[3]/gr)*Bz - (g[5]/gr)*By);
    dg[3] = eq/me/c*(Ey + (g[5]/gr)*Bx - (g[1]/gr)*Bz);
    dg[5] = eq/me/c*(Ez + (g[1]/gr)*By - (g[3]/gr)*Bx);

    return GSL_SUCCESS;
}

int main() {
    // 1) Get all initial conditions in one shot
    auto initList = initializeParticles();
    int nParticles = initList.size();

    // 4) Open one trace file per particle
    vector<ofstream> traceFile(nParticles);
    for (int pid = 0; pid < nParticles; ++pid) {
        ostringstream fn;
        fn << "/Users/josephzhang/Desktop/Coherence/semiclassical_radition/data/track_pid_" << pid << ".dat";
        traceFile[pid].open(fn.str());
        traceFile[pid] << "# t  x  y  z  vx  vy  vz\n";
    }

    // 2) & 3) GSL solver setup: one solver per particle
    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    vector<gsl_odeiv_step*>    steps(nParticles);
    vector<gsl_odeiv_control*> controls(nParticles);
    vector<gsl_odeiv_evolve*>  evolves(nParticles);
    for (int pid = 0; pid < nParticles; ++pid) {
        steps[pid]    = gsl_odeiv_step_alloc(T, 6);
        controls[pid] = gsl_odeiv_control_y_new(1e-6, 1e-8);
        evolves[pid]  = gsl_odeiv_evolve_alloc(6);
    }
    gsl_odeiv_system sys{func, nullptr, 6, nullptr};

    // Allocate state arrays for each particle
    vector<array<double,6>> g(nParticles);
    vector<double> t_pid(nParticles, get_t0());
    vector<double> h(nParticles, get_dt());
    double dt = get_dt();

    // Initialize each g to its ParticleInit
    for (int pid = 0; pid < nParticles; ++pid) {
        g[pid][0] = initList[pid].x;
        g[pid][1] = initList[pid].px;
        g[pid][2] = initList[pid].y;
        g[pid][3] = initList[pid].py;
        g[pid][4] = initList[pid].z;
        g[pid][5] = initList[pid].pz;
    }

    // Allocate and zero output arrays for coherent sum
    size_t total = (size_t)i_m1 * j_m1 * freq_m1;
    vector<double> TIrx(total,0.0), TIry(total,0.0), TIrz(total,0.0);
    vector<double> TIix(total,0.0), TIiy(total,0.0), TIiz(total,0.0);
    vector<double> TIr_J(total,0.0), TIi_J(total,0.0);

    // Global time loop
    double t = get_t0();
    long long stepCount = 0;
    while (t <= get_tf()) {
        double tNext = t + dt;

        // advance each particle from t->tNext
        for (int pid = 0; pid < nParticles; ++pid) {
            while (t_pid[pid] < tNext) {
                int status = gsl_odeiv_evolve_apply(
                        evolves[pid], controls[pid], steps[pid],
                        &sys,
                        &t_pid[pid], tNext, &h[pid], g[pid].data()
                );
                if (status != GSL_SUCCESS) {
                    cerr << "GSL error pid="<<pid<<" at t="<<t_pid[pid]<< "\n";
                    break;
                }
            }
            // reset solver for next step
            gsl_odeiv_evolve_reset(evolves[pid]);
            gsl_odeiv_step_reset(steps[pid]);
            h[pid] = dt;
        }

        // accumulate integrals — coherent sum over particles first, then ×dt
        for (int i = 0; i < i_m1; ++i) {
            double th_x = x_list[i];
            for (int j = 0; j < j_m1; ++j) {
                double th_y = y_list[j];
                for (int freq_i = 0; freq_i < freq_m1; ++freq_i) {
                    size_t idx = IDX(i,j,freq_i);

                    // sum contributions over all particles
                    double sum_fx_r = 0, sum_fy_r = 0, sum_fz_r = 0;
                    double sum_fx_i = 0, sum_fy_i = 0, sum_fz_i = 0;
                    double sum_Jr = 0,   sum_Ji = 0;

                    for (int pid = 0; pid < nParticles; ++pid) {
                        // current state
                        auto& gv = g[pid];
                        // kinematic quantities
                        double gamma_fac = sqrt(gv[1]*gv[1] + gv[3]*gv[3] + gv[5]*gv[5] + 1.0);
                        double vx = gv[1]/gamma_fac, vy = gv[3]/gamma_fac, vz = gv[5]/gamma_fac;
                        // fields & radiation reaction
                        double Ex,Ey,Ez,Bx,By,Bz;
                        get_fields(tNext, gv.data(), Ex,Ey,Ez,Bx,By,Bz);
                        auto rem = Wr(gv.data(), tNext, Ex,Ey,Ez,Bx,By,Bz);
                        double ax = rem.ax, ay = rem.ay, az = rem.az;

                        double Beta0 = 1.0 - 0.5/pow(gamma_fac, 2);

                        double delta_z = gv[4]/c - Beta0*tNext;
                        double delta_vz = vz - Beta0;

                        double g1 = 1.0/pow(
                                1.0/2.0/el_gamma_ini/el_gamma_ini + (th_x*th_x+th_y*th_y)/2.0 - delta_vz - th_x*vx - th_y*vy,
                                2.0
                        );      // g1 ~ (1 - n dot v)^-2
                        // integrand
                        double fx = g1 * ( th_y*((th_x-vx)*ay - (th_y-vy)*ax)
                                           - (1.0/2.0/el_gamma_ini/el_gamma_ini - (th_x*th_x+th_y*th_y)/2.0 - delta_vz)*ax
                                           + (th_x-vx)*az );
                        double fy = g1 * ( (th_y-vy)*az
                                           - (1.0/2.0/el_gamma_ini/el_gamma_ini - (th_x*th_x+th_y*th_y)/2.0 - delta_vz)*ay
                                           - th_x*((th_x-vx)*ay - (th_y-vy)*ax) );
                        double fz = g1 * ( th_x*((1.0/2.0/el_gamma_ini/el_gamma_ini - (th_x*th_x+th_y*th_y)/2.0 - delta_vz)*ax - (th_x-vx)*az)
                                           - th_y*((th_y-vy)*az
                                                 - (1.0/2.0/el_gamma_ini/el_gamma_ini - (th_x*th_x+th_y*th_y)/2.0 - delta_vz)*ay) );
                        double small_J = g1*(th_x*ax + th_y*ay + az);

                        double k_mu_x_mu = (el_gamma_ini*freq_list[freq_i])/(el_gamma_ini - freq_list[freq_i]) *
                                        ((1.0/2.0/el_gamma_ini/el_gamma_ini + (th_x*th_x+th_y*th_y)/2.0)*tNext - delta_z
                                         - th_x*gv[0]/c - th_y*gv[2]/c);

                        sum_fx_r += fx * cos(k_mu_x_mu);
                        sum_fy_r += fy * cos(k_mu_x_mu);
                        sum_fz_r += fz * cos(k_mu_x_mu);

                        sum_fx_i += fx * sin(k_mu_x_mu);
                        sum_fy_i += fy * sin(k_mu_x_mu);
                        sum_fz_i += fz * sin(k_mu_x_mu);

                        sum_Jr  += small_J * cos(k_mu_x_mu);
                        sum_Ji  += small_J * sin(k_mu_x_mu);
                    } // pid

                    // multiply by dt once
                    TIrx[idx] += sum_fx_r * dt;
                    TIry[idx] += sum_fy_r * dt;
                    TIrz[idx] += sum_fz_r * dt;

                    TIix[idx] += sum_fx_i * dt;
                    TIiy[idx] += sum_fy_i * dt;
                    TIiz[idx] += sum_fz_i * dt;

                    TIr_J[idx] += sum_Jr * dt;
                    TIi_J[idx] += sum_Ji * dt;
                }
            }
        }

        // trajectory output for each particle
        for (int pid = 0; pid < nParticles; ++pid) {
            auto& gv = g[pid];
            double gamma_fac = sqrt(gv[1]*gv[1] + gv[3]*gv[3] + gv[5]*gv[5] + 1.0);
            double vx = gv[1]/gamma_fac, vy = gv[3]/gamma_fac, vz = gv[5]/gamma_fac;
            traceFile[pid]
                    << tNext/tp0 << ' '
                    << gv[0]   << ' '
                    << gv[2]   << ' '
                    << gv[4]   << ' '
                    << vx     << ' '
                    << vy     << ' '
                    << vz     << '\n';
        }

        // 6) Every 200 dt increments print to log (using dt count)
        ++stepCount;
        if (stepCount % 200 == 0) {
            cout << "[LOG] steps=" << stepCount
                 << " t/tp0=" << tNext/tp0
                 << "\n";
        }

        t = tNext;
    }

    // clean up
    for (int pid = 0; pid < nParticles; ++pid) {
        gsl_odeiv_evolve_free(evolves[pid]);
        gsl_odeiv_control_free(controls[pid]);
        gsl_odeiv_step_free(steps[pid]);
        traceFile[pid].close();
    }
    return 0;
}