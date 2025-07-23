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
static const double t0 = 0.0;
static const double tf = 0.5 * tp0;    // tp0 == lambda0/c
static const double dt = 0.001 * tp0;

// -----------------------------------------------------------------------------
// 1) ParticleInit holds one particle's editable initial state.
struct ParticleInit {
    double x, y, z;    // initial position
    double px, py, pz; // initial momentum
};

// 1) Define all your particles here.  Directly edit these entries.
//    Add or remove lines to change the number and initial conditions.
vector<ParticleInit> initializeParticles() {
    return {
            //   x     y     z      px    py    pz
            { 0.0,  0.0,  0.0,   0.0,  0.0,  1.0 },   // Particle 0
            { 0.0,  0.0, -0.99*lambda0,   0.0,  0.0,  1.0 },   // Particle 1
            //{ 0.2,  0.0, -0.1,   0.0,  0.0,  1.0 }    // Particle 2
            // … add more as needed
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
    auto em = Wr(g, t, Ex, Ey, Ez, Bx, By, Bz);

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
    // Get all initial conditions in one shot
    auto initList = initializeParticles();
    int nParticles = initList.size();

    // 4) Open one trace file per particle
    vector<ofstream> traceFile(nParticles);
    for (int pid = 0; pid < nParticles; ++pid) {
        ostringstream fn;
        fn << "data/track_pid_" << pid << ".dat";
        traceFile[pid].open(fn.str());
        traceFile[pid] << "# t  x  y  z  vx  vy  vz\n";
    }

    // 2) & 3) GSL solver setup (shared)
    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    gsl_odeiv_step    *step    = gsl_odeiv_step_alloc(T, 6);
    gsl_odeiv_control *control = gsl_odeiv_control_y_new(1e-6, 1e-8);
    gsl_odeiv_evolve  *evolve  = gsl_odeiv_evolve_alloc(6);
    gsl_odeiv_system   sys{func, nullptr, 6, nullptr};


    // Loop over each particle
    for (int pid = 0; pid < nParticles; ++pid) {
        // Load initial state into g[6]
        double g[6] = {
                initList[pid].x,
                initList[pid].px,
                initList[pid].y,
                initList[pid].py,
                initList[pid].z,
                initList[pid].pz
        };

        double t = t0;
        double h = dt;
        int stepCount = 0;

        // Time integration
        while (t <= tf) {
            double tNext = t + dt;
            // advance from t->tNext
            while (t < tNext) {
                int status = gsl_odeiv_evolve_apply(
                        evolve, control, step, &sys, &t, tNext, &h, g
                );
                if (status != GSL_SUCCESS) {
                    cerr << "Error: GSL solve failed at pid="
                         << pid << " t=" << t << "\n";
                    break;
                }
            }
            gsl_odeiv_evolve_reset(evolve);
            gsl_odeiv_step_reset(step);
            h = dt;

            // compute velocity
            double gamma = sqrt(g[1]*g[1] + g[3]*g[3] + g[5]*g[5] + 1.0);
            double vx = g[1] / gamma;
            double vy = g[3] / gamma;
            double vz = g[5] / gamma;

            // write trace line
            traceFile[pid]
                    << t/tp0 << ' '
                    << g[0]   << ' '
                    << g[2]   << ' '
                    << g[4]   << ' '
                    << vx     << ' '
                    << vy     << ' '
                    << vz     << '\n';

            // 6) Every 100 steps print to log
            if (++stepCount % 100 == 0) {
                cout << "[LOG] PID=" << pid
                     << " Step=" << stepCount
                     << " t/tp0=" << fixed << setprecision(3)
                     << (t/tp0)
                     << "\n";
            }
        }
    }

    // clean up
    gsl_odeiv_evolve_free(evolve);
    gsl_odeiv_control_free(control);
    gsl_odeiv_step_free(step);
    for (auto& f : traceFile) f.close();
    return 0;
}