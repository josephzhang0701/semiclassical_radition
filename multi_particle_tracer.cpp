// Created by Joseph Zhang on 7/23/25.
// Refactored for parallel performance using OpenMP; original comments preserved.

#include "const.h"

#include <iostream>
#include <cstdio>
#include <cmath>
#include <memory>
#include <vector>
#include <array>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <omp.h>         // Added for OpenMP

#include "linespace.h"
#include "vector3d.h"
#include "chi_e.h"

using namespace std;
using namespace SCATMECH;

// =================== 文件名后缀定义 ===================
const std::string output_suffix = "dual_thin_";  // 可根据需要设为任何字符串，比如"_caseA"

// =================== Verbose Level 控制 ===================
// 0: No output except fatal errors
// 1: Main steps only
// 2: + Each i
// 3: + Each (i,j)
// 4: + Each (i,j,freq_i)  (非常多！仅用于极小网格debug)
const int verbose = 1;

// ========== Simulation Constants ==========
inline double get_t0() { return 0.0; }
inline double get_dt() { return 0.001 * tp0; }
inline double get_tf() { return 7.5   * tp0; }

constexpr int i_m1 = 300, j_m1 = 300, freq_m1 = 300;
constexpr size_t N_GRID = static_cast<size_t>(i_m1) * j_m1 * freq_m1;

// ========== Grid & Physics Parameters ==========
static double el_gamma_ini = (Ekinmean / me_c2) + 1.0;
static const auto x_list = linespace(-7.5/el_gamma_ini, 7.5/el_gamma_ini, i_m1).generate();
static const auto y_list = linespace(-7.5/el_gamma_ini, 7.5/el_gamma_ini, j_m1).generate();
static const auto freq_list = linespace(
        1e-7*el_gamma_ini * me_c2/hb,
        (1-1e-7)*el_gamma_ini * me_c2/hb,
        freq_m1
).generate();

static const double dtx = x_list[1] - x_list[0];
static const double dty = y_list[1] - y_list[0];
static const double dw  = freq_list[1] - freq_list[0];

// 3D index helper\inline
size_t IDX(int i, int j, int the_k) {
    return (static_cast<size_t>(i) * j_m1 + j) * freq_m1 + the_k;
}

// ========== Particle Initialization ==========
struct ParticleInit {
    double x, y, z;    // initial position
    double px, py, pz; // initial momentum
};

// 先声明一个空的全局列表
static std::vector<ParticleInit> particle_init_list;

// 静态构造器，在 main 之前就把两个粒子 push 进去
struct ParticleInitLoader {
    ParticleInitLoader() {
        double p0 = std::sqrt(el_gamma_ini * el_gamma_ini - 1.0);
        particle_init_list.reserve(2);
        particle_init_list.push_back({ 0.0,  0.0,  0.0,           0.0, 0.0, p0 });
        particle_init_list.push_back({ 0.0,  0.0, -0.99*lambda0,    0.0, 0.0, p0 });
    }
} particle_init_loader;


// ========== Field Calculation ==========
inline void get_fields(double t, const double g[6],
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

// ========== ODE RHS ==========
int func(double t, const double g[], double dg[], void*) {
    double Ex, Ey, Ez, Bx, By, Bz;
    get_fields(t, g, Ex, Ey, Ez, Bx, By, Bz);

    double px=g[1], py=g[3], pz=g[5];
    double gr = sqrt(px*px + py*py + pz*pz + 1.0);

    double inv_gr = 1.0/gr, coef = eq/(me*c);

    dg[0] = px * c * inv_gr;
    dg[2] = py * c * inv_gr;
    dg[4] = pz * c * inv_gr;

    dg[1] = coef*(Ex + (py*inv_gr)*Bz - (pz*inv_gr)*By);
    dg[3] = coef*(Ey + (pz*inv_gr)*Bx - (px*inv_gr)*Bz);
    dg[5] = coef*(Ez + (px*inv_gr)*By - (py*inv_gr)*Bx);

    return GSL_SUCCESS;
}

// ========== Main ==========
int main() {

    // --- 1. Particle Initial State ---
    std::cout << "[LOG] Initializing particle states ..." << std::endl;
    const auto& plist = particle_init_list;

    // --- 2. Open Trace Files using FILE* ---
    int N_PARTICLES = plist.size();
    std::vector<FILE*> traceFile(N_PARTICLES);

    std::cout << "[LOG] Opening particle trace files ..." << std::endl;
    for (int pid = 0; pid < N_PARTICLES; ++pid) {
        std::string filename = "/Users/josephzhang/Desktop/Coherence/semiclassical_radition/data/"
                               + output_suffix
                               + "track_pid_"
                               + std::to_string(pid)
                               + ".dat";
        traceFile[pid] = fopen(filename.c_str(), "w");
        if (!traceFile[pid]) {
            std::cerr << "Error: cannot open " << filename << " for writing\n";
            std::exit(1);
        }
        fputs("# t  x  y  z  vx  vy  vz\n", traceFile[pid]);
    }

    // --- 3. GSL ODE Solvers ---
    std::cout << "[LOG] Setting up GSL ODE solvers ..." << std::endl;
    const gsl_odeiv_step_type* T = gsl_odeiv_step_rkf45;
    std::vector<gsl_odeiv_step*>    steps(N_PARTICLES);
    std::vector<gsl_odeiv_control*> controls(N_PARTICLES);
    std::vector<gsl_odeiv_evolve*>  evolves(N_PARTICLES);
    for (int pid = 0; pid < N_PARTICLES; ++pid) {
        steps[pid]    = gsl_odeiv_step_alloc(T, 6);
        controls[pid] = gsl_odeiv_control_y_new(1e-6, 1e-8);
        evolves[pid]  = gsl_odeiv_evolve_alloc(6);
    }
    gsl_odeiv_system sys{func, nullptr, 6, nullptr};

    // --- 4. Particle State Arrays ---
    std::vector<std::array<double,6>> g(N_PARTICLES);
    std::vector<double> t_pid(N_PARTICLES, get_t0());
    std::vector<double> h    (N_PARTICLES, get_dt());
    for (int pid = 0; pid < N_PARTICLES; ++pid) {
        // 注意顺序：{ x, px, y, py, z, pz }
        g[pid] = {
                plist[pid].x,  plist[pid].px,
                plist[pid].y,  plist[pid].py,
                plist[pid].z,  plist[pid].pz
        };
        std::cerr << "[DBG] g["<<pid<<"] = { "
                  << g[pid][0]<<", "
                  << g[pid][1]<<", "
                  << g[pid][2]<<", "
                  << g[pid][3]<<", "
                  << g[pid][4]<<", "
                  << g[pid][5]
                  << " }\n";
    }
    double dt = get_dt();

    // --- 5. Allocate Output Arrays ---
    // ===== memory allocate with alignment & restrict =====
    double * __restrict TIrx, * __restrict TIry, * __restrict TIrz;
    double * __restrict TIix, * __restrict TIiy, * __restrict TIiz;
    double * __restrict TJr, * __restrict TJi;
    auto alloc_aligned = [&](size_t n){
        void *p;
        posix_memalign(&p, 64, n * sizeof(double));
        memset(p, 0, n * sizeof(double));
        return (double*)p;
    };
    TIrx   = alloc_aligned(N_GRID);
    TIry   = alloc_aligned(N_GRID);
    TIrz   = alloc_aligned(N_GRID);
    TIix   = alloc_aligned(N_GRID);
    TIiy   = alloc_aligned(N_GRID);
    TIiz   = alloc_aligned(N_GRID);
    TJr  = alloc_aligned(N_GRID);
    TJi  = alloc_aligned(N_GRID);

    // --- 6. Main Simulation Loop ---
    std::cout << "[LOG] Starting main simulation loop ..." << std::endl;

    double t = get_t0();
    int stepCount = 0;

    omp_set_num_threads(omp_get_max_threads());

    while (t <= get_tf()) {
        double tNext = t + dt;

        // --- 6.1 Advance Each Particle ---
        for (int pid = 0; pid < N_PARTICLES; ++pid) {
            while (t_pid[pid] < tNext) {
                int status = gsl_odeiv_evolve_apply(
                        evolves[pid], controls[pid], steps[pid],
                        &sys, &t_pid[pid], tNext, &h[pid], g[pid].data());
                if (status != GSL_SUCCESS) break;
            }
            gsl_odeiv_evolve_reset(evolves[pid]);
            gsl_odeiv_step_reset(steps[pid]);
            h[pid] = dt;
        }

        // --- 6.2 Accumulate Integrals on Coherent Grid (Parallelized) ---
        #pragma omp parallel for collapse(3) schedule(dynamic)
        for (int i = 0; i < i_m1; ++i) {
            double th_x = x_list[i];
            for (int j = 0; j < j_m1; ++j) {
                double th_y = y_list[j];
                for (int freq_i = 0; freq_i < freq_m1; ++freq_i) {
                    size_t idx = IDX(i,j,freq_i);
                    double sum_fx_r=0, sum_fy_r=0, sum_fz_r=0;
                    double sum_fx_i=0, sum_fy_i=0, sum_fz_i=0;
                    double sum_Jr=0, sum_Ji=0;
                    double the_freq = freq_list[freq_i];


                    for (int pid = 0; pid < N_PARTICLES; ++pid) {
                        double* gv = g[pid].data();
                        double px=gv[1], py=gv[3], pz=gv[5];
                        double gamma_fac = sqrt(px*px + py*py + pz*pz + 1.0);
                        double vx = px/gamma_fac, vy = py/gamma_fac, vz = pz/gamma_fac;

                        double Ex,Ey,Ez,Bx,By,Bz;
                        get_fields(tNext, gv, Ex,Ey,Ez,Bx,By,Bz);
                        auto rem = Wr(gv, tNext, Ex,Ey,Ez,Bx,By,Bz);
                        double ax = rem.ax, ay = rem.ay, az = rem.az;

                        double Beta0 = 1.0 - 0.5/(gamma_fac*gamma_fac);
                        double delta_z = gv[4]/c - Beta0*tNext;
                        double delta_vz = vz - Beta0;

                        double tmp = 0.5/(el_gamma_ini*el_gamma_ini) +
                                     0.5*(th_x*th_x+th_y*th_y) -
                                     delta_vz - th_x*vx - th_y*vy;
                        double inv_tmp2 = 1.0/(tmp*tmp);

                        double fx = inv_tmp2 * ( th_y*((th_x-vx)*ay - (th_y-vy)*ax)
                                                 - (0.5/(el_gamma_ini*el_gamma_ini) - 0.5*(th_x*th_x+th_y*th_y) - delta_vz)*ax
                                                 + (th_x-vx)*az );
                        double fy = inv_tmp2 * ( (th_y-vy)*az
                                                 - (0.5/(el_gamma_ini*el_gamma_ini) - 0.5*(th_x*th_x+th_y*th_y) - delta_vz)*ay
                                                 - th_x*((th_x-vx)*ay - (th_y-vy)*ax) );
                        double fz = inv_tmp2 * ( th_x*((0.5/(el_gamma_ini*el_gamma_ini) - 0.5*(th_x*th_x+th_y*th_y) - delta_vz)*ax - (th_x-vx)*az)
                                                 - th_y*((th_y-vy)*az
                                                         - (0.5/(el_gamma_ini*el_gamma_ini) - 0.5*(th_x*th_x+th_y*th_y) - delta_vz)*ay) );
                        double small_J = inv_tmp2*(th_x*ax + th_y*ay + az);

                        double k_mu_x_mu = (el_gamma_ini * the_freq)/(el_gamma_ini - the_freq*hb/me_c2) *
                                           ((0.5/(el_gamma_ini*el_gamma_ini) + 0.5*(th_x*th_x+th_y*th_y))*tNext
                                            - delta_z - th_x*gv[0]/c - th_y*gv[2]/c);

                        double cs = cos(k_mu_x_mu), sn = sin(k_mu_x_mu);

                        sum_fx_r += fx * cs;
                        sum_fy_r += fy * cs;
                        sum_fz_r += fz * cs;
                        sum_fx_i += fx * sn;
                        sum_fy_i += fy * sn;
                        sum_fz_i += fz * sn;
                        sum_Jr   += small_J * cs;
                        sum_Ji   += small_J * sn;
                    }

                    TIrx[idx] += sum_fx_r * dt;
                    TIry[idx] += sum_fy_r * dt;
                    TIrz[idx] += sum_fz_r * dt;
                    TIix[idx] += sum_fx_i * dt;
                    TIiy[idx] += sum_fy_i * dt;
                    TIiz[idx] += sum_fz_i * dt;
                    TJr[idx] += sum_Jr * dt;
                    TJi[idx] += sum_Ji * dt;
                }
            }
        }

        // --- 6.3 Particle Trajectory Output ---
        ++stepCount;
        if (stepCount % 200 == 1) {
            std::cout << "[LOG] steps=" << stepCount
                      << " t/tp0=" << tNext / tp0
                      << std::endl;
        }

        for (int pid = 0; pid < N_PARTICLES; ++pid) {
            double* gv = g[pid].data();
            double gamma_fac = sqrt(gv[1]*gv[1] + gv[3]*gv[3] + gv[5]*gv[5] + 1.0);
            double vx = gv[1]/gamma_fac, vy = gv[3]/gamma_fac, vz = gv[5]/gamma_fac;
            fprintf(traceFile[pid], "%f %e %e %e %e %e %e\n",
                    tNext/tp0, gv[0], gv[2], gv[4], vx, vy, vz);
        }

        t = tNext;
    }

    // --- 7. Clean Up ---
    for (int pid = 0; pid < N_PARTICLES; ++pid) {
        gsl_odeiv_evolve_free(evolves[pid]);
        gsl_odeiv_control_free(controls[pid]);
        gsl_odeiv_step_free(steps[pid]);
        fclose(traceFile[pid]);
    }

    // =================== 7.1 能谱 (Spectra) ===================
    std::cout << "[LOG] Calculating spectra ..." << std::endl;
    std::string spectra_filename = "/Users/josephzhang/Desktop/Coherence/semiclassical_radition/data/" + output_suffix + "spectra.dat";
    FILE* fout1 = fopen(spectra_filename.c_str(), "w");
    if (!fout1) {
        std::cerr << "Error: Cannot open spectra.dat for writing!" << std::endl;
        return 1;
    }
    for (int ki = 0; ki < freq_m1; ++ki) {
        double PP = 0.0;
//        #pragma omp parallel for collapse(2) schedule(static)
        for (int i = 0; i < i_m1; ++i) {
            for (int j = 0; j < j_m1; ++j) {
                size_t idx = IDX(i, j, ki);
                double rr02 = el_gamma_ini - hb * freq_list[ki] / me_c2;
                double C = eq * eq / 4.0 / pi / pi / hb;
                double D1 = (rr02*rr02 + el_gamma_ini*el_gamma_ini) / (2.0 * el_gamma_ini*el_gamma_ini);
                double D2 = (hb * freq_list[ki] / me_c2);
                D2 = D2*D2 / (2.0 * el_gamma_ini*el_gamma_ini*el_gamma_ini*el_gamma_ini);
                Vector3D<double> Ireal{TIrx[idx], TIry[idx], TIrz[idx]};
                Vector3D<double> Iimag{TIix[idx], TIiy[idx], TIiz[idx]};
                double Jreal = TJr[idx], Jimag = TJi[idx];
                PP += C * (D1 * (dot(Ireal, Ireal) + dot(Iimag, Iimag)) + D2 * (Jreal*Jreal + Jimag*Jimag)) * dtx * dty;
            }
        }
        fprintf(fout1, "%.15e %.15e\n", hb * freq_list[ki] / me_c2 / el_gamma_ini, PP);
    }
    fclose(fout1);
    std::cout << "[LOG] Spectra calculation finished." << std::endl;

    // =================== 7.2 角分布 (Angular Distribution) ===================
    std::cout << "[LOG] Calculating angular distribution ..." << std::endl;
    std::string angdist_filename = "/Users/josephzhang/Desktop/Coherence/semiclassical_radition/data/" + output_suffix + "angular_dist.dat";
    FILE* fout2 = fopen(angdist_filename.c_str(), "w");
    if (!fout2) {
        std::cerr << "Error: Cannot open angular_dist.dat for writing!" << std::endl;
        return 1;
    }
//    #pragma omp parallel for collapse(2) schedule(dynamic)
    for (int i = 0; i < i_m1; ++i) {
        for (int j = 0; j < j_m1; ++j) {
            double PP = 0.0;
            for (int ki = 0; ki < freq_m1; ++ki) {
                size_t idx = IDX(i, j, ki);
                double rr02 = el_gamma_ini - hb * freq_list[ki] / me_c2;
                double C = eq * eq / 4.0 / pi / pi;
                double D1 = (rr02*rr02 + el_gamma_ini*el_gamma_ini) / (2.0 * el_gamma_ini*el_gamma_ini);
                double D2 = (hb * freq_list[ki] / me_c2);
                D2 = D2*D2 / (2.0 * el_gamma_ini*el_gamma_ini*el_gamma_ini*el_gamma_ini);
                Vector3D<double> Ireal{TIrx[idx], TIry[idx], TIrz[idx]};
                Vector3D<double> Iimag{TIix[idx], TIiy[idx], TIiz[idx]};
                double Jreal = TJr[idx], Jimag = TJi[idx];
                PP += C * (D1 * (dot(Ireal, Ireal) + dot(Iimag, Iimag)) + D2 * (Jreal*Jreal + Jimag*Jimag)) * dw;
            }
            fprintf(fout2, "%.8f %.8f %.15e\n", x_list[i], y_list[j], PP);
        }
    }
    fclose(fout2);
    std::cout << "[LOG] Angular distribution calculation finished." << std::endl;

    return 0;
}
