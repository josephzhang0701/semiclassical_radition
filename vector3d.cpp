//******************************************************************************
//** SCATMECH: Polarized Light Scattering C++ Class Library
//**
//** File: vector.cpp
//**
//** Thomas A. Germer
//** Sensor Science Division, National Institute of Standards and Technology
//** 100 Bureau Dr. Stop 8443; Gaithersburg, MD 20899-8443
//** Phone: (301) 975-2876
//** Email: thomas.germer@nist.gov
//**
//** Version: 7.00 (January 2015)
//**
//******************************************************************************
#include <iostream>
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
   #include "vector3d.h"

   using namespace std;
namespace SCATMECH {
    Vector3D<double>
    perpto(const Vector3D<double>& a,const Vector3D<double>& b)
    {
        // Modified by TAG 14 MAR 2005 by adding these two lines, and comparing
        // _norm to a small number, rather than zero...
        Vector3D<double> _a = unit(a);
        Vector3D<double> _b = unit(b);

        Vector3D<double> temp=cross(_a,_b);
        double _norm=Norm(temp);
        if (_norm>1E-10) return temp/_norm;
        else {
            static Vector3D<double> x(1.,0.,0.);
            Vector3D<double> temp=cross(_a,x);
            double _norm = Norm(temp);
            if (_norm>1E-10) return temp/_norm;
            else {
                static Vector3D<double> y(0.,1.,0.);
                Vector3D<double> temp=cross(_a,y);
                double _norm = Norm(temp);

                return temp/_norm;
            }
        }
    }


} // namespace SCATMECH




