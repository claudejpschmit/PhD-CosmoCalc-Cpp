#pragma once

#include "CosmologyCalculatorClass.hpp"
#include <iostream>
#include <fstream>

using namespace std;

class CosmoWrite : public CosmoCalc {
    public:
 
        // ------------ Functions -------------- //
        CosmoWrite(map<string, double> params);

        ~CosmoWrite();

        void calculate_Ml(int l, double k_fixed, double k2_low, double k2_high,\
                          int step, double stepsize = 0);
        void calculate_Nl(int l, double k_fixed, double k2_low, double k2_high,\
                          int step, double stepsize = 0);
        void calculate_distances(double zmax);
        void calculate_densities_rho(double zmax);
        void calculate_densities_Omega(double zmax);
        void calculate_H(double zmax);
        void calculate_P(double k_min, double k_max, int step,\
                         string unit_k, string unit_P);
        void calculate_P_CLASS(double kmin, double kmax, double z, int step);
        void calculate_dTb(double zmin, double zmax, int step);
        void calculate_xHI(double zmin, double zmax, int step);
        void calculate_Ts(double zmin, double zmax, int step);
        void calculate_Tk(double zmin, double zmax, int step);
        void calculate_P_compare(double k_low, double k_high, int kstep, double z_low, double z_high, int zstep);
        void calculate_bessels(int l);
        void calculate_bessels_basic(int l);
        void calculate_bessels_exact(int l);

};
