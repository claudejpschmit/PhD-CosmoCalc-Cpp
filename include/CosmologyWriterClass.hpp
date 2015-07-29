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
        void calculate_comoving(double zmax);
        void calculate_inverse_comoving(double rmax);
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
        void calculate_bessels_cubic(int l);

        void calculate_integrandMM(int l, double k1, double k2, int step);
        void calculate_integrandMN(int l, double k1, double k2, int step);
        void calculate_integrandNN(int l, double k1, double k2, int step);
        
        void calculate_integrandlong(int l, double k1, double k2, int step);
        void calculate_integrandsimple(int l, double k1, double k2, int step);

        void calculate_qdot();
        void calculate_q();

        void generate_movie(int l);
        void calculate_Cl_simple(int l, double k, double k_min, double k_max, double k_stepsize);
        void calculate_Cl_full(int l, double k, double k_min, double k_max, double k_stepsize);
        void calculate_Cl_simple_rsd(int l, double k, double k_min, double k_max, double k_stepsize);
        void calculate_Cl_full_rsd(int l, double k, double k_min, double k_max, double k_stepsize);

        void generate_movie_Cl(int l_min, int l_max, double k, double k_min,\
                double k_max, double k_stepsize);
        void calculate_Ml(int lmin, int lmax, double k, double kappa);
        void calculate_Nl(int lmin, int lmax, double k, double kappa);
};
