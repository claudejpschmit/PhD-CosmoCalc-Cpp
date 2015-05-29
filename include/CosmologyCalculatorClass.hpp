#pragma once

#include "CosmoBasis.hpp"
#include <string>
#include <map>
#include <cmath>
#include <vector>
#include "ClassEngine.hpp"
#include "stdafx.h"
#include "interpolation.h"

using namespace std;
using namespace alglib;

class CosmoCalc : public CosmoBasis {
    
    public:

 
        // ------------ Functions -------------- //
        CosmoCalc(map<string, double> params);

        ~CosmoCalc();

        void show_cosmo_calcs();
        void write_pks(string filename, double z);
        void updateClass(map<string, double> params);

        double hubble_time();
        double hubble_dist();
        double comoving_radial_dist(double z);
        double D_C(double z);
        double D_now(double z);
        double comoving_dist_transverse(double z);
        double D_M(double z);
        double angular_diam_dist(double z, double z2 = -1);
        double D_A(double z, double z2 = -1);
        double luminosity_dist(double z);
        double D_L(double z);
        double comoving_volume(double z);
        double V_C(double z);
        double age_of_universe(double z);
        double light_travel_time(double z);
        double D_ltt(double z);
        double H(double z);
        double H_SI(double z);
        double rho_crit(double z);
        double rho_M(double z);
        double rho_R(double z);
        double rho_V(double z);
        double Omega_M(double z);
        double Omega_R(double z);
        double Omega_V(double z);
        double num_baryons();
        double T(double z);
        double n_H_tot(double z);
        double n_b(double z);
        double n_H(double z);
        double n_p(double z);
        double n_e(double z);

        void update_q();
        void update_Pk_interpolator(map<string, double> params);
        void create_bessel_interpolant(int lmin, int lmax);
        double bessel_j_interp(int l, double x);
        double Pk_interp(double k, double z);
        double Pkz_calc(double k, double z);
        double P_growth(double z);
        double D1(double z);
        double P_delta(double k, string units_k = "default",\
                       string units_P = "default");
        double transfer(double x);

        double Cl(int l, double k1, double k2, double k_low, double k_high);
        double corr_Tb(int l, double k1, double k2, double k_low, double k_high);
        double corr_Tb_rsd(int l, double k1, double k2,\
                           double k_low, double k_high);
        double M(int l, double k1, double k2);
        double delta_Tb_bar(double z);
        double T_S(double z);
        double x_HI(double z);
        double T_K(double z);
        double N_bar(int l, double k1, double k2);

    protected:

        // ------------ Variables -------------- //
        int k_steps, zsteps_Ml, Pk_steps;
        int lmin_bess;
        double zmin_Ml, zmax_Ml, stepsize_Ml, prefactor_Ml;
        vector<double> q_Ml, r_Ml, H_f;
        real_1d_array matterpowerspectrum_k, matterpowerspectrum_z, matterpowerspectrum_P;
        spline2dinterpolant Pk_interpolator;

        vector<spline1dinterpolant> bessel_interp_list;


        ClassParams pars;
        ClassEngine *CLASS;
       

};
