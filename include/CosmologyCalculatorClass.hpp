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

        /**
         * Constructor that initializes all the parameters 
         * for the CosmoCalc working class.
         *
         * @param params is a dictionary containing the programme parameters. 
         *               The key to the dictionary is a string, 
         *               and the value a double.
         */
        CosmoCalc(map<string, double> params);

        /**
         * Standard destructor, destroying CLASS object used.
         */
        ~CosmoCalc();

        void update_Pk_interpolator(map<string, double> params);

        /**
         * Function outputs some standard cosmological calculations to the user.
         */
        void show_cosmo_calcs();
        void updateClass(map<string, double> params);
        void update_q();

        /**
         * Writing the Pk's at a redshift to a file.
         * The Pk's range from k = 0.0001 to k = 10.
         * Two files are writte, the original filename contains CLASS values and a second file,
         * with the prefix "new_" contains the interpolated values from CAMB.
         *
         * @param filename is the name of the file to which the results will be written.
         *
         * @param z is the redshift at which the power spectrum is calculated.
         */
        void write_pks(string filename, double z);
        
        /**
         * Determines the Hubble Time in [s * Mpc/km]:
         *  t_H = 1/H_0
         */
        double hubble_time();

        /**
         * Determines the Hubble Distance in [Mpc]:
         *  D_H = c/H_0
         */
        double hubble_dist();
        
        /**
         * Determines the comoving distance (line of sight) in [Mpc],
         * aka. Comoving radial distance, D_now:
         *  D_C = D_H int(dz / E(z), 0, z) = D_H * Z(z) = cZ(z)/H_0
         *
         * @params z is the redshift up to which the distance is calculated.
         */
        double comoving_radial_dist(double z);

        /**
         * Alternative for comoving_radial_dist.
         *
         * @params z is the redshift up to which the distance is calculated.
         */
        double D_C(double z);
        
        /**
         * Alternative for comoving_radial_dist.
         *
         * @params z is the redshift up to which the distance is calculated.
         */
        double D_now(double z);

        /**
         * Determines the comoving distance (transverse) in [Mpc],
         * aka. Proper motion distance:
         *  D_M = D_H/sqrt(1-O_tot) S_k(sqrt(1-O_tot) D_C/D_H)
         *
         * @params z is the redshift up to which the distance is calculated.
         */
        double comoving_dist_transverse(double z);
        
        /**
         * Alternative for comoving_dist_transverse.
         *
         * @params z is the redshift up to which the distance is calculated.
         */
        double D_M(double z);

        /**
         * Determines the angular diameter distance [Mpc]:
         *  D_A = D_M / (1+z)
         * Second form:
         * Angular diameter distance between 2 objects at redshifts z & z2 [Mpc],
         * formula only holds for O_tot <= 1:
         *  D_A12 = 1/(1+z2) * (D_M2 sqrt(1+(1-O_tot) D_M1^2/D_H^2) -
         *                      H_M1 sqrt(1+(1-O_tot) D_M2^2/D_H^2))
         *
         * @params z is the redshift up to which the distance is calculated.
         *
         * @params z2 is the redshift of the second object. This is optional.
         */
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

        // ------------ Functions -------------- //
        
        void create_bessel_interpolant(int lmin, int lmax);

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
