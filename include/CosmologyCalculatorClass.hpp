#pragma once

#include "CosmoBasis.hpp"
#include <string>
#include <map>
#include <cmath>
#include <vector>
#include "ClassEngine.hpp"
#include "stdafx.h"
#include "interpolation.h"
#include "CAMB_interface.hpp"

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
        void update_Pk_interpolator_direct(map<string, double> params);

        /**
         * Function outputs some standard cosmological calculations to the user.
         */
        void show_cosmo_calcs();
        void updateClass(map<string, double> params);
        void update_q();

        /**
         * Writing the Pk's at a redshift to a file.
         * The Pk's range from k = 0.0001 to k = 10.
         * Two files are writte, the original filename contains CLASS values 
         * and a second file, with the prefix "new_" contains the interpolated 
         * values from CAMB.
         *
         * @param filename is the name of the file to which the 
         *                 results will be written.
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
         * @param z is the redshift up to which the distance is calculated.
         */
        double comoving_radial_dist(double z);

        /**
         * Alternative for comoving_radial_dist.
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double D_C(double z);
        
        /**
         * Alternative for comoving_radial_dist.
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double D_now(double z);

        /**
         * Determines the comoving distance (transverse) in [Mpc],
         * aka. Proper motion distance:
         *  D_M = D_H/sqrt(1-O_tot) S_k(sqrt(1-O_tot) D_C/D_H)
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double comoving_dist_transverse(double z);
        
        /**
         * Alternative for comoving_dist_transverse.
         *
         * @param z is the redshift up to which the distance is calculated.
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
         * @param z is the redshift up to which the distance is calculated.
         *
         * @param z2 is the redshift of the second object. This is optional.
         */
        double angular_diam_dist(double z, double z2 = -1);

        /**
         * Alternative for angular_diam_dist.
         *
         * @param z is the redshift up to which the distance is calculated.
         *
         * @param z2 is the redshift of the second object. This is optional
         */
        double D_A(double z, double z2 = -1);

        /**
         * Determines the luminosity distance [Mpc]:
         *  D_L = (1+z)^2 * D_A
         *
         * @param z redshift of the object.
         */
        double luminosity_dist(double z);

        /**
         * Alternative for luminosity_dist [Mpc].
         *
         * @param z redshift of the object.
         */
        double D_L(double z);

        /**
         * Determines the comoving volume up to a redshift in [Mpc^3]:
         *  V_C = int D_H (1+z)^2 * D_A^2 / E(z) dOmega dz
         *
         * @param z is the redshift up to which the volume is calculated.
         */
        double comoving_volume(double z);
        
        /**
         * Alternative for comoving_volume [Mpc^3].
         *
         * @param z is the redshift up to which the volume is calculated.
         */
        double V_C(double z);

        /**
         * Determines the age of the universe at the redshift z [s * Mpc/km]:
         *  t(z) = t_H int(dz/((1+z) E(z)), z, inf)
         *
         * @param z redshift at which the age is calculated.
         */
        double age_of_universe(double z);

        /**
         * Determines the light travel time [s * Mpc/km]:
         *  ltt = t(0) - t(z)
         *
         * @param z is the redshift up to which the light travel time 
         *           is calculated.
         */
        double light_travel_time(double z);

        /** 
         * Determines the distance based on the light travel time [Mpc]:
         *  D_ltt = c * (t(0) - t(z))
         *
         * @param z is the redshift up to which the distance is calculated.
         */
        double D_ltt(double z);

        /**
         * Determines the Hubble Constant at a redshift [km * s^-1 * Mpc^-1]
         *
         * @param z is the redshift at which the Hubble constant is calculated.
         */
        double H(double z);

        /**
         * Determines the Hubble Constant in SI units [s^-1]
         *
         * @param z is the redshift at which the Hubble constant is calculated.
         */
        double H_SI(double z);

        /** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */
        double rho_crit(double z);
        
        /** 
         * Determines the matter density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the matter density is calculated.
         */
        double rho_M(double z);
        
        /** 
         * Determines the radiation density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the radiation density is calculated.
         */
        double rho_R(double z);
        
        /** 
         * Determines the vacuum density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the vacuum density is calculated.
         */
        double rho_V(double z);
        
        /** 
         * Determines the relative matter density at a redshift.
         *
         * @param z is the redshift at which the relative matter density 
         *          is calculated.
         */
        double Omega_M(double z);
        
        /** 
         * Determines the relative radiation density at a redshift.
         *
         * @param z is the redshift at which the relative radiation  density 
         *          is calculated.
         */

        double Omega_R(double z);
        
        /** 
         * Determines the relative vacuum density at a redshift.
         *
         * @param z is the redshift at which the relative vacuum  density 
         *          is calculated.
         */
        double Omega_V(double z);
        
        /** 
         * Determines an estimate for the number of baryons in the Universe.
         */
        double num_baryons();
        
        /** 
         * Determines the (CMB) temperature of the universe at a redshift [K].
         *
         * @param z is the redshift at which the temperature is calculated.
         */

        double T(double z);
        
        /** 
         * Determines the total hydrogen density at a redshift [m^-3].
         *
         * @param z is the redshift at which the hydrogen density is calculated.
         */
        double n_H_tot(double z);
        
        /** 
         * Determines the baryon number density at a redshift [m^-3].
         *
         * @param z is the redshift at which the baryon density is calculated.
         */
        double n_b(double z);
        
        /** 
         * Determines the hydrogen density at a redshift [m^-3].
         *
         * @param z is the redshift at which the hydrogen density is calculated.
         */
        double n_H(double z);
        
        /** 
         * Determines the (free) proton density at a redshift [m^-3].
         *
         * @param z is the redshift at which the proton density is calculated.
         */
        double n_p(double z);
        
        /** 
         * Determines the (free) electron density at a redshift [m^-3].
         *
         * @param z is the redshift at which the electron density is calculated.
         */
        double n_e(double z);
        
        /** 
         * Determines the spherical bessel function j_l(x) from the CAMB 
         * interpolation.
         *
         * @param l is the index of the spherical bessel function.
         *
         * @param x is the value at which the bessel function is evaluated.
         */
        double bessel_j_interp(int l, double x);
        double bessel_j_interp_basic(int l, double x);
        double bessel_j_interp_cubic(int l, double x);


        /** 
         * Determines the matter power spectrum P(k,z) from the CAMB 
         * interpolation.
         *
         * @param k is the scale at whcih the power spectrum is evaluated.
         *
         * @param z is the redshift at which the power spectrum is evaluated.
         */
        double Pk_interp(double k, double z);

        /** 
         * TODO: determine units...
         * Determines the matter power spectrum P(k,z) from a basic analytic
         * calculation.
         *
         * @param k is the scale at which the power spectrum is evaluated.
         *
         * @param z is the redshift at which the power spectrum is evaluated.
         */
        double Pkz_calc(double k, double z);
        
        /** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */
        double P_growth(double z);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double D1(double z);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double P_delta(double k, string units_k = "default",\
                       string units_P = "default");/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double transfer(double x);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */


        double Cl(int l, double k1, double k2, double k_low, double k_high);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double corr_Tb(int l, double k1, double k2, double k_low, double k_high);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double corr_Tb_rsd(int l, double k1, double k2,\
                           double k_low, double k_high);
        
        /** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double M(int l, double k1, double k2);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double delta_Tb_bar(double z);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double T_S(double z);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double x_HI(double z);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double T_K(double z);/** 
         * Determines the critical density at a redshift [kg/m^3].
         *
         * @param z is the redshift at which the critical density is calculated.
         */

        double N_bar(int l, double k1, double k2);

        double integrandMM(int l, double k1, double k2, double k);
        double integrandMN(int l, double k1, double k2, double k);
        double integrandNN(int l, double k1, double k2, double k);

    protected:

        // ------------ Functions -------------- //
        
        void create_bessel_interpolant_ALGLIB(int lmin, int lmax);
        void create_bessel_interpolant_OWN(int lmax);
        // ------------ Variables -------------- //
        int zsteps_Ml, Pk_steps;
        int lmin_bess;
        double zmin_Ml, zmax_Ml, stepsize_Ml, prefactor_Ml, k_stepsize;
        vector<double> q_Ml, r_Ml, H_f;
        spline2dinterpolant Pk_interpolator;

        vector<spline1dinterpolant> bessel_interp_list;
        
        vector<vector<double>> bessel_values;

        ClassParams pars;
        ClassEngine *CLASS;
        CAMB_CALLER *CAMB;
       
};
