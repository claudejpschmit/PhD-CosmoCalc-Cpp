#pragma once

#include <string>
#include <map>
#include <cmath>
#include <boost/math/constants/constants.hpp>

using namespace std;

class CosmoBasis {
    
    public:
        
        // ------------ Functions -------------- //

        /**
         * Constructor that initializes all the parameters for the CosmoCalc backbone.
         *
         * @param params is a dictionary containing the programme parameters. 
         *               The key to the dictionary is a string, and the value a double.
         */
        CosmoBasis(map<string,double> params);

        /**
         * Standard destructor
         */
        ~CosmoBasis();

        /**
         * Tester function to output some parameters for inspection.
         */
        void show_params();
        
        /**
         * This functions return the spherical bessel j as calculated by boost.
         *
         * @param l is the integer index of the spherical bessel function.
         *
         * @param x is the argument of the spherical bessel function.
         */
        double sph_bessel(unsigned int l, double x);


    protected:
        
        // ------------ Functions -------------- //
        
        /**
         * This function verifies whether the fiducial parameters contain
         * all the necessary keys, if not, standard values are added.
         */
        void check_params();
        
        /**
         * This function takes the list of cosmological parameters and translates it 
         * to the parameters which will actually be used in cosmological computations.
         *
         * Note on Neutrinos:
         * Here we add the neutino density to the matter density as well as to the radiation 
         * density. For high redshifts, radiation is dominating so having the neutrinos in
         * has a great effect, since O_nu is of the same order as O_gamma.
         * At late times O_R becomes very unimportant so adding O_nu in has a very minor effect.
         * As for adding O_nu to O_M, O_b/O_nu = 100 and O_CDM/O_nu = 1000 so adding it in
         * or not should not have a great effect at any epoch.
         * 
         * @param params is a dictionary containing the parameters to be translated.
         *
         */
        void generate_params(map<string,double> params);

        /**
         * Helper function to determine the denominator for various integrals, E(z) = H(z)/H_0.
         *
         * @param z redshift at which the function is evaluated.
         */
        double E(double z);

        /**
         * Helper function to integrate 1/E(z) between 0 and z.
         *
         * @param z redshift up to which E^-1 is to be integrated.
         */
        double Z(double z);

        /**
         * This function determines the curvature term in the metric. It is either sinhx, x or sin x.
         *
         * @param x where the curvature term is to be evaluated.
         */
        double S_k(double x);

        /**
         * Helper function to conver Mpc to m
         *
         * @param x distance in Mpc to be converted.
         */
        double mpc_to_m(double x);

        /**
         * Helper function to convert m to Mpc
         *
         * @param x distance in m to be converted.
         */
        double m_to_mpc(double x);

        /**
         * Function which sets the parameters defined in the given object to
         * the newest Planck 2015 results.
         *
         * @param params
         */
        void params_to_planck15(map<string, double> params);
 
        // ------------ Variables -------------- //

        /// \brief params contains the working parameters.
        map<string, double> current_params;
        /**
         * \brief fiducial_params contains the fiducial parameters, 
         *        which are not meant to be changed.
         */
        map<string, double> fiducial_params;
        
        
        /// \brief The CMB temperature today [K].
        double T_CMB;
        /// \brief T_gamma is the photon temperature today [K].
        double T_gamma;
        /// \brief The Hubble constant in standard units today [km s^-1 Mpc^-1].
        double H_0;
        /// \brief The Hubble parameter today.
        double h;
        /// \brief The relative baryon density today.
        double O_b;
        /// \brief The relative Cold Dark Matter density today.
        double O_cdm;
        /// \brief The relative photon density today.
        double O_gamma;
        /// \brief The relative relativistic neutrino density today.
        double O_nu_rel;
        /// \brief The non-relative relativistic neutrino density today.
        double O_nu;
        /// \brief The relative radiation density today.
        double O_R;
        /// \brief The relative curvature density today.
        double O_k;
        /// \brief The relative matter density today.
        double O_M;
        /// \brief The total relative density today.
        double O_tot;
        /// \brief The relative Vacuum density today.
        double O_V;
        /// \brief The hubble distance today [Mpc].
        double D_H;
        /// \brief The hubble time.
        double t_H;
        
        /// \brief Cosmological Bias factor..
        double b_bias;
        /// \brief Scales at horizon crossing [Mpc^-1].
        double k_eq;
        /// \brief T_* = hc/k_b Lambda_21cm [K].
        double T_star;

        // ------------ Constants -------------- //
        
        const double pi = boost::math::constants::pi<double>();
        /// \brief The speed of light in [m/s].
        const double c = 299792458.0;
        /// \brief The Boltzmann Constant in [m^2 kg s^-2 K^-1].
        const double k_b = 1.3806488*pow(10,-23);
        /// \brief The Baryon mass in [kg].
        const double m_b = 1.674927351*pow(10,-27);
        /// \brief The Electron mass in [kg].
        const double m_e = 9.10938291*pow(10,-31);
        /// \brief The Charge of an electron in [Coulomb].
        const double e = 1.60217657*pow(10,-19);
        /// \brief Plancks Constant in [m^2 kg s^-1].
        const double h_planck = 6.62606957*pow(10,-34);
        /// \brief The Gravitational Constant in [m^3 kg^-1 s^2].
        const double G = 6.67384*pow(10,-11);
        
        /// \brief Spontaneous decay rate of the spin flip transition [s^-1].
        const double A_10 = 2.85*pow(10,-15);
        /// \brief Logarithmic tilt of the spectrum of fluctuations.
        const double n_s = 0.95;
        /// \brief Variance of the matter fluctuations today smoothed on a scale of 8 h^-1 Mpc.
        const double sigma_8 = 0.8;

        /// \brief Width of the Reionization regime.
        const double delta_z_rei = 4;
        /// \brief Center of Reionization.
        const double z_rei = 10;
        /// \brief Redshift for Recombination..
        const double z_CMB = 1100;

        /// \brief Beta factor.
        const double beta = 0.7;
        
};
