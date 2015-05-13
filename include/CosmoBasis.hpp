#pragma once

#include <string>
#include <map>
#include <cmath>

using namespace std;

class CosmoBasis {
    
    public:
        /**
         * Constructor that initializes all the parameters for the CosmoCalc backbone.
         *
         * @param params is a dictionary containing the programme parameters. 
         *               The key to the dictionary is a string, and the value a double.
         */
        CosmoBasis(map<char*,double> params);

        void show_params();

    protected:

        /// \brief params contains the working parameters.
        map<char*, double> params;
        /**
         * \brief fiducial_params contains the fiducial parameters, 
         *        which are not meant to be changed.
         */
        map<char*, double> fiducial_params;
        /// \brief The speed of light in [m/s].
        double c;
        /// \brief The Boltzmann Constant in [m^2 kg s^-2 K^-1].
        double k_b;
        /// \brief The Baryon mass in [kg].
        double m_b;
        /// \brief The Electron mass in [kg].
        double m_e;
        /// \brief The Charge of an electron in [Coulomb].
        double e;
        /// \brief Plancks Constant in [m^2 kg s^-1].
        double h_planck;
        /// \brief The Gravitational Constant in [m^3 kg^-1 s^2].
        double G;
        /// \brief T_* = hc/k_b Lambda_21cm [K].
        double T_star;

};
