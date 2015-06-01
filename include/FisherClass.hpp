#pragma once

#include "CosmologyCalculatorClass.hpp"
#include <armadillo>

using namespace arma;

class Fisher {
    public:
        Fisher(map<string, double> params);
        ~Fisher();

        void update_Model(map<string, double> new_params);
        void compute_Cl(int l);
        void compute_Cl_inv();
        double Cl_derivative(int l, string param_key, double k1, double k2);
        double Cl_loglog_derivative(int l, string param_key,\
                double k1, double k2);
        vector<vector<double>> Cl_derivative_matrix(int l, string param_key);
        double compute_Fl(int l, string param_key1, string param_key2);
        double F(string param_key1, string param_key2);
        void write_logder(string param_key, double param_val,\ 
                double stepsize_low, double stepsize_high,\
                int steps, int l, double k1, double k2,\
                string suffix = "");


    private:
        CosmoCalc *CALC;
        map<string, double> current_params, fiducial_params, var_params;
        mat Cl, Cl_inv;
        vector<double> krange;
        double kmin, kmax;
        vector<double> abcisses_done, logderivs_calculated,\
            abcisses_done_simple, derivs_calculated;
};
