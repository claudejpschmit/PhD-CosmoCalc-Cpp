#pragma once

#include "CosmologyCalculatorClass.hpp"
#include <armadillo>
#include <fstream>
#include <string>
#include "stdafx.h"
#include "interpolation.h"

using namespace alglib;
using namespace arma;

class Fisher {
    public:
        Fisher(map<string, double> params, string Fl_filename);
        ~Fisher();

        void update_Model(map<string, double> new_params, int *Pk_index, int *Tb_index, int *q_index);
        mat compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange);
        mat compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange,\
                spline1dinterpolant bessels);
        mat compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange,\
                spline1dinterpolant bessels, spline1dinterpolant bessels_lminus1);

        double Cl_derivative(int l, string param_key, double k1, double k2,\
                int *Pk_index, int *Tb_index, int *q_index);
        double Cl_loglog_derivative(int l, string param_key, double k1, double k2,\
                int *Pk_index, int *Tb_index, int *q_index);
        vector<vector<double>> Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
                int *Tb_index, int *q_index, vector<double> krange);
        vector<vector<double>> Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
                int *Tb_index, int *q_index, vector<double> krange, spline1dinterpolant bessels);
        vector<vector<double>> Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
                int *Tb_index, int *q_index, vector<double> krange, spline1dinterpolant bessels,\
                spline1dinterpolant bessels_lminus1);

        double compute_Fl(int l, string param_key1, string param_key2,\
                int *Pk_index, int *Tb_index, int *q_index);
        double compute_Fl(int l, string param_key1, string param_key2, int *Pk_index,\
                int *Tb_index, int *q_index, spline1dinterpolant bessels);
        double compute_Fl(int l, string param_key1, string param_key2, int *Pk_index,\
                int *Tb_index, int *q_index, spline1dinterpolant bessels,\
                spline1dinterpolant bessels_lminus1);

        void initializer(string param_key, int *Pk_index, int *Tb_index, int *q_index);

        double F(string param_key1, string param_key2);

    private:
        vector<double> give_kmodes(int l, double kmax, int steps);

        ofstream Fl_file;
        CosmoCalc *CALC;
        map<string, double> current_params, fiducial_params, var_params;
        double kmin, kmax;
        vector<double> abcisses_done, logderivs_calculated,\
            abcisses_done_simple, derivs_calculated;

};
