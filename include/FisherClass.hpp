#pragma once

#include "CosmologyCalculatorClass.hpp"
#include "Helper.hpp"
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
        
        double Cl_derivative(int l, string param_key, double k1, double k2,\
                int *Pk_index, int *Tb_index, int *q_index);
        double Cl_loglog_derivative(int l, string param_key, double k1, double k2,\
                int *Pk_index, int *Tb_index, int *q_index);

        mat Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
                int *Tb_index, int *q_index, vector<double> krange);
            
        double compute_Fl(int l, string param_key1, string param_key2, int ksteps_Cl,\
                double *cond_num, int *Pk_index, int *Tb_index, int *q_index);
        double compute_Fl(int l, string param_key1, string param_key2, double kstepsize,\
                double *cond_num, int *Pk_index, int *Tb_index, int *q_index);

        void initializer(string param_key, int *Pk_index, int *Tb_index, int *q_index);

        double F(string param_key1, string param_key2);
        double F_fixed_kstepsize(string param_key1, string param_key2);
        double F_fixed_kstepsize(int lmin, int lmax, int lsteps);

        void Fl_varying_ksteps(int l, string param_key1, string param_key2, int min_ksteps_Cl,\
                int max_ksteps_Cl, int ksteps_spacing);
        void Fl_varying_ksteps_smart(int l, string param_key1, string param_key2, int min_ksteps_Cl,\
                int max_ksteps_Cl);

        void write_matrix(mat matrix, string filename);
        mat read_matrix(string filename, int n_rows, int n_cols);
        bool check_file(string filename);
        string update_runinfo(int lmin, int lmax,\
                int lstepsize, double kstepsize);
        Fisher_return_pair build_Fisher_inverse(vector<string> filenames_Fl);
        Ellipse find_error_ellipse(Fisher_return_pair fisher, string param1, string param2);


    private:
        vector<double> give_kmodes(int l, double kmax, int steps);
        vector<double> give_kmodes(int l, double k_max, double kstepsize);
        vector<string> model_params_keys;
        CosmoCalc *CALC;
        ofstream Fl_file;
        map<string, double> current_params, fiducial_params, var_params;
        double kmin, kmax;
        vector<double> abcisses_done, logderivs_calculated,\
            abcisses_done_simple, derivs_calculated;
        bool noise, rsd;
};
