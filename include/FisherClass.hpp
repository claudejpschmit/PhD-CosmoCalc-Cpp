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
        
        /**
         * Constructor for the Fisher object. 
         * This constructor initializes the following objects:
         *      - the CALC pointer, which initializes the CosmoCalc object to 
         *        be used to calculate the Cl's.
         *      - the fiducial parameters of the simulation.
         *      - the parameters that are to be varied by the run
         *      - stores information about whether RSD and Noise are to be 
         *        calculated.
         * @param params is a dictionary that determines the initial/fiducial
         *          values of the parameters of the run.
         * @param Fl_filename determines the file into which the sum over all
         *          the Fl's are to be written to. TODO: Check whether this is 
         *          necessary or not.
         * @param params_keys_considered is a list with the parameter keys that
         *          are to be varied in the run.
         */
        Fisher(map<string, double> params, string Fl_filename,\
                vector<string> params_keys_considered);

        /**
         * The destructor releases the CALC pointer as well as closing the
         * Fl_file.
         */
        ~Fisher();
        
        /**
         *
         */
        double F(string param_key1, string param_key2);
        
        /**
         *
         */
        double F_fixed_kstepsize(string param_key1, string param_key2);
        
        /**
         * This function determines the Fl's in some l range in a number of 
         * parallel threads. The values are then stored in output/Fisher/.
         * This version uses a fixed k_stepsize of 0.0178.
         * @param lmin is the minimum of the l range.
         * @param lmax is the maximum of the l range.
         * @param n_points_per_thread is the number of points each thread is
         *          determining. It is used to obtain the total number of 
         *          lmodes calculated between the range and thus the l
         *          stepsize.
         * @param n_threads is the number of threads that should be used.
         *
         * TODO:
         * @return Returns 0 atm this is superfluous 
         */
        double F_fixed_kstepsize(int lmin, int lmax, int n_points_per_thread,\
                int n_threads);
        
        /**
         *
         */
        void Fl_varying_ksteps(int l, string param_key1,\
                string param_key2, int min_ksteps_Cl, int max_ksteps_Cl,\
                int ksteps_spacing);

        /**
         *
         */
        void Fl_varying_ksteps_smart(int l, string param_key1,\
                string param_key2, int min_ksteps_Cl, int max_ksteps_Cl);
        
    private:

        // **** Functions **** //

        /**
         * This function updates the model used to perform the calculations
         * given a set of new parameters. It updates:
         *      - comoving distance q(z),
         *      - the power spectrum Pk(z),
         *      - and the G21 function Tb(z).
         * @param new_params contains the dictionary of new parameters that
         *          the model should be updated to.
         * @param Pk_index points to the index this new model should be 
         *          identified as for the matter power spectrum.
         * @param Tb_index points to the index this new model should be 
         *          identified as for the 21cm signal.
         * @param Pk_index points to the index this new model should be 
         *          identified as for the comoving distance.

         */
        void update_Model(map<string, double> new_params, int *Pk_index,\
                int *Tb_index, int *q_index);

        /**
         * This function computes a covariance matrix containing all the 
         * Cl(k1,k2)'s for a given model, set of (k1,k2)'s and l. 
         * If this exact covariance matrix has already been calculated
         * before, the function finds the file where it was stored and just
         * reads it from there. Otherwise it calculates it and stores it. 
         *
         * @param l is the angular mode for which the Covariance matrix is 
         *          calculated.
         * @param Pk_index is the matter power spectrum model index used.
         * @param Tb_index is the 21cm signal model index used.
         * @param q_index is the comoving distance model index used.
         * @param krange is a list of k's that determine the modes which should
         *          be included in the covariance matrix.
         * @return Returns an aramdillo mat object, which is the covariance 
         *          matrix considering a certain model as a function of the
         *          angular mode l.
         */
        mat compute_Cl(int l, int Pk_index, int Tb_index, int q_index,\
                vector<double> krange);
        
        /**
         * 
         */
        double Cl_derivative(int l, string param_key, double k1, double k2,\
                int *Pk_index, int *Tb_index, int *q_index);
        /**
         *
         */
        double Cl_loglog_derivative(int l, string param_key, double k1,\
                double k2, int *Pk_index, int *Tb_index, int *q_index);

        /**
         * This function computes the derivative of the covariance matrix e 
         * with respect to some parameter. A four point derivative is used.
         * If this exact derivative matrix has already been calculated
         * before, the function finds the file where it was stored and just
         * reads it from there. Otherwise it calculates it and stores it. 
         *
         * @param l is the angular mode for which the Covariance matrix is 
         *          calculated.
         * @param param_key is the parameter key identifying the parameter 
         *          respect to which the derivative is to be taken.
         * @param Pk_index points to the matter power spectrum model index.
         * @param Tb_index points to the 21cm signal model index.
         * @param q_index points to the comoving distance model index.
         * @param krange is a list of k's that determine the modes which should
         *          be included in the covariance matrix.
         * @return Returns an aramdillo mat object, which is the derivative 
         *          matrix considering a certain model as a function of the
         *          angular mode l.
         */
        mat Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
                int *Tb_index, int *q_index, vector<double> krange);
        
        /**
         *
         */
        double compute_Fl(int l, string param_key1, string param_key2,\ 
                int ksteps_Cl, double *cond_num, int *Pk_index, int *Tb_index,\
                int *q_index);

        /**
         *
         */
        double compute_Fl(int l, string param_key1, string param_key2,\
                double kstepsize, double *cond_num, int *Pk_index,\
                int *Tb_index, int *q_index);

        /**
         * This function initializes the necessary interpolator objects needed
         * for the run. This is necessary so that the interpolators are ready 
         * when needed from each thread. Otherwise, each thread would start 
         * calculating the interpolators itself, which is undesirable and
         * can easily lead to bugs.
         *
         * @param param_key is the key of the parameter that is to be varied
         *          in the run.
         * @param Pk_index points to the matter power spectrum model index.
         * @param Tb_index points to the 21cm signal model index.
         * @param q_index points to the comoving distance model index.
         */
        void initializer(string param_key, int *Pk_index, int *Tb_index,\
                int *q_index);

        /**
         * This function writes a matrix to a binary file.
         *
         * @param matrix is the matrix that is to be written to file.
         * @param filename is the name of the file to which the matrix is 
         *          to be written.
         */
        void write_matrix(mat matrix, string filename);
        
        /**
         * This function reads a matrix from a binary file.
         *
         * @param filename is the name of the file from which the matrix is
         *          read.
         * @param n_rows is the number of rows in the matrix.
         * @param n_cols is the number of columns in the matrix.
         *
         * @return Returns an armadillo matrix object of the matrix read.
         */
        mat read_matrix(string filename, int n_rows, int n_cols);
        
        /**
         * This function checks whether a file with a given name already
         * exists or not.
         *
         * @param filename is the name of the file to be checked.
         *
         * @return Returns true if the file exists and false if not.
         */
        bool check_file(string filename);
        
        /**
         * This function records the run parameters to a file called 
         * RUN_INFO.dat so that one can later identify which parameters
         * were used in each run.
         *
         * @param lmin is the minimum l mode considered in the run.
         * @param lmax is the maximal l mode considered in the run.
         * @param lstepsize is the stepsize between adjacent l modes.
         * @param kstepsize is the stepsize between adjacent k modes.
         *
         * @return Returns filename prefix to be used to store resulting 
         *          Fl values for this run.
         */
        string update_runinfo(int lmin, int lmax,\
                int lstepsize, double kstepsize);

        /**
         *
         */
        vector<double> give_kmodes(int l, double kmax, int steps);
        
        /**
         *
         */
        vector<double> give_kmodes(int l, double k_max, double kstepsize);

        // **** Variables **** //

        vector<string> model_params_keys;
        CosmoCalc *CALC;
        ofstream Fl_file;
        map<string, double> current_params, fiducial_params, var_params;
        double kmin, kmax;
        vector<double> abcisses_done, logderivs_calculated,\
            abcisses_done_simple, derivs_calculated;
        bool noise, rsd;

};
