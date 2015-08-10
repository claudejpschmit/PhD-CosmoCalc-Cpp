#pragma once

//This class should include all the functions to test the CosmoCalc class.
//ie. ultimately the best function will be determined and this class will 
//be superfluous... ha... like that is ever going to happen.
// This is mainly so I don't lose the overview of CosmoCalc.

#include "CosmologyCalculatorClass.hpp"

class SanityChecker : public CosmoCalc {
    
    public:
        SanityChecker(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index);
        ~SanityChecker();
    
        void kappa_integral(int l, double z, double zp, double *out1, double k_low, double k_high, int n);
        double kappa_integrand(int l, double z, double zp, double kappa);
        void Compare_Cl(int l, double k1, double k2, double k_low,\
                double k_high, int n_levin, double *ratio, double *time_r);
        void plot_integrad_z(int l, double k1, double k2, int zsteps, string filename);
        void plot_intjj(int l, double zp, int zsteps, string filename);
};
