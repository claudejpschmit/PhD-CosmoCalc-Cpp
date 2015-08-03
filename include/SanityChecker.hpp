#pragma once

//This class should include all the functions to test the CosmoCalc class.
//ie. ultimately the best function will be determined and this class will 
//be superfluous... ha... like that is ever going to happen.
// This is mainly so I don't lose the overview of CosmoCalc.

#include "CosmologyCalculatorClass.hpp"

class SanityChecker : public CosmoCalc {
    
    public:
        SanityChecker(map<string, double> params);
        ~SanityChecker();
    
        void compare_interp(int l, double k1, double k2, double z, double *out1, double *out2);
        void compare_new(int l, double k1, double k2, double kappa, double *out1, double *out2);

        void kappa_integrand(int l, double k1, double k2, double z, double kappa, double zp, double *out1);
        void kappa_integral(int l, double k1, double k2, double z, double zp, double *out1, double k_low, double k_high);
        double limber3(int l, double z);
        void Cl_compare(int l, double k1, double k2);
        void zp_integrand(int l, double k1, double k2, double z, double zp, double *out1, double k_low, double k_high);
        void zp_integrand2(int l, double k1, double k2, double z, double zp, double *out1, double k_low, double k_high);
        void z_integrand_with_jq(int l, double k1, double kappa, double z, double *res);
        void z_integrand_no_jq(int l, double k1, double kappa, double z, double *res);
        double f(int l, double kappa, double z);
        double integral_z_jq(int l, double k1, double kappa);
        double integral_z_nojq(int l, double k1, double kappa);
        double integral_limber(int l, double k1, double kappa);


};
