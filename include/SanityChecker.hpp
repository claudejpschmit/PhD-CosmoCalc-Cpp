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
    
        void kappa_integral(int l, double z, double zp, double *out1, double k_low, double k_high);
        double kappa_integrand(int l, double z, double zp, double kappa);
};
