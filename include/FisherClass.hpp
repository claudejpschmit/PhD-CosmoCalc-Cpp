#pragma once

#include "CosmologyCalculatorClass.hpp"

class Fisher {
    public:
        Fisher(map<string, double> params);
        ~Fisher();
    
    private:
        CosmoCalc *CALC;
        map<string, double> current_params, fiducial_params, var_params;
        vector<double> Cl, Cl_inv, krange;
        double kmin, kmax;
};
