#pragma once

#include "CosmologyCalculatorClass.hpp"
#include <iostream>
#include <fstream>

using namespace std;

class CosmoWrite : public CosmoCalc {
    public:
 
        // ------------ Functions -------------- //
        CosmoWrite(map<string, double> params);

        ~CosmoWrite();

        void calculate_Ml(int l, double k_fixed, double k2_low, double k2_high,\
                          int step, double stepsize = 0);
        void calculate_Nl(int l, double k_fixed, double k2_low, double k2_high,\
                          int step, double stepsize = 0);
        void calculate_distances(double zmax);
        void calculate_densities_rho();
        void calculate_densities_Omega();
        void calculate_H();
        void calculate_P();
        void calculate_P_CLASS();
        void calculate_dTb();
        void calculate_xHI();
        void calculate_Ts();
        void calculate_Tk();

};
