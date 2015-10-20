#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include "CosmoBasis.hpp"
#include <fstream>
#include <sys/stat.h>
#include "Integrator.hpp"
#include "CosmologyCalculatorClass.hpp"
#include "CosmologyWriterClass.hpp"
#include <time.h>
#include "FisherClass.hpp"
#include <armadillo>
#include "stdafx.h"
#include "interpolation.h"
#include "CAMB_interface.hpp"
#include "Global21cmInterface.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include "SanityChecker.hpp"
#include "LevinIntegrator.hpp"
//#include <gsl/gsl_integration.h>


using namespace std;
using namespace arma;
using namespace alglib;

int main(int argc, char* argv[])
{

    map<string,double> params;    
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
     
    params.insert(pair<string,double>("kmax",1));
    params.insert(pair<string,double>("zmax",8));
    params.insert(pair<string,double>("zsteps",500));
    params.insert(pair<string,double>("w_DE",-0.9));
    //params.insert(pair<string,double>("omnuh2", 0.0007));
    vector<string> keys = {"ombh2", "omch2", "hubble", "fesc", "fstar"};
    Fisher fish(params, "delete_me.dat", keys);
    
    //fish.F_fixed_kstepsize(1000,5000,5,7);
    //vector<double> krange;
    //krange.push_back(0.1);
    //krange.push_back(0.2);
    //krange.push_back(0.3);

    //mat B = fish.Cl_derivative_matrix(1000, "fesc", &Pk_index,\
        &Tb_index, &q_index, krange);
    //cout << B << endl;
    //cout << setprecision(10) << CALC.Tb_interp(10,Tb_index) << endl;
    //cout << CALC.Pk_interp(0.1,8,Pk_index) << endl; 
    //cout << fish.F_fixed_kstepsize("ombh2", "ombh2") << endl;

    //mat A = fish.read_matrix("output/matrices/Cla_fesc_2995_0.325404_1.0018_39_7_8_nrnn.bin",39,39);
    //cout << A << endl;
    //CosmoWrite writer(params, &Pk_index, &Tb_index, &q_index);
    //writer.calculate_Cl_full(100, 0.3, 0.001, 1, 0.001);
    
    //Fisher fish(params,"Fisher_singlel_2000_kmax_1.dat");
    //Fisher fish(params, "bla2.dat");
    //mat A  = fish.read_matrix("output/matrices/Cl_1000_0.108649_1_65_7_8_nrnn.bin",65,65);
    //mat B = pinv(A,0.001);
    //cout << A * B << endl;
    //mat B  = fish.read_matrix("output/matrices/Cl_1000_1_5_7_8_nrnn.bin",5,5);
    //cout << A << endl;
    //cout << "-------------" << endl;
    //cout << B << endl;
   
    //int l = 2000;
    //int min_ksteps_Cl = 34;
    //int max_ksteps_Cl = 60;
    //int ksteps_spacing = 5;
   
    //fish.Fl_varying_ksteps_smart(l, "ombh2", "ombh2", min_ksteps_Cl, max_ksteps_Cl);
    //fish.Fl_varying_ksteps(l, "ombh2", "ombh2", min_ksteps_Cl, max_ksteps_Cl, ksteps_spacing);
    return 0;
}


