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
    /*
    clock_t t1, t2;
    map<string,double> params;
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 

    SanityChecker check(params, &Pk_index, &Tb_index, &q_index);
    int l = 100;
    double k1 = 0.3;
    double k2 = 0.4;
    double k_low = 0.0001;
    double k_high = 2.0;
    t1 = clock();
    double res1 = check.corr_Tb(l, k1, k2, k_low, k_high, Pk_index, Tb_index, q_index);
    t2 = clock();
    float d1 = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
    t1 = clock();
    double res2 = check.Cl_MC(l, k1, k2, k_low, k_high,Pk_index, Tb_index, q_index);
    t2 = clock();
    float d2 = ((float)t2-(float)t1)/CLOCKS_PER_SEC;
    cout << res1 << " " << res2 << endl;
    cout << "time: " << d1 << " " << d2 << endl;
    */
    
    
    map<string,double> params;    
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
    
    //Fisher fish(params, "aaa.dat");
    //cout << fish.F("ombh2", "ombh2") << endl;
    params.insert(pair<string,double>("kmax",1));
    params.insert(pair<string,double>("zmax",8));
    params.insert(pair<string,double>("zsteps",500));

    //Fisher fish(params,"Fisher_singlel_kmax_1_new.dat");
    Fisher fish(params, "bla.dat");
    int l = 1000;
    int min_ksteps_Cl = 4;
    int max_ksteps_Cl = 6;
    int ksteps_spacing = 5;
    fish.Fl_varying_ksteps(l, "ombh2", "ombh2", min_ksteps_Cl, max_ksteps_Cl, ksteps_spacing);
    mat A = fish.read_matrix("output/matrices/Cl_1000_1_5_7_8_nrnn.bin",5,5);  
    mat B = fish.read_matrix("output/matrices/Cl_1000_1_5_7_8_nrn.bin",5,5); 
    cout << A << endl;
    cout << "-------------------------" << endl;
    cout << B << endl;
    return 0;
}


