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
#include <gsl/gsl_integration.h>

using namespace std;
using namespace arma;
using namespace alglib;

int main(int argc, char* argv[])
{

    map<string,double> params;
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0; 
    
    SanityChecker check(params, &Pk_index, &Tb_index, &q_index);
    int l = 100;
    double k1 = 0.3;
    double kappa = 0.5;
    double res1 = check.M(l, k1, kappa, Pk_index, Tb_index, q_index);
    double res2 = check.M_gsl(l, k1, kappa, Pk_index, Tb_index, q_index);
    cout << res1 << " " << res2 << endl;
    /*
    map<string,double> params;    
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
  
    //Fisher fish(params, "Fisher_norsd.dat");
    //cout << fish.F("ombh2", "ombh2") << endl;
    
    params.insert(pair<string,double>("kmax",1));
    params.insert(pair<string,double>("zmax",8));
    params.insert(pair<string,double>("zsteps",500));

    Fisher fish(params,"Fisher_singlel_kmax_1_wNoise.dat");
    int l = 1000;
    int min_ksteps_Cl = 4;
    int max_ksteps_Cl = 100;
    int ksteps_spacing = 5;
    fish.Fl_varying_ksteps(l, "ombh2", "ombh2", min_ksteps_Cl, max_ksteps_Cl, ksteps_spacing);
    */
    return 0;
}


