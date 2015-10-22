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
    params.insert(pair<string,double>("kmax",1));
    params.insert(pair<string,double>("zmax",8));
    params.insert(pair<string,double>("zsteps",500));
    vector<string> keys = {"ombh2", "omch2", "hubble", "fesc"};
    Fisher fish(params, "delete_me.dat", keys);
    
    fish.F_fixed_kstepsize(1000,5000,5,7);
    
    return 0;
}


