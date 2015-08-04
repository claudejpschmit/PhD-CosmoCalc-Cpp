#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include "CosmoBasis.hpp"
#include <fstream>
#include <sys/stat.h>
#include "Integrator.hpp"
#include "CosmologyCalculatorClass.hpp"
//#include "ClassEngine.hpp"
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

using namespace std;
using namespace arma;
using namespace alglib;

int main(int argc, char* argv[])
{
    map<string, double> params;
    SanityChecker check(params);

    int l = 100;
    double k1 = 0.3;
    double k2 = 0.4;
    double z = 7.8;
    double zp = 8.0;
    double k_low = 0.1;
    double k_high = 1.0;

    double res;
    check.kappa_integral(l, k1, k2, z, zp, &res, k_low, k_high);
    
    cout << res << endl;
    return 0;
}
