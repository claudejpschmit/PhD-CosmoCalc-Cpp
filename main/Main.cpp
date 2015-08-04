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
    double k_low = 0.001;
    double k_high = 1.0;

    double res;
    
    for (int i = 0; i < 100; i++)
    {
        z = 7 + i*0.02;
        check.kappa_integral(l, z, zp, &res, k_low, k_high);
        if (res < 0.9 || res > 1.1) 
            cout << z << " " << res << endl;
    } 
    
    ofstream file;
    file.open("integrand.dat");
    z = 8.9;
    for (int i = 0; i < 10000; i++)
    {
        double kappa = 0.001 + i*0.0001;
        double res2 = check.kappa_integrand(l, z, zp, kappa);
        file << kappa << " " << res2 << endl;
    } 
    file.close();

    return 0;
}
