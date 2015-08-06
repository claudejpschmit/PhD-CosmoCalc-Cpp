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
    
    // Uncomment for Cl_compare
    /*

    SanityChecker check(params);

        int l = 100;
    double k1 = 0.3;
    double k2 = 0.6;
    double k_low = 0.01;
    double k_high = 1.0;
    int n = 8;
    double res;
    int counter = 0;
    ofstream file;
    stringstream name;
    name << "output/Cl_compare_" << n << "_new2.dat";
    double ratio, time_r;
    file.open(name.str());
    for (int i = 0; i < 50; i++)
    {
        l = 100 + i * 100;
        cout << "----- calc for l = " << l << endl;
        check.Compare_Cl(l, k1, k2, k_low, k_high,n,&ratio,&time_r);
        file << l << " " << ratio << " " << time_r << endl;
    }
    file.close();
    cout << "This was for n = " << n << endl;
    file.open("integrand.dat");
    double z = 7.0;
    double zp = 8.0;
    for (int i = 0; i < 20000; i++)
    {
        double kappa = 0.001 + i*0.0001;
        double res2 = check.kappa_integrand(l, z, zp, kappa);
        file << kappa << " " << res2 << endl;
    } 
    file.close();
    */
   
    Fisher fish(params,"Fisher3.dat");
    fish.F("ombh2","ombh2");

    return 0;
}


