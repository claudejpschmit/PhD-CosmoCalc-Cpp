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

using namespace std;
using namespace arma;
using namespace alglib;

int main(int argc, char* argv[])
{
    map<string, double> params;
    
    // Uncomment for Cl_compare
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;

    //SanityChecker check(params, &Pk_index, &Tb_index, &q_index);
    //check.plot_integrad_z(2000, 0.8, 0.75, 1000, "integrand_z.dat");
/*
    int l = 100;
    double k1 = 0.3;
    double k2 = 0.6;
    double k_low = 0.01;
    double k_high = 2.0;
    int n = 8;
    double res;
    int counter = 0;
    ofstream file;
    stringstream name;
    name << "output/Cl_compare_" << n << "_2.dat";
    double ratio, time_r;
    file.open(name.str());
    for (int i = 0; i < 1; i++)
    {
        l = 1600 + i * 100;
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
        double res2 = check.kappa_integrand(1600, z, zp, kappa);
        file << kappa << " " << res2 << endl;
    } 
    file.close();
    
*/   
    Fisher fish(params,"Fisher3.dat");
    cout << "The result is = " << fish.F("ombh2","ombh2")<< endl;


    return 0;
}


