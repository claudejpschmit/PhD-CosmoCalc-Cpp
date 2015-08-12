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

    //CosmoWrite writer(params, &Pk_index, &Tb_index, &q_index);
    /*SanityChecker check(params, &Pk_index, &Tb_index, &q_index);
    //check.plot_integrand_zp(1200, 8.0, 1.0, 10000, "integrand_zp.dat");
    
    //check.plot_integrand_z(2000, 0.8, 0.75, 1000, "integrand_z.dat");

    //int l = 100;
    
    double k1 = 1.0;
    double k2 = 0.4;
    double k_low = 0.0001;
    double k_high = 2.0;
    int n = 8;
    double res;
    int counter = 0;
    ofstream file;
    stringstream name;
    name << "output/Cl_compare_gauss_Levin_highzsteps.dat";
    double ratio, time_r;
    file.open(name.str());
    for (int i = 0; i < 50; i++)
    {
        int l = 1 + i*100;
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
    
    Fisher fish(params,"Fisher.dat");
    cout << "The result is = " << fish.F("ombh2","ombh2")<< endl;
    /*
    double k1 = 0.5;
    double k2 = 0.7;
    ofstream file;
    file.open("test.dat");
    int l = 9;
    double res = fish.Cl_loglog_derivative(l, "ombh2", k1, k2, &Pk_index, &Tb_index, &q_index);   
    file << l << " " << res << endl;
    cout << l << " " << res << endl;

    #pragma omp parallel num_threads(6) private(Pk_index, Tb_index, q_index) 
    {
        #pragma omp for
        for (int i = 10; i < 50; i++)
        {
            l = i;
            Tb_index = 0;
            Pk_index = 0;
            q_index = 0;
            double res = fish.Cl_loglog_derivative(l, "ombh2", k1, k2, &Pk_index, &Tb_index, &q_index);   
            file << l << " " << res << endl;
            cout << l << " " << res << endl;
        }
    }
    file.close();
    */
    /*
    clock_t t1,t2;
    float diff1,diff2,diff3;
    SanityChecker check(params, &Pk_index, &Tb_index, &q_index);
    check.plot_intjj(2000, 8.0, 500, "test.dat");
    //check.plot_integrand_z(2000, 0.3, 0.5, 500, "integrand_z_new.dat");
    double k1 = 0.6;
    double k2 = 0.6;
    for (int l = 4000; l < 4025; l++){
        double cl1 = check.corr_Tb_rsd(l, k1, k2, 0.0001, 2, 0, 0, 0);
        t1 = clock();
        double cl2 = check.corr_Tb_MN(l, k1, k2, 0.0001, 2, 0, 0, 0);
        t2 = clock();
        diff1 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
        cout << l << " Rsd: " << cl1 << " " << cl2 << endl; 
        t1 = clock();
        double cl3 = check.corr_Tb(l, k1, k2, 0.0001, 2, 0, 0, 0);
        t2 = clock();
        diff2 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
        t1 = clock();
        double cl4 = check.corr_Tb_MM(l, k1, k2, 0.0001, 2, 0, 0, 0);
        t2 = clock();
        diff3 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;

        cout << "no Rsd: " << cl3 << " " << cl4 << " ratio = " << cl1/cl3  << endl;
        cout << "Time: " << diff1 << " " << diff2 << " " << diff3 << " " << diff1/diff2 << endl;
    }
    */
    //check.plot_integrad_z(221, 0.4, 0.4, 1000, "test.dat");
    //for (int l = 10; l < 1000; l++) 


    //CosmoWrite write(params, &Pk_index, &Tb_index, &q_index);
    //write.calculate_bessels_exact(1200);

    return 0;
}


