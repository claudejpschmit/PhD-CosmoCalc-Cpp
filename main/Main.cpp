#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include "CosmoBasis.hpp"
#include <fstream>
#include <sys/stat.h>
#include "Integrator.hpp"
#include "CosmologyCalculatorClass.hpp"
#include "ClassEngine.hpp"
#include "CosmologyWriterClass.hpp"
#include <time.h>
#include "FisherClass.hpp"
#include <armadillo>
#include "stdafx.h"
#include "interpolation.h"
#include "CAMB_interface.hpp"

using namespace std;
using namespace arma;
using namespace alglib;

double f (double x)
{
    return x*x;
}

int main(int argc, char* argv[])
{
     
    map<string,double> params;
    //CosmoWrite writer(params);
    //writer.calculate_P_compare(0.0001, 10, 10000, 7, 9, 3);

   /* 
    params.insert(pair<string,double>("ombh2",0.0226));
    params.insert(pair<string,double>("omch2",0.112));
    params.insert(pair<string,double>("omnuh2",0.00064));
    params.insert(pair<string,double>("omk",0.0));
    params.insert(pair<string,double>("hubble",70.0));

    params.insert(pair<string,double>("zmin",7.0));
    params.insert(pair<string,double>("zmax",9.0));
    params.insert(pair<string,double>("Pk_steps",3));
    CAMB_CALLER *camb;
    camb = new CAMB_CALLER;
    camb->call(params);
    vector<double> k, vP;
    vector<vector<double>> Pz;
    Pz = camb->get_Pz_values();
    k = camb->get_k_values();
    for (int i = 0; i < Pz.size(); ++i) {
        vP.insert(vP.end(), Pz[i].begin(), Pz[i].end());
    }
    cout <<vP.size()<<endl;
    cout << Pz.size() << endl;
    cout << k[1] << endl;
    delete camb;
    */
    //CosmoCalc calc(params);

    //cout << "test" << endl;
    
    /*** Bessel interpolation timer ***/
    /*
    CosmoWrite writer(params);
    
    clock_t t1, t2;
    t1 = clock();
    writer.calculate_bessels(5);    
    t2 = clock();
    float diff = (float)t2 - (float)t1;
    cout << "runtime for fancy interp was " << diff/CLOCKS_PER_SEC << endl;
    
    t1 = clock();
    writer.calculate_bessels_basic(5);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime for basic interp was " << diff/CLOCKS_PER_SEC << endl;
    
    t1 = clock();
    writer.calculate_bessels_cubic(5);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime for cubic interp was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_bessels_exact(5);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime for exact interp was " << diff/CLOCKS_PER_SEC << endl;
    */
    //writer.update_Pk_interpolator(writer.give_current_params());
    //writer.calculate_P_compare(0.0001, 10, 10000, 7, 9, 3);
    
    /* 
    clock_t t1, t2;

    t1 = clock();
    Fisher fish(params);
    cout << "Result is " << fish.F("ombh2", "ombh2")<< endl;
    t2 = clock();
    float diff ((float)t2 - (float)t1);
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    
    //fish.compute_Cl(10);
    //fish.show_Cl_mat();
    //fish.compute_Cl_inv();
    //fish.show_Cl_inv_mat();
    //double res = fish.Cl_loglog_derivative(142, "ombh2", 0.1, 0.1);
    //cout << res << endl;
    */
    (void) argc;
    (void) argv;
/*
    mat A = randu<mat>(2,2);
    mat B = randu<mat>(2,2);

    A(0,0) = 1;
    A(0,1) = 0.1;
    A(1,0) = 0;
    A(1,1) = 1.1;
    cout << A << endl;
    cout << B << endl;
    cout << A*B  << endl;
  */
    /*
       const char output_path[] = "output";
       if (stat(output_path, &sb) == 0 && S_ISDIR(sb.st_mode)) {
       cout << "output directory already exists." << endl;
       } else {
       mkdir(output_path, 0700);
       cout << "output directory generated!" << endl; 
       }

       auto g = [](double x) 
       {
       return 2*f(x);
       };
       cout << "integration yields "<< integrate(g, 0.0, 1.0, 100, simpson()) << endl;
       cout << "New integration yields " << integrate_simps(g, 0.0, 1.0, 100) << endl;
       map<string,double> params;

       params["O"] = 0.02;
       params["ombh2"] = 0.02;
       cout << params["O"] << endl;
       CosmoBasis base(params);
    //base.show_params(); 
    //cout << 1.34 * pow(10,2) << endl;

    ofstream output;
    const char filename[] = "bessels.dat";
    char path[50] = "";
    strcat(path, output_path);
    strcat(path, "/");
    strcat(path, filename);
    cout << path << endl;
    output.open(path);
    for (int n = 0; n < 100000; ++n) {
    output << n << " " << base.sph_bessel(10,double(n)/100.0) << endl;
    }
    output.close();

    //CosmoCalc calc(params);
    //calc.show_cosmo_calcs();
    vector<double> v;
    v.clear();
    cout << v.size()<< endl;

    //CosmoWrite writer(params);
    */
    /* 
       clock_t t1, t2;

       Fisher fish(params);
       t1 = clock();
       double res = fish.Cl_loglog_derivative(142, "ombh2", 0.01, 0.01);
       cout << res << endl;

    //fish.write_logder("ombh2", 0.0226, 0.0001, 0.01, 99, 142, 0.01, 0.01,\
    "_l142_0-01_0-01");
    t2 = clock();
    float diff ((float)t2 - (float)t1);
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    */
    //t1 = clock();
    //writer.calculate_Ml(5, 0.1, 0.01, 0.5, 10000); 
    //writer.calculate_distances(10);
    //t2 = clock();
    //float diff ((float)t2 - (float)t1);
    //cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    /* 
       t1 = clock();
    //writer.calculate_densities_rho(5000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    //writer.calculate_densities_Omega(10000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    //writer.calculate_H(1000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_P(0.001, 10, 10000, "default", "default");
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_P_CLASS(0.001, 10, 0, 10000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_dTb(5, 20, 100);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_xHI(0, 20, 100);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_Ts(0, 20, 100);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_Tk(0, 1000, 1000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

*/


    return 0;
}
