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

double fl (double x)
{
    return 0.551742 * x + 9778.15;
}

int main(int argc, char* argv[])
{ 
    map<string,double> params;
    ofstream outfile;
    //outfile.open("run_history.dat", ios::out | ios::app);

    //CosmoWrite writer(params);
    //writer.calculate_Cl_simple(2, 0.1, 0.00001, 1, 0.00001);
    //writer.generate_movie_Cl(1, 100, 0.5, 0.001, 1, 0.0001);
    //writer.generate_movie(30);
    //writer.calculate_integrandsimple(30, 1, 1, 100);
    //writer.calculate_integrandlong(30, 1, 1, 100);
    
    //writer.calculate_qdot();
    //writer.calculate_q();
    //writer.calculate_dTb(5, 20, 100);
    /* writer.calculate_integrandMM(198, 0.03, 0.03, 1000000);
   
    writer.calculate_integrandMN(142, 1, 1, 1000000);
    
    writer.calculate_integrandMN(142, 0.3, 0.3, 1000000);
    writer.calculate_integrandMN(142, 0.3, 1, 1000000);
    writer.calculate_integrandMN(142, 0.03, 1, 1000000);
    writer.calculate_integrandMN(142, 1.32, 2.5, 1000000);
    writer.calculate_integrandMN(142, 0.012, 0.012, 1000000);
    writer.calculate_integrandMN(142, 0.003, 0.003, 1000000);
    writer.calculate_integrandMN(142, 0.015, 0.015, 1000000);

    writer.calculate_integrandMM(199, 0.5, 0.5, 1000000);
    */ 


         
    clock_t t1, t2;
    string Fl_filepath = "output/Fls_rsd_simple_noq.dat"; 
    Fisher fish(params, Fl_filepath);
    
    t1 = clock();
    double res = fish.F("ombh2", "ombh2");
    t2 = clock();
    
    //outfile << " ##################### " << endl;
    //outfile << "kstep_Cl = 3" << endl;
    //outfile << "lmax = 15" << endl; 
    cout << "Result is " << res << endl;
    //outfile << "Result is " << res << endl;
    float diff ((float)t2 - (float)t1);
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    //outfile << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    
    /*     
    ifstream filesimlpe, filelong;
    filesimlpe.open("output/Fls_k4_simple.dat");
    filelong.open("output/Fls_k4_long.dat");
    vector<double> flsimple, flslong;
    double a, b;
    while (filesimlpe >> a >> b)
    {
        flsimple.push_back(b);
    }
    while (filelong >> a >> b)
    {
        flslong.push_back(b);
    }
    ofstream errorfile("output/error_fls_k4.dat");
    double err;

    for (int i = 0; i < flsimple.size(); ++i)
    {
        err = abs((flsimple[i] - flslong[i]) / flsimple[i]);
        errorfile << i << " " << err << endl;
    }

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
