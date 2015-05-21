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

using namespace std;
using namespace arma;

double f (double x)
{
    return x*x;
}

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    struct stat sb;

    mat A = randu<mat>(2,2);
    mat B = randu<mat>(2,2);
    
    //matrix multiplication not working.
    //cout << A * B << endl;



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

    Fisher fish(params);
    //fish.write_logder("ombh2", 0.0226, 0.0001, 0.01, 99, 143, 0.1, 0.1,\
                      "_l143_0-1_0-1");

    /*
    clock_t t1, t2;
    CosmoWrite writer(params);
    t1 = clock();
    //writer.calculate_Ml(5, 0.1, 0.01, 0.5, 10000); 
    //writer.calculate_distances(10);
    t2 = clock();
    float diff ((float)t2 - (float)t1);
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    
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
