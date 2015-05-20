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

using namespace std;


double f (double x)
{
    return x*x;
}

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    struct stat sb;

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
    base.show_params(); 
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
    
    clock_t t1, t2;
    CosmoWrite writer(params);
    t1 = clock();
    //writer.calculate_Ml(5, 0.1, 0.01, 0.5, 10000); 
    writer.calculate_distances(10);
    t2 = clock();
    float diff ((float)t2 - (float)t1);
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    
    return 0;
}
