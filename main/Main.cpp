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
#include "stdafx.h"
#include "interpolation.h"

using namespace std;
using namespace alglib;

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

    CosmoCalc calc(params);
    calc.show_cosmo_calcs();
    
    // testing interpolator
    real_1d_array x;
    double _x[] = {0.0, 0.0, 0.0, 0, 0};
    x.setcontent(5,_x);
    x(0) = 0.0;
    x(1) = 0.25;
    x(2) = 0.5;
    x(3) = 0.75;
    x(4) = 1.0;
    real_1d_array y = "[0.0, 0.5, 1.0]";
    real_1d_array f = "[0.00, 0.0625, 0.25, 0.5625, 1.0, 0.5, 0.5625, 0.75, 1.0625, 1.5, 2.00, 2.0625, 2.25, 2.5625, 3.00]";
    double vx = 0.4214;
    double vy = 0.621;
    double v, dx, dy, dxy;
    spline2dinterpolant s;

    spline2dbuildbicubicv(x, 5, y, 3, f, 1, s);
    v = spline2dcalc(s, vx, vy);
    cout << v << endl;



    return 0;
}
