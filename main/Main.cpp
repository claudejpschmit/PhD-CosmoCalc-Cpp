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
#include "Global21cmInterface.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include "SanityChecker.hpp"
#include "LevinIntegrator.hpp"

using namespace std;
using namespace arma;
using namespace alglib;



double fl (double x)
{
    return boost::math::cyl_bessel_j(0,100*x);
}
double f(double x)
{
    return x;
}

double integr(double r1, double r2, int l, double a, double b, int nsteps)
{
    auto integrand = [&](double x)
    {
        return x*boost::math::sph_bessel(l, r1*x)*boost::math::sph_bessel(l, r2*x);
    };

    double res = integrate_simps(integrand, a, b, nsteps);
    return res;
}

int main(int argc, char* argv[])
{
    clock_t t1, t2;
    double r1 = 1000;
    double r2 = 1000;
    double nu = 12;
    int l = 1;
    int nsteps = 15;
    double a = 1;
    double b = 2;
    /*
    Levin levin(a,b);
    double res1 = levin.integrate_2sphj(&f, r1,r2, l, nsteps);
    double res2 = integr(r1,r2, l, a, b, 100000);
    cout << res1 << endl;
    cout << res2<< endl;
    */
    
    ofstream file;
    file.open("Levin2");

    double res2 = integr(r1,r2,nu, a, b, 100000);
    for (int i = 4; i < 100; i++)
    {
        nsteps = i;
        Levin levin(a,b);
        double res1 = levin.integrate_2sphj(&f, r1, r2, nu, nsteps);

        cout << res1 << endl;
        file << nsteps << " " << res1 << " " << res2 << " " << res1/res2 << endl;
    }
    file.close();
    
    /*
    t1 = clock();
    double res2 = integr(r1,r2, l, a, b, 300000);
    t2 = clock();
    float time2 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    cout << time2 << endl;
    cout << res2 << endl;
    */


    

    return 0;
}
