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

using namespace std;
using namespace arma;
using namespace alglib;

double f (double x)
{
    if (x == 0.0)
        return 0.0;
    else
        return boost::math::sph_bessel(1,100*x) * boost::math::sph_bessel(1,101*x);
}

double fl (double x)
{
    return boost::math::cyl_bessel_j(0,100*x);
}

double u(double x, int k, double d)
{
    return pow(x-d, k-1);
}

double up(double x, int k, double d)
{
    return (k-1)*pow(x-d, k-2);
}


void fun()
{
    mat matrix, rhs, c, point;
    int n = 9;
    double a = 1;
    double b = 2;
    double r = 100;
    int l = 0;
    matrix = randu<mat>(2*n,2*n);
    rhs = randu<mat>(2*n,1);
    point = randu<mat>(n,1);
    double d = (a+b)/2.0 + 0.00000000000001;
    for (int i = 0; i < n; i++)
    {
        point(i,0) = a + (i-1.0)*(b-a)/(double)(n-1.0);
    }

    for (int i = 0; i < n; i++)
    {
        double x = point(i,0);
        rhs(i,0) = 1.0; //f(x)
        rhs(i+n,0) = 0; //g(x)

        for (int k = 0; k < n; k++)
        {
            matrix(i, k) = up(x,k,d) + (double)l * u(x,k,d)/x;
            matrix(i, k + n) = r * u(x,k,d);
            matrix(i + n, k) = -r * u(x,k,d);
            matrix(i + n, k + n) = up(x,k,d) - (double)(l+1)*u(x,k,d)/x;
        }
    }
    matrix = matrix.i();
    c = randu<mat>(n*2,1);
    c = matrix * rhs;

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;

    for (int k = 0; k < n; k++)
    {
        sum1+=c(k,0) * u(b,k,d);
        sum2+=c(k,0) * u(a,k,d);
        sum3+=c(k+n,0) * u(b,k,d);
        sum4+=c(k+n,0) * u(a,k,d);
    }
    double j1 = boost::math::cyl_bessel_j(l,r*b);
    double j2 = boost::math::cyl_bessel_j(l,r*a);
    double j3 = boost::math::cyl_bessel_j(l+1,r*b);
    double j4 = boost::math::cyl_bessel_j(l+1,r*a);
    double result = sum1 * j1 - sum2 * j2 + sum3 * j3 - sum4 * j4;

    cout << result << endl;
}

void fun2()
{
    mat matrix, rhs, c, point;
    int n = 9;
    double a = 1;
    double b = 2;
    double r1 = 100;
    double r2 = 101;
    int l = 1;
    matrix = randu<mat>(4*n,4*n);
    rhs = randu<mat>(4*n,1);
    point = randu<mat>(n,1);
    double d = (a+b)/2.0 + 0.00000000000001;
    for (int i = 0; i < n; i++)
    {
        point(i,0) = a + (i-1.0)*(b-a)/(double)(n-1.0);
    }

    for (int i = 0; i < n; i++)
    {
        double x = point(i,0);
        rhs(i,0) = 1.0; //f(x)
        rhs(i+n,0) = 0; //g(x)
        rhs(i+2*n,0) = 0;
        rhs(i+3*n,0) = 0;

        for (int k = 0; k < n; k++)
        {
            matrix(i, k) = up(x,k,d) -2*(double)(l+1) * u(x,k,d)/x;
            matrix(i, k + n) = r1 * u(x,k,d);
            matrix(i, k + 2*n) = r2 * u(x,k,d);
            matrix(i, k + 3*n) = 0;

            matrix(i + n, k) = -r1 * u(x,k,d);
            matrix(i + n, k + n) = up(x,k,d) - 2*u(x,k,d)/x;
            matrix(i + n, k + 2*n) = 0;
            matrix(i + n, k + 3*n) = r2*u(x,k,d);

            matrix(i + 2*n, k) = -r2 * u(x,k,d);
            matrix(i + 2*n, k + n) = 0;
            matrix(i + 2*n, k + 2*n) = up(x,k,d)-2*u(x,k,d)/x;
            matrix(i + 2*n, k + 3*n) = r1*u(x,k,d);

            matrix(i + 3*n, k) = 0;
            matrix(i + 3*n, k + n) = -r2 * u(x,k,d);
            matrix(i + 3*n, k + 2*n) = -r1 * u(x,k,d);
            matrix(i + 3*n, k + 3*n) = up(x,k,d) + 2 *(double)(l-1)*u(x,k,d)/x;
        }
    }
    matrix = matrix.i();
    c = randu<mat>(4*n,1);
    c = matrix * rhs;

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;
    double sum7 = 0;
    double sum8 = 0;

    for (int k = 0; k < n; k++)
    {
        sum1+=c(k,0) * u(b,k,d);
        sum2+=c(k,0) * u(a,k,d);
        sum3+=c(k+n,0) * u(b,k,d);
        sum4+=c(k+n,0) * u(a,k,d);
        sum5+=c(k+2*n,0) * u(b,k,d);
        sum6+=c(k+2*n,0) * u(a,k,d);
        sum7+=c(k+3*n,0) * u(b,k,d);
        sum8+=c(k+3*n,0) * u(a,k,d);

    }
    double j1 = boost::math::sph_bessel(l,r1*b);
    double j2 = boost::math::sph_bessel(l,r1*a);
    double j3 = boost::math::sph_bessel(l-1,r1*b);
    double j4 = boost::math::sph_bessel(l-1,r1*a);

    double j5 = boost::math::sph_bessel(l,r2*b);
    double j6 = boost::math::sph_bessel(l,r2*a);
    double j7 = boost::math::sph_bessel(l-1,r2*b);
    double j8 = boost::math::sph_bessel(l-1,r2*a);


    double result = sum1 * j1 * j5 - sum2 * j2 * j6 +\
                    sum3 * j3 * j5 - sum4 * j4 * j6 +\
                    sum5 * j1 * j7 - sum6 * j2 * j8 +\
                    sum7 * j3 * j7 - sum8 * j4 * j8;

    cout << result << endl;
}
void fun3()
{
    mat matrix, rhs, c, point;
    int n = 8;
    double a = 1;
    double b = 2;
    double r = 10;
    int l = 1; // calculate integral for l-1
    matrix = randu<mat>(3*n,3*n);
    rhs = randu<mat>(3*n,1);
    point = randu<mat>(n,1);
    double d = (a+b)/2.0 + 0.00000000000001;
    for (int i = 0; i < n; i++)
    {
        point(i,0) = a + (i-1.0)*(b-a)/(double)(n-1.0);
    }

    for (int i = 0; i < n; i++)
    {
        double x = point(i,0);
        rhs(i,0) = 1.0; //f(x)
        rhs(i+n,0) = 0; //g(x)
        rhs(i+2*n,0) = 0;

        for (int k = 0; k < n; k++)
        {
            matrix(i, k) = up(x,k,d) + 2*(double)(l-1) * u(x,k,d)/x;
            matrix(i, k + n) = -2*r * u(x,k,d);
            matrix(i, k + 2*n) = 0;

            matrix(i + n, k) = r * u(x,k,d);
            matrix(i + n, k + n) = up(x,k,d) - u(x,k,d)/x;
            matrix(i + n, k + 2*n) = -r * u(x,k,d);

            matrix(i + 2*n, k) = 0;
            matrix(i + 2*n, k + n) = 2.0*r * u(x,k,d);
            matrix(i + 2*n, k + 2*n) = up(x,k,d)-2.0*(double)l*u(x,k,d)/x;
        }
    }
    matrix = matrix.i();
    c = randu<mat>(3*n,1);
    c = matrix * rhs;

    double sum1 = 0;
    double sum2 = 0;
    double sum3 = 0;
    double sum4 = 0;
    double sum5 = 0;
    double sum6 = 0;

    for (int k = 0; k < n; k++)
    {
        sum1+=c(k,0) * u(b,k,d);
        sum2+=c(k,0) * u(a,k,d);
        sum3+=c(k+n,0) * u(b,k,d);
        sum4+=c(k+n,0) * u(a,k,d);
        sum5+=c(k+2*n,0) * u(b,k,d);
        sum6+=c(k+2*n,0) * u(a,k,d);
    }
    double j1 = boost::math::cyl_bessel_j(l,r*b);
    double j2 = boost::math::cyl_bessel_j(l,r*a);
    double j3 = boost::math::cyl_bessel_j(l-1,r*b);
    double j4 = boost::math::cyl_bessel_j(l-1,r*a);

    double result = sum1 * j3 * j3 - sum2 * j4 * j4 +\
                    sum3 * j1 * j3 - sum4 * j2 * j4 +\
                    sum5 * j1 * j1 - sum6 * j2 * j2;

    cout << result << endl;
}


int main(int argc, char* argv[])
{ 
    map<string,double> params;
    SanityChecker check(params);
    ofstream file;

    double res1, res2;
    int l = 1000;
    double k1 = 0.11;
    int nmax = 10;
    double zmin = 7;
    double zmax = 9;
    double zstepsize = 0.001;
    int zsteps = (zmax - zmin)/zstepsize;
    double kappa = 0.1;
    double I1 = check.integral_z_jq(l, k1, kappa);
    double I2 = check.integral_z_nojq(l, k1, kappa);
    double limber = check.integral_limber(l,k1,kappa);

    file.open("test.dat");
    
    for (int i = 0; i< 1000; i++)
    {
        kappa += 0.001;
        double I1 = check.integral_z_jq(l, k1, kappa);
        double I2 = check.integral_z_nojq(l, k1, kappa);
        double limber = check.integral_limber(l,k1,kappa);

        file << kappa << " " << I1 << " " << I2 << " " << limber << endl;
    }
    
    file.close();
  
    /*
    for (int n = 0; n < nmax; n++)
    {
        double kappa = 0.1 + n * 0.1;
        stringstream filename;
        filename << "f_" << kappa << ".dat"; 
        file.open(filename.str());

        for (int i = 0; i < zsteps; i++)
        {
            double z = zmin + i*zstepsize;
            double f = check.f(l, kappa, z);
            double j = check.sph_bessel_camb(l, k1 * check.D_C(z));
        
            file << z << " " << f << " " << j << endl;
        }
        file.close();
    }
    */
    return 0;
}
