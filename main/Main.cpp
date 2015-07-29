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
/*
template<typename T>
double integrate_levin2(T &f, const int nterm = 10)
{
    const double pi = boost::math::constants::pi<double>();
    double beta = 1.0, a = 0.0, b = 0.0, sum = 0.0;
    double ans;
    if (nterm > 100)
    {
        cout << "nterm too large" << endl;
        throw("nterm too large");
    }
    else {
        Levin series(100,0.0);
        cout << setw(5) << "N" << setw(19) << "Sum (direct)" << setw(21) << "Sum (Levin)" << endl;
        cout << "what?" << endl;
        for (int n = 0; n<=nterm;n++) {
            b+=pi;
            cout << " qromb " << endl;
            double s = qromb(f, a, b, 1.0E-8);
            cout << " qromb done " << endl;
            a=b;
            sum += s;
            double omega = (beta+n)*s;
            ans = series.next(sum, omega, beta);
            cout << setw(5) << n << fixed << setprecision(14) << setw(21) << sum << setw(21) << ans << endl;
        }
    }
    return ans;
}
*/

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
    //fun();
    //cout << integrate_simps(f,1,2,100000) << endl;
    //cout << qromb(f,1,2,10E-10) << endl; 
    /*
    CosmoCalc cosmo(params);
    cout << boost::math::sph_bessel(100,1000) << " " << cosmo.sph_bessel_camb(100,1000) << endl;
    auto integrand1 = [&](double x)
    {
        if (x == 0.0)
            return 0.0;
        else
            return cosmo.sph_bessel_camb(10, x) * cosmo.sph_bessel_camb(10, x);
    };
    //cout << qromb(f, 0, 1000000, 1.0E-8) << endl;
    cout << qromb(integrand1, 0, 1000000, 1.0E-8) << endl;

    //cout << integrate_levin(f,0,100000) << endl;
    //cout << integrate_levin(integrand1,0,100000)<< endl;
    //cout << cosmo.limber(100,1) << endl;
    double res = cosmo.corr_Tb(100,0.4,0.5, 0.0001, 1);
    cout << res << endl;
    //double res1 = cosmo.corr_Tb_new(100,0.4,0.5, 0.0001, 1);
    double res2 = cosmo.corr_Tb_new2(100,0.4,0.5, 0.0001, 1);
    cout << res2 << endl;
    double res_s = cosmo.Cl_simplified(100, 0.4, 0.5);
    double res_sl = cosmo.Cl_simplified_levin(100, 0.4, 0.5);

    cout << res2 << " " << res_s << " " << res_sl << " "  << endl;
    */
    /*
    ofstream fout;
    fout.open("output/interp.dat");
    double z = 6.5;
    while (z <=9.5){
        fout << z << " " << cosmo.delta_Tb_bar_G21(z) << " " << cosmo.delta_Tb_bar(z)/1000.0 <<\
            " " << cosmo.Tb_interp(z)/1000.0 << endl;
        z += 0.01;
    }
    fout.close();
    */
    //ofstream outfile;
    //outfile.open("run_history.dat", ios::out | ios::app);
    
    clock_t t1, t2;

    CosmoWrite writer(params);   
    writer.calculate_inverse_comoving(15000);
    //SanityChecker check(params);
    //ofstream file;
    //file.open("zp_integrand3.dat");
    //double res1, res2;
    //check.compare_interp(1000, 0.6, 0.5, 8.0,  &res1, &res2);
    //check.compare_interp(1000, 0.4, 0.5, 8.3,  &res1, &res2);
    //check.compare_interp(1000, 0.4, 0.5, 8.6341,  &res1, &res2);
    //check.compare_interp(1000, 0.4, 0.5, 7.0,  &res1, &res2);
    //check.compare_interp(1000, 0.4, 0.5, 7.1231,  &res1, &res2);
    //check.compare_interp(1000, 0.4, 0.5, 7.5611,  &res1, &res2);
    //check.compare_interp(1000, 0.2, 0.3, 8.0,  &res1, &res2);
    //check.compare_interp(1000, 0.2, 0.3, 7.5,  &res1, &res2);
    //check.Cl_compare(1000, 0.4, 0.5);
    /*for (int i = 0; i < 2000; i++) {
        //double k1 = 0+i * 0.0001;
        double zp = 7.8 + i * 0.0002;
        //cout << " for kappa = " << kappa << endl;
        check.zp_integrand2(1000, 0.6, 0.5, 8.0, zp, &res1, 0.1, 1.0);
        //check.kappa_integrand(1000, 0.4, 0.2, 8.0, kappa, 8.2, &res1);
        //check.compare_interp(1000, k1, 0.5, 8.0, &res1, &res2);
        //file << k1 << " " << res1 << " " << res2 << endl;
        file << zp << " " << res1 << endl;
    }*/
    /* 
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 0.3);
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 0.6);
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 1.0);
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 1.5);
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 2);
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 4);
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 20);
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 100);
    check.kappa_integral(1000, 0.6, 0.5, 7.5, 8.5, &res2, 0.1, 200);
    */
   // file.close();
    
    /*
    writer.compare(1000,0.4,0.5);

    for (int i = 0; i < 100; i++) {
        double k1 = 0+i * 0.01;
        cout << " for k1 = " << k1 << endl;
        writer.compare(1000, k1,0.5);
    }
    */
    //writer.calculate_Cl_simple(1000, 0.5, 0.08, 1, 0.0001);

    //cout << "simple done" << endl;
    //params["zsteps"] = 10;
    params.insert(pair<string,double>("k_stepsize",0.0001));
    //CosmoWrite writer2(params);
    //writer.calculate_Cl_full(1000, 0.5, 0.1, 1, 0.0001);
    //writer.calculate_bessels_exact(1000);
    //writer.calculate_bessels_cubic(1000);
    
    //cosmo.compare(1000, 0.5, 0.5);
    //cout << cosmo.limber2(500,24) << endl;
    //cout << cosmo.corr_Tb_new(100, 0.5, 0.5, 0.4, 0.6) << endl;
    //cout << cosmo.Cl_simplified(100, 0.5, 0.5)<< endl;
    /*  ofstream file;
    file.open("output/limber.dat");
    double q = 0.99;
    int l = 3000;
    for (int i = 0; i < 2000; i++) {
        q += 0.00001;
        file << q << " " << 2*(l+0.5)*(l+0.5)*cosmo.limber(l, 1,q)/3.14159265359<< endl;
    }
    file.close();*/
    //writer.calculate_Ml(1,2000,0.01,0.2);
    //writer.calculate_Nl(1,200,0.01,0.02);
    //writer.calculate_Cl_full(100, 0.5, 0.08, 1, 0.0001);
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

/*
    clock_t t1, t2;
    string Fl_filepath = "output/G21_included/Fls_Noise.dat"; 
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
*/
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
