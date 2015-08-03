#pragma once

#include <iostream>
#include <boost/math/constants/constants.hpp>

using namespace std;

class Levin 
{
    public:
        Levin(double a, double b);
        ~Levin();
        
        /** 
         *  this gives int_a^b f(x)J_nu(r*x) dx
         */
        //template<typename Function>
        double integrate_singleJ(double (*f)(double), double r, double nu, int n);
        double integrate_doubleJ(double (*f)(double), double r, double nu, int n);
        double integrate_2sphj(double (*f)(double), double r1, double r2, int l, int n);


    protected:
        double bessel_J(double nu, double x);
        double sph_bessel(int l, double x);
        double u(int k, double x);
        double up(int k, double x);

        double a, b, d;
        const double pi = boost::math::constants::pi<double>();

};
