#pragma once
// This class will try to recover the results for Santos & Cooray 2006
#include "CosmoBasis.hpp"
#include <string>
#include <map>
#include <cmath>
#include <vector>
#include "stdafx.h"
#include "interpolation.h"
#include "Helper.hpp"

using namespace std;
using namespace alglib;

class Tomography : public CosmoBasis {
    public:
        Tomography(map<string, double> params);
        ~Tomography();
        double Cl(int l, double f1, double f2);
        double Cl_s(int l, double f1, double f2);
        double Cl_n(int l, double f1, double f2);
        double Cl_f(int l, double f1, double f2);
    private:
        // ------------ Functions -------------- //
        double T21(double f);
        double F(double k, double f);
        double I(int l, double k, double f);
        double J(int l, double k, double f);
        double f(double z);
        double W(double f0, double f);
        double redshift_z(double f);
        double T(double z);
        double Tk(double z);
        double ytot(double z);
        double Tc(double z);
        double yc(double z);
        double kappa10(double Tk);
        double nH(double z);
        double yalpha(double z);
        double xHI(double z);
        double b_b(double k, double z);
        double b_alpha(double k, double z);
        double b_xH(double k, double z);
        double g(double z);
        double P_dcdc(double k, double f1, double f2);
        double r(double z);
        double Cl_f(int i, int l, double f1, double f2);
        double _g(double z); 
        double Pkz_calc(double k, double z);
        double P_growth(double z);
        double D1(double z);
        double P_delta(double k, string units_k = "default",\
                string units_P = "default");
        double transfer(double x);
        
        // ------------ Variables -------------- //
        double f21 = 1420; // 21cm resframe frequency
        spline1dinterpolant kappa10_interpolator;
        vector<double> Ai, betai, alphai, xii;
};
