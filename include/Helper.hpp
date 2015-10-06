#pragma once

#include "stdafx.h"
#include "interpolation.h"
#include <vector>
#include <armadillo>
#include <string>
using namespace std;
using namespace alglib;
using namespace arma;

struct Pk_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, tcmb, w_DE;
    spline2dinterpolant interpolator;
};

struct Tb_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, s8, T_CMB, n_s,w_DE, fstar,fesc,nion,fx,flya;
    spline1dinterpolant interpolator;
};

struct q_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, t_cmb,w_DE;
    double h;
    spline1dinterpolant interpolator;
    spline1dinterpolant interpolator_Hf;
};

struct Tb_analytic_interpolator {
    double ombh2, omnuh2, omch2, hubble, t_cmb;
    spline1dinterpolant interpolator;
};

struct Fisher_return_pair {
    mat matrix;
    vector<vector<vector<string>>> matrix_indecies;
};

struct Ellipse {
    double a2, b2, theta;
};

