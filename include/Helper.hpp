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
    double ombh2, omnuh2, omch2, omk, hubble, tcmb;
    spline2dinterpolant interpolator;
};

struct Tb_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, s8, T_CMB, n_s, fstar,fesc,nion,fx,flya;
    spline1dinterpolant interpolator;
};

struct q_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, t_cmb;
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

