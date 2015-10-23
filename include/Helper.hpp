#pragma once

#include "stdafx.h"
#include "interpolation.h"
#include <vector>
#include <armadillo>
#include <string>
using namespace std;
using namespace alglib;
using namespace arma;

/**
 * Pk_interpolator includes a value for each parameter affecting the power
 * spectrum Pk(z), and an interpolation object corresponding to those 
 * parameter values.
 */
struct Pk_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, tcmb, w_DE;
    spline2dinterpolant interpolator;
};

/**
 * Tb_interpolator includes a value for each parameter affecting the 21cm
 * brightness temperature fluctuations Tb(z), and an interpolation object 
 * corresponding to those parameter values. This version is to be used
 * with the Global21cm Code.
 */
struct Tb_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, s8, T_CMB, n_s, w_DE,\
        fstar, fesc, nion, fx, flya;
    spline1dinterpolant interpolator;
};

/**
 * Tb_interpolator includes a value for each parameter affecting the 21cm
 * brightness temperature fluctuations Tb(z), and an interpolation object 
 * corresponding to those parameter values. This version is to be used
 * with the Global21cm Code.
 */
struct Tb_interpolator_ares {
    double ombh2, omnuh2, omch2, omk, hubble, s8, T_CMB, n_s, w_DE,\
        fstar, fesc, nion, fX, Tmin, Nlw, cX, HeByMass;
    spline1dinterpolant interpolator;
};

/**
 * q_interpolator includes a value for each parameter affecting the comoving
 * distance r(z) = q(z), and an interpolation object corresponding to those 
 * parameter values for both the distance as well as the Hubble function.
 */
struct q_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, t_cmb, w_DE;
    double h;
    spline1dinterpolant interpolator;
    spline1dinterpolant interpolator_Hf;
};

/**
 * Tb_analytic_interpolator includes a value for each parameter affecting the
 * 21cm brightness temperature fluctuations Tb(z), and an interpolation object 
 * corresponding to those parameter values. This version is for the analytic
 * determination of Tb(z).
 */
struct Tb_analytic_interpolator {
    double ombh2, omnuh2, omch2, hubble, t_cmb;
    spline1dinterpolant interpolator;
};

/**
 * A Fisher_return pair includes a matrix of Fisher values, as well as a
 * matrix of parameter key pairs so to identify which matrix elements 
 * correspond to which parameters.
 */

struct Fisher_return_pair {
    mat matrix;
    vector<vector<vector<string>>> matrix_indecies;
};

/**
 * The ellipse structure includes all necessary information to define an 
 * ellipse.
 */
struct Ellipse {
    double a2, b2, theta, cx, cy;
};

