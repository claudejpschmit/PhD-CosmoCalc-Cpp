#pragma once

#include "stdafx.h"
#include "interpolation.h"
#include <vector>

using namespace std;
using namespace alglib;

struct Pk_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble;
    spline2dinterpolant interpolator;
};

struct Tb_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble, s8, T_CMB, n_s, fstar,fesc,nion,fx,flya;
    spline1dinterpolant interpolator;
};


