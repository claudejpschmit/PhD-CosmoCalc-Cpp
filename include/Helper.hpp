#pragma once

#include "stdafx.h"
#include "interpolation.h"

using namespace std;
using namespace alglib;

struct Pk_interpolator {
    double ombh2, omnuh2, omch2, omk, hubble;
    spline2dinterpolant interpolator;
};
