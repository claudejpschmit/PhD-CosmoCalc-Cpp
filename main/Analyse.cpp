#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <armadillo>
#include "stdafx.h"
#include "interpolation.h"
#include "Analyser.hpp"

using namespace std;
using namespace arma;
using namespace alglib;

int main(int argc, char* argv[])
{
    vector<string> keys = {"ombh2", "omch2", "hubble","fesc", "fstar"};
    Analyser analyse;
    Fisher_return_pair finv = analyse.build_Fisher_inverse(keys, "11", "output/Fisher/");
    analyse.draw_error_ellipses(finv, keys, 11, "output/Fisher/");

    return 0;
}


