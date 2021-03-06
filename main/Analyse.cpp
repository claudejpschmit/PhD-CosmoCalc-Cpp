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
    int run_number;
    vector<string> keys;
    if (argc < 2){
        run_number = 11;
        keys = {"ombh2", "omch2", "hubble","fesc","w_DE"};
    }
    else {
        run_number = atoi(argv[1]);
        for (int i = 2; i < argc; i++)
            keys.push_back(argv[i]);
    }
    stringstream prefix;
    if (run_number < 10)
        prefix << "0" << run_number;
    else 
        prefix << run_number;

    Analyser analyse;
    Fisher_return_pair finv = analyse.build_Fisher_inverse(keys, prefix.str(), "output/Fisher/");
    analyse.draw_error_ellipses(finv, keys, run_number, "output/Fisher/");

    return 0;
}


