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
    //Fisher_return_pair finv = analyse.build_Fisher_inverse(keys, prefix.str(), "output/Fisher/");
    Fisher_return_pair finv;
    mat A = randu<mat>(3,3);
    A(0,0) = 1;
    A(0,1) = 0.25;
    A(0,2) = 0;
    A(1,0) = 0.25;
    A(1,1) = 0.25;
    A(1,2) = -0.025;
    A(2,0) = 0;
    A(2,1) = -0.025;
    A(2,2) = 0.01;
    cout << A << endl;
    vector<vector<vector<string>>> matrix_indecies;
    vector<string> params;
    vector<vector<string>> row;
    
    // row 1
    params.clear();
    params.push_back("ombh2");
    params.push_back("ombh2");
    row.push_back(params);
    
    params.clear();
    params.push_back("ombh2");
    params.push_back("omch2");
    row.push_back(params);
    
    params.clear();
    params.push_back("ombh2");
    params.push_back("hubble");
    row.push_back(params);
    matrix_indecies.push_back(row);
    row.clear();

    
    // row 2
    params.clear();
    params.push_back("omch2");
    params.push_back("ombh2");
    row.push_back(params);
    
    params.clear();
    params.push_back("omch2");
    params.push_back("omch2");
    row.push_back(params);
        
    params.clear();
    params.push_back("omch2");
    params.push_back("hubble");
    row.push_back(params);
    matrix_indecies.push_back(row);
    row.clear();

    // row 3
    params.clear();
    params.push_back("hubble");
    params.push_back("ombh2");
    row.push_back(params);
    
    params.clear();
    params.push_back("hubble");
    params.push_back("omch2");
    row.push_back(params);
    
    params.clear();
    params.push_back("hubble");
    params.push_back("hubble");
    row.push_back(params);
    matrix_indecies.push_back(row);
    row.clear();
        
    finv.matrix = A;
    finv.matrix_indecies = matrix_indecies;
    analyse.draw_error_ellipses(finv, keys, run_number, "output/Fisher/");

    return 0;
}


