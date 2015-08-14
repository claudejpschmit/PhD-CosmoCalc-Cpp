#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include "CosmoBasis.hpp"
#include <fstream>
#include <sys/stat.h>
#include "Integrator.hpp"
#include "CosmologyCalculatorClass.hpp"
#include "CosmologyWriterClass.hpp"
#include <time.h>
#include "FisherClass.hpp"
#include <armadillo>
#include "stdafx.h"
#include "interpolation.h"
#include "CAMB_interface.hpp"
#include "Global21cmInterface.hpp"
#include <boost/math/special_functions/bessel.hpp>
#include "SanityChecker.hpp"
#include "LevinIntegrator.hpp"

using namespace std;
using namespace arma;
using namespace alglib;

int main(int argc, char* argv[])
{
    map<string,double> params;
    
    int Pk_index = 0;
    int Tb_index = 0;
    int q_index = 0;
    clock_t t1,t2;
    float d1;
    Fisher fish(params,"Fisher_norsd.dat");
    t1 = clock();
    cout << "The result is = " << fish.F("ombh2", "ombh2")<< endl;
    t2 = clock();
    d1 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    cout << "time1 : " << d1 << endl;

    Fisher fish2(params,"Fisher_norsd.dat");

    t1 = clock();
    cout << "The result is = " << fish2.new_F("ombh2", "ombh2")<< endl;
    t2 = clock();
    d1 = ((float)t2 - (float)t1)/CLOCKS_PER_SEC;
    cout << "New Time = " << d1 << endl;


    return 0;
}


