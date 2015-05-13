#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include "CosmoBasis.hpp"

using namespace std;

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    map<string,double> params;
     
    params["O"] = 0.02;
    params["ombh2"] = 0.02;
    cout << params["O"] << endl;
    CosmoBasis base(params);
    base.show_params(); 
    //cout << 1.34 * pow(10,2) << endl;
    return 0;
}
