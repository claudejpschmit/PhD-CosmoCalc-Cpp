#pragma once

#include "CosmoBasis.hpp"
#include <string>
#include <map>
#include <cmath>


using namespace std;
    
class CosmoCalc : public CosmoBasis
{
    public:

        // ------------ Functions -------------- //
 
        CosmoCalc(map<string, double> params);

        ~CosmoCalc();

};
