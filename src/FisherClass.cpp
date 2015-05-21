#include "FisherClass.hpp"

Fisher::Fisher(map<string, double> params)
    :
        kmin(0.001),
        kmax(5)
{
    CALC = new CosmoCalc(params);
    
    this->current_params = CALC->give_current_params();
    this->fiducial_params = CALC->give_fiducial_params();

    string model_params_keys[] = {"ombh2", "omch2", "omnuh2", "omk", "hubble"};
    for (int i = 0; i < 5; ++i) {
        string key = model_params_keys[i];
        if (current_params[key] == 0.0)
            var_params.insert(pair<string,double>(key,0.0001));
        else
            var_params.insert(pair<string,double>(key,current_params[key]/100));
    }
    int ksteps = 100;
    double kstepsize = (kmax - kmin)/(double)ksteps;
    for (int n = 0; n <= ksteps; ++n) 
        krange.push_back(kmin + n * kstepsize);
    cout << "Fisher initialized" << endl;

}

Fisher::~Fisher()
{
    delete CALC;
}

