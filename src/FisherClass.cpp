#include "FisherClass.hpp"
#include <fstream>

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

void Fisher::update_Model(map<string, double> new_params)
{
    this->CALC->generate_params(new_params);
    this->CALC->update_q();
    this->CALC->updateClass(new_params);
}

void Fisher::compute_Cl(int l)
{

}

void Fisher::compute_Cl_inv()
{

}

double Fisher::Cl_derivative(int l, string param_key, double k1, double k2)
{
    return 0;
}

double Fisher::Cl_loglog_derivative(int l, string param_key,\
                                    double k1, double k2)
{
    double h = this->var_params[param_key];
    double x = this->current_params[param_key];
    double f1,f2,f3,f4;
    bool do_calc = true;
    int index;

    this->current_params[param_key] = x + 2*h;
    for (int i = 0; i < this->abcisses_done.size(); ++i) { 
        if (this->abcisses_done[i] == this->current_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }
    
    if (do_calc) {
        this->update_Model(this->current_params);
        f1 = log(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax));
        abcisses_done.push_back(current_params[param_key]);
        logderivs_calculated.push_back(f1); 
    } else 
        f1 = logderivs_calculated[index];
    do_calc = true;

    this->current_params[param_key] = x + h;
    for (int i = 0; i < this->abcisses_done.size(); ++i) { 
        cout << this->abcisses_done[i] << " =? " <<\
                this->current_params[param_key] << endl;

        if (this->abcisses_done[i] == this->current_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }
    
    if (do_calc) {
        this->update_Model(this->current_params);
        f2 = log(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax));
        abcisses_done.push_back(current_params[param_key]);
        logderivs_calculated.push_back(f2); 
    } else 
        f2 = logderivs_calculated[index];
    do_calc = true;

    
    this->current_params[param_key] = x - h;
    for (int i = 0; i < this->abcisses_done.size(); ++i) { 
        if (this->abcisses_done[i] == this->current_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }
    
    if (do_calc) {
        this->update_Model(this->current_params);
        f3 = log(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax));
        abcisses_done.push_back(current_params[param_key]);
        logderivs_calculated.push_back(f3); 
    } else 
        f3 = logderivs_calculated[index];
    do_calc = true;

    this->current_params[param_key] = x - 2*h;
    for (int i = 0; i < this->abcisses_done.size(); ++i) { 
        if (this->abcisses_done[i] == this->current_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }
    
    if (do_calc) {
        this->update_Model(this->current_params);
        f4 = log(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax));
        abcisses_done.push_back(current_params[param_key]);
        logderivs_calculated.push_back(f4); 
    } else 
        f4 = logderivs_calculated[index];
    do_calc = true;

    this->current_params[param_key] = x;

    double num = -f1 + 8*f2 - 8*f3 + f4;
    double res = num /(12*h);
    
    return x*res;
}

vector<vector<double>> Fisher::Cl_derivative_matrix(int l, string param_key)
{
    vector<vector<double>> res;
    return res;
}

double Fisher::compute_Fl(int l, string param_key1, string param_key2)
{
    return 0;
}

double Fisher::F(string param_key1, string param_key2)
{
    return 0;
}

void Fisher::write_logder(string param_key, double param_val,\ 
        double stepsize_low, double stepsize_high,\
        int steps, int l, double k1, double k2,\
        string suffix)
{
    if (this->current_params[param_key] != param_val){
        this->current_params[param_key] = param_val;
        update_Model(this->current_params);
    }

    ofstream file;
    string filename = "output/fisherlogder_" + param_key + suffix + ".dat";

    file.open(filename);
    file << "# This file contains data log derivative of the Cl's wrt " <<\
         param_key << " vs derivative stepsize." << endl;
    file << "# Column 1: derivative stepsize h" << endl;
    file << "# Column 2: Logderivative dlnCl/dx." << endl;
  
    abcisses_done.clear();
    logderivs_calculated.clear();
    double y1;
    double stepsize = (stepsize_high - stepsize_low)/(double)steps;
    for (int i = 0; i <= steps; i++) {
        
        this->var_params[param_key] = stepsize_low + i * stepsize;
        y1 = this->Cl_loglog_derivative(l, param_key, k1, k2);
       
        file << this->var_params[param_key] << " " << y1 << endl;
        cout << i << "th point written" << endl;
    }
    file.close();
 
}

