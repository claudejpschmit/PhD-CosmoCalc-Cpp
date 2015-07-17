#include "FisherClass.hpp"
#include <fstream>
#include <time.h>

Fisher::Fisher(map<string, double> params, string Fl_filename)
{
    cout << "... Beginning to build FisherClass ..." << endl;
    CALC = new CosmoCalc(params);
    this->current_params = CALC->give_current_params();
    this->fiducial_params = CALC->give_fiducial_params();
    kmin = this->fiducial_params["kmin"];
    kmax = this->fiducial_params["kmax"];
    string model_params_keys[] = {"ombh2", "omch2", "omnuh2", "omk", "hubble"};
    for (int i = 0; i < 5; ++i) {
        string key = model_params_keys[i];
        if (current_params[key] == 0.0)
            var_params.insert(pair<string,double>(key,0.0001));
        else
            var_params.insert(pair<string,double>(key,current_params[key]/100));
    }
    
    Fl_file.open(Fl_filename);

    cout << "... Fisher built ..." << endl;
}

Fisher::~Fisher()
{
    Fl_file.close();
    delete CALC;
}

void Fisher::update_Model(map<string, double> new_params)
{
    this->CALC->generate_params(new_params);
    this->CALC->update_q();
    this->CALC->update_q_prime();
    this->CALC->update_Pk_interpolator_direct(new_params);
    this->CALC->update_G21(new_params);
    //this->CALC->updateClass(new_params);
}

void Fisher::compute_Cl(int l)
{
    for (unsigned int i = 0; i < this->krange.size(); ++i) {
        double k1 = this->krange[i];
        for (unsigned int j = i; j < this->krange.size(); ++j) {
            double k2 = this->krange[j];
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax);
            Cl(i,j) = res;
            Cl(j,i) = res;
        }
    }
}

void Fisher::show_Cl_mat()
{
    cout << Cl << endl;
}

void Fisher::show_Cl_inv_mat()
{
    cout << Cl_inv << endl;
}

void Fisher::compute_Cl_inv()
{
    this->Cl_inv = this->Cl.i();
}

double Fisher::Cl_derivative(int l, string param_key, double k1, double k2)
{
    double h = this->var_params[param_key];
    double x = this->current_params[param_key];

    double f1,f2,f3,f4;
    bool do_calc = true;
    int index;

    this->current_params[param_key] = x + 2*h;
    for (unsigned int i = 0; i < this->abcisses_done_simple.size(); ++i) { 
        if (this->abcisses_done_simple[i] == this->current_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }

    if (do_calc) {
        this->update_Model(this->current_params);
        f1 = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax);
        abcisses_done_simple.push_back(current_params[param_key]);
        derivs_calculated.push_back(f1); 
    } else 
        f1 = derivs_calculated[index];
    do_calc = true;

    this->current_params[param_key] = x + h;
    for (unsigned int i = 0; i < this->abcisses_done_simple.size(); ++i) { 
        if (this->abcisses_done_simple[i] == this->current_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }

    if (do_calc) {
        this->update_Model(this->current_params);
        f2 = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax);
        abcisses_done_simple.push_back(current_params[param_key]);
        derivs_calculated.push_back(f2); 
    } else 
        f2 = derivs_calculated[index];
    do_calc = true;


    this->current_params[param_key] = x - h;
    for (unsigned int i = 0; i < this->abcisses_done_simple.size(); ++i) { 
        if (this->abcisses_done_simple[i] == this->current_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }

    if (do_calc) {
        this->update_Model(this->current_params);
        f3 = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax);
        abcisses_done_simple.push_back(current_params[param_key]);
        derivs_calculated.push_back(f3); 
    } else 
        f3 = derivs_calculated[index];
    do_calc = true;

    this->current_params[param_key] = x - 2*h;
    for (unsigned int i = 0; i < this->abcisses_done_simple.size(); ++i) { 
        if (this->abcisses_done_simple[i] == this->current_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }

    if (do_calc) {
        this->update_Model(this->current_params);
        f4 = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax);
        abcisses_done_simple.push_back(current_params[param_key]);
        derivs_calculated.push_back(f4); 
    } else 
        f4 = derivs_calculated[index];
    do_calc = true;

    this->current_params[param_key] = x;

    double num = -f1 + 8*f2 - 8*f3 + f4;
    double res = num /(12*h);

    return x*res;
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
    for (unsigned int i = 0; i < this->abcisses_done.size(); ++i) { 
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
    for (unsigned int i = 0; i < this->abcisses_done.size(); ++i) { 
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
    for (unsigned int i = 0; i < this->abcisses_done.size(); ++i) { 
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
    for (unsigned int i = 0; i < this->abcisses_done.size(); ++i) { 
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
    double h = this->var_params[param_key];
    double x = this->current_params[param_key];

    vector<vector<double>> res, f1matrix, f2matrix, f3matrix, f4matrix;
    vector<double> row;
    this->current_params[param_key] = x + 2 * h;
    this->update_Model(this->current_params);
    for (unsigned int i = 0; i < this->krange.size(); ++i) {
        double k1 = this->krange[i];
        row.clear();
        for (unsigned int j = 0; j < this->krange.size(); ++j) {
            double k2 = this->krange[j];
            double cl = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax);
            row.push_back(cl);
        }
        f1matrix.push_back(row);
    }
    this->current_params[param_key] = x + h;
    this->update_Model(this->current_params);
    for (unsigned int i = 0; i < this->krange.size(); ++i) {
        double k1 = this->krange[i];
        row.clear();
        for (unsigned int j = 0; j < this->krange.size(); ++j) {
            double k2 = this->krange[j];
            row.push_back(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax));
        }
        f2matrix.push_back(row);
    }

    this->current_params[param_key] = x - h;
    this->update_Model(this->current_params);
    for (unsigned int i = 0; i < this->krange.size(); ++i) {
        double k1 = this->krange[i];
        row.clear();
        for (unsigned int j = 0; j < this->krange.size(); ++j) {
            double k2 = this->krange[j];
            row.push_back(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax));
        }
        f3matrix.push_back(row);
    }

    this->current_params[param_key] = x - 2 * h;
    this->update_Model(this->current_params);
    for (unsigned int i = 0; i < this->krange.size(); ++i) {
        double k1 = this->krange[i];
        row.clear();
        for (unsigned int j = 0; j < this->krange.size(); ++j) {
            double k2 = this->krange[j];
            row.push_back(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax));
        }
        f4matrix.push_back(row);
    }

    this->current_params[param_key] = x;
    this->update_Model(this->current_params);

    double num;
    for (unsigned int i = 0; i < this->krange.size(); ++i) {
        row.clear();
        for (unsigned int j = 0; j < this->krange.size(); ++j) {
            num = -f1matrix[i][j] + 8*f2matrix[i][j] - 8*f3matrix[i][j] +\
                  f4matrix[i][j];
            num = num / (12.0 * h);    
            row.push_back(num);
        }
        res.push_back(row);
    }

    return res;
}

double Fisher::compute_Fl(int l, string param_key1, string param_key2)
{
    vector<vector<double>> Cl_alpha, Cl_beta;
    //This determines the size of the Cl matrices.
    int ksteps_Cl = 4;
    krange = give_kmodes(l, this->fiducial_params["kmax"], ksteps_Cl); 
    //for (int i = 0; i < krange.size(); i++)
    //    cout << krange[i] << " ";
    //cout << endl;

    Cl = randu<mat>(krange.size(),krange.size());
    Cl_inv = Cl;

    cout << "... derivative matrix calulation started" << endl;
    Cl_alpha = this->Cl_derivative_matrix(l, param_key1);
    if (param_key1 == param_key2)
        Cl_beta = Cl_alpha;
    else
        Cl_beta = this->Cl_derivative_matrix(l, param_key2);

    cout << "-> The derivative matrices are done for l = " << l << endl;
    cout << "... The Cl and Cl_inv matrices will be calculated for l = " << l << endl;

    this->compute_Cl(l);
    //cout << Cl << endl;

    this->compute_Cl_inv();
    cout << "-> Cl & Cl_inv are done for l = " << l << endl;
    mat Cl_a, Cl_b;
    Cl_a = randu<mat> (this->krange.size(), this->krange.size());
    Cl_b = Cl_a;
    for (unsigned int i = 0; i < this->krange.size(); ++i) {
        for (unsigned int j = 0; j < this->krange.size(); ++j) {
            Cl_a(i,j) = Cl_alpha[i][j];
            Cl_b(i,j) = Cl_beta[i][j];

        }
    } 
   
    //cout << "with det = " << det(Cl) << endl;
    //cout << Cl_a << endl;
    //cout << Cl_inv << endl;
    //cout << Cl_b << endl;
    mat product = Cl_a * this->Cl_inv;
    //cout << product << endl;
    product = product * Cl_b;
    //cout << product << endl;
    product = product * this->Cl_inv;
    //cout << product << endl;

    return 0.5 * trace(product);
}

double Fisher::F(string param_key1, string param_key2)
{
    int lmax = 5000;
    double sum = 0;
    // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!
    for (int l = 1; l <= lmax; ++l) {
        cout << "Computation of Fl starts for l = " << l << endl;
        double fl = this->compute_Fl(l, param_key1, param_key2);
        cout << "fl with l = " << l << " is: " << fl << endl;
        Fl_file << l << " " << fl << endl;
        sum += (2*l + 1) * fl;
    }
    return sum;
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
    file << "# Column 3: Real runtime per point in second." << endl;

    abcisses_done.clear();
    logderivs_calculated.clear();
    double y1;
    double stepsize = (stepsize_high - stepsize_low)/(double)steps;
    clock_t t1, t2;

    for (int i = 0; i <= steps; i++) {
        t1 = clock();

        this->var_params[param_key] = stepsize_low + i * stepsize;
        y1 = this->Cl_loglog_derivative(l, param_key, k1, k2);

        t2 = clock();
        float diff ((float)t2 - (float)t1);

        file << this->var_params[param_key] << " " << y1 << " " << diff/CLOCKS_PER_SEC << endl;
        cout << i << "th point written" << endl;
        cout << "Runtime for the " << i << "th point was: " << diff/CLOCKS_PER_SEC  << " seconds." << endl;
    }

    file.close(); 
}

vector<double> Fisher::give_kmodes(int l, double k_max, int steps)
{
    double k_min; 
    if (l < 100) {
        k_min = (double)l/20000.0;
        if (k_min < 0.0001)
            k_min = 0.0001;
    } else if (l < 500) {
        k_min = (double)l/15000.0;
    } else {
        k_min = (double)l/10000.0;
    }
    double stepsize = (k_max - k_min)/(double)steps;
    vector<double> range;
    for (int i = 0; i <= steps; ++i)
    {
        range.push_back(k_min + i * stepsize); 
    }
    return range;
}
