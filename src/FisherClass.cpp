#include "FisherClass.hpp"
#include <fstream>
#include <time.h>

// IMPORTANT: Every parameter that I want to vary needs to be added to model_params_keys[]!!!

Fisher::Fisher(map<string, double> params, string Fl_filename, vector<string> params_keys_considered)
{
    cout << "... Beginning to build FisherClass ..." << endl;
    int Pk_index, Tb_index, q_index;
    CALC = new CosmoCalc(params, &Pk_index, &Tb_index, &q_index);
    this->current_params = CALC->give_current_params();
    this->fiducial_params = CALC->give_fiducial_params();
    kmin = this->fiducial_params["kmin"];
    kmax = this->fiducial_params["kmax"];
    model_params_keys = params_keys_considered;
    for (int i = 0; i < model_params_keys.size(); ++i) {
        string key = model_params_keys[i];
        if (current_params[key] == 0.0)
            var_params.insert(pair<string,double>(key,0.0001));
        else
            var_params.insert(pair<string,double>(key,current_params[key]/100));
    }
    noise = false;
    rsd = false;
    if (fiducial_params["noise"] == 1.0)
        noise = true;
    if (fiducial_params["rsd"] == 1.0)
        rsd = true;
    Fl_file.open(Fl_filename);
    cout << "... Fisher built ..." << endl;
}

Fisher::~Fisher()
{
    Fl_file.close();
    delete CALC;
}

void Fisher::update_Model(map<string, double> new_params, int *Pk_index, int *Tb_index, int *q_index)
{
    cout << "model is being updated" << endl;
    //generate params should not be necessary anymore...
    this->CALC->update_q(new_params, q_index);
    this->CALC->update_Pk_interpolator_direct(new_params, Pk_index);
    this->CALC->update_G21(new_params, Tb_index);
    //this->CALC->update_Tb_analytic(new_params, Tb_index);
}

mat Fisher::compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange)
{
    mat Cl = randu<mat>(krange.size(),krange.size());
    
    // Remove the lines below and the if/else statement when not reading/writing matrix
    stringstream matrix_filename;
    string suffix;
    if (noise && rsd)
        suffix = "rn";
    else if (noise && !rsd)
        suffix = "nrn";
    else if (!noise && rsd)
        suffix = "rnn";
    else if (!noise && !rsd)
        suffix = "nrnn";
    matrix_filename << "output/matrices/Cl_" << l << "_"<<\
        krange[0] << "_" << krange[krange.size()-1] << "_"<< krange.size() << "_"<<\
        fiducial_params["zmin"] << "_"<< fiducial_params["zmax"] << "_" << suffix << ".bin";
    if (check_file(matrix_filename.str()))
    {
        cout << "///reading matrix from file///" << endl;
        Cl = read_matrix(matrix_filename.str(),krange.size(),krange.size());
    }
    else
    {
        cout << "///calculating matrix///" << endl;

        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, Pk_index, Tb_index, q_index);
                Cl(i,j) = res;
                Cl(j,i) = res;
            }
        }
        
        //!!!!!!!!!!! this line also needs to be removed if not writing.
        write_matrix(Cl, matrix_filename.str());
    }

    return Cl;
}

double Fisher::Cl_derivative(int l, string param_key, double k1, double k2, int *Pk_index, int *Tb_index, int *q_index)
{
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];

    double f1,f2,f3,f4;
    bool do_calc = true;
    int index;

    working_params[param_key] = x + 2*h;
    /*
     * CAN'T do this in parallel.
     for (unsigned int i = 0; i < this->abcisses_done_simple.size(); ++i) { 
     if (this->abcisses_done_simple[i] == this->current_params[param_key]) {
     do_calc = false;
     index = i;
     break;
     }
     }
     */

    if (do_calc) {
        this->update_Model(working_params, Pk_index, Tb_index, q_index);
        f1 = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
        abcisses_done_simple.push_back(current_params[param_key]);
        derivs_calculated.push_back(f1); 
    } else 
        f1 = derivs_calculated[index];
    do_calc = true;

    working_params[param_key] = x + h;
    for (unsigned int i = 0; i < this->abcisses_done_simple.size(); ++i) { 
        if (this->abcisses_done_simple[i] == working_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }

    if (do_calc) {
        this->update_Model(working_params, Pk_index, Tb_index, q_index);
        f2 = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
        abcisses_done_simple.push_back(current_params[param_key]);
        derivs_calculated.push_back(f2); 
    } else 
        f2 = derivs_calculated[index];
    do_calc = true;


    working_params[param_key] = x - h;
    for (unsigned int i = 0; i < this->abcisses_done_simple.size(); ++i) { 
        if (this->abcisses_done_simple[i] == working_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }

    if (do_calc) {
        this->update_Model(working_params, Pk_index, Tb_index, q_index);
        f3 = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
        abcisses_done_simple.push_back(current_params[param_key]);
        derivs_calculated.push_back(f3); 
    } else 
        f3 = derivs_calculated[index];
    do_calc = true;

    working_params[param_key] = x - 2*h;
    for (unsigned int i = 0; i < this->abcisses_done_simple.size(); ++i) { 
        if (this->abcisses_done_simple[i] == working_params[param_key]) {
            do_calc = false;
            index = i;
            break;
        }
    }

    if (do_calc) {
        this->update_Model(working_params, Pk_index, Tb_index, q_index);
        f4 = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
        abcisses_done_simple.push_back(current_params[param_key]);
        derivs_calculated.push_back(f4); 
    } else 
        f4 = derivs_calculated[index];
    do_calc = true;

    working_params[param_key] = x;

    double num = -f1 + 8*f2 - 8*f3 + f4;
    double res = num /(12*h);

    return x*res;
}

double Fisher::Cl_loglog_derivative(int l, string param_key,\
        double k1, double k2, int *Pk_index, int *Tb_index, int *q_index)
{
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];

    double f1,f2,f3,f4;

    working_params[param_key] = x + 2*h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);
    f1 = log(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index));

    working_params[param_key] = x + h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);
    f2 = log(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index));

    working_params[param_key] = x - h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);
    f3 = log(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index));

    working_params[param_key] = x - 2*h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);
    f4 = log(this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index));

    working_params[param_key] = x;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);

    double num = -f1 + 8*f2 - 8*f3 + f4;
    double res = num /(12*h);

    return x*res;
}

mat Fisher::Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
        int *Tb_index, int *q_index, vector<double> krange)
{
    mat res = randu<mat>(krange.size(),krange.size());

    // Remove the lines below and the if/else statement when not reading/writing matrix
    stringstream matrix_filename;
    string suffix;
    if (rsd)
        suffix = "r";
    else 
        suffix = "nr";
    matrix_filename << "output/matrices/Cla_" << param_key << "_"<< l << "_" <<\
        krange[0] << "_" << krange[krange.size()-1] << "_"<< krange.size() << "_"<<\
        fiducial_params["zmin"] << "_"<< fiducial_params["zmax"] << "_" << suffix << ".bin";
    if (check_file(matrix_filename.str()))
    {
        cout << "///reading matrix from file///" << endl;
        res = read_matrix(matrix_filename.str(),krange.size(),krange.size());
    }
    else
    {
        cout << "///calculating matrix///" << endl;

        map<string,double> working_params = fiducial_params;
        double h = this->var_params[param_key];
        double x = working_params[param_key];
        mat f1matrix = randu<mat>(krange.size(),krange.size());
        mat f2matrix = randu<mat>(krange.size(),krange.size());
        mat f3matrix = randu<mat>(krange.size(),krange.size());
        mat f4matrix = randu<mat>(krange.size(),krange.size());
        working_params[param_key] = x + 2 * h;
        this->update_Model(working_params, Pk_index, Tb_index, q_index);
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
                f1matrix(i,j) = res;
                f1matrix(j,i) = res;
            }
        }

        working_params[param_key] = x + h;
        this->update_Model(working_params, Pk_index, Tb_index, q_index);
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
                f2matrix(i,j) = res;
                f2matrix(j,i) = res;
            }
        }

        working_params[param_key] = x - h;
        this->update_Model(working_params, Pk_index, Tb_index, q_index);
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
                f3matrix(i,j) = res;
                f3matrix(j,i) = res;
            }
        }

        working_params[param_key] = x - 2 * h;
        this->update_Model(working_params, Pk_index, Tb_index, q_index);
        for (unsigned int i = 0; i < krange.size(); ++i) {
            double k1 = krange[i];
            for (unsigned int j = i; j < krange.size(); ++j) {
                double k2 = krange[j];
                double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index, *Tb_index, *q_index);
                f4matrix(i,j) = res;
                f4matrix(j,i) = res;
            }
        }

        working_params[param_key] = x;
        this->update_Model(working_params, Pk_index, Tb_index, q_index);

        double num;
        for (unsigned int i = 0; i < krange.size(); ++i) {
            for (unsigned int j = 0; j < krange.size(); ++j) {
                num = -f1matrix(i,j) + 8*f2matrix(i,j) - 8*f3matrix(i,j) +\
                      f4matrix(i,j);
                num = num / (12.0 * h);    
                res(i,j) = num;
            }
        }


        //!!!!!!!!!!! this line also needs to be removed if not writing.
        write_matrix(res, matrix_filename.str());
    }
    return res;
}

double Fisher::compute_Fl(int l, string param_key1, string param_key2, int ksteps_Cl, double *cond_num,\
        int *Pk_index, int *Tb_index, int *q_index)
{
    vector<double> krange = give_kmodes(l, this->fiducial_params["kmax"], ksteps_Cl); 

    mat Cl = randu<mat>(krange.size(),krange.size());
    mat Cl_inv = Cl;

    cout << "... derivative matrix calulation started" << endl;
    mat Cl_a = this->Cl_derivative_matrix(l, param_key1, Pk_index, Tb_index, q_index, krange);
    mat Cl_b = randu<mat>(krange.size(),krange.size());
    if (param_key1 == param_key2)
        Cl_b = Cl_a;
    else
        Cl_b = this->Cl_derivative_matrix(l, param_key2, Pk_index, Tb_index, q_index, krange);

    cout << "-> The derivative matrices are done for l = " << l << endl;
    cout << "... The Cl and Cl_inv matrices will be calculated for l = " << l << endl;

    Cl = compute_Cl(l, *Pk_index, *Tb_index, *q_index, krange);
    *cond_num = cond(Cl);
    //Cl_inv = Cl.i();
    Cl_inv = pinv(Cl, 0.0001);
    cout << "-> Cl & Cl_inv are done for l = " << l << endl;

    mat product = Cl_a * Cl_inv;
    product = product * Cl_b;
    product = product * Cl_inv;

    return 0.5 * trace(product);
}

double Fisher::compute_Fl(int l, string param_key1, string param_key2, double kstepsize, double *cond_num,\
        int *Pk_index, int *Tb_index, int *q_index)
{
    vector<double> krange = give_kmodes(l, this->fiducial_params["kmax"], kstepsize); 

    mat Cl = randu<mat>(krange.size(),krange.size());
    mat Cl_inv = Cl;

    cout << "... derivative matrix calulation started" << endl;
    mat Cl_a = this->Cl_derivative_matrix(l, param_key1, Pk_index, Tb_index, q_index, krange);
    mat Cl_b = randu<mat>(krange.size(),krange.size());
    if (param_key1 == param_key2)
        Cl_b = Cl_a;
    else
        Cl_b = this->Cl_derivative_matrix(l, param_key2, Pk_index, Tb_index, q_index, krange);

    cout << "-> The derivative matrices are done for l = " << l << endl;
    cout << "... The Cl and Cl_inv matrices will be calculated for l = " << l << endl;

    Cl = compute_Cl(l, *Pk_index, *Tb_index, *q_index, krange);
    *cond_num = cond(Cl);
    //Cl_inv = Cl.i();
    Cl_inv = pinv(Cl, 0.0001);
    cout << "-> Cl & Cl_inv are done for l = " << l << endl;

    mat product = Cl_a * Cl_inv;
    product = product * Cl_b;
    product = product * Cl_inv;

    return 0.5 * trace(product);
}

void Fisher::initializer(string param_key, int *Pk_index, int *Tb_index, int *q_index)
{
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
    working_params[param_key] = x + 2 * h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x + h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x - h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x - 2 * h;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);

    working_params[param_key] = x;
    this->update_Model(working_params, Pk_index, Tb_index, q_index);
}

double Fisher::F(string param_key1, string param_key2)
{
    int Pk_index, Tb_index, q_index;
    int ksteps_Cl = 4;
    if (param_key1 == param_key2)
        initializer(param_key1, &Pk_index, &Tb_index, &q_index);
    else {
        initializer(param_key1, &Pk_index, &Tb_index, &q_index);
        initializer(param_key2, &Pk_index, &Tb_index, &q_index);
    }
    int l0 = 1000;
    int lmax = 2000;
    double sum = 0;
    // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!

    // The following line parallelizes the code
    // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
    // to define how many threads should be used.

#pragma omp parallel num_threads(1) private(Pk_index, Tb_index, q_index) 
    {
#pragma omp for reduction (+:sum)
        for (int l = l0; l < lmax; ++l) {
            stringstream ss, ss2, res;
            double cond_num = 0;
            ss << "Computation of Fl starts for l = " << l << "\n";
            cout << ss.str();
            double fl = this->compute_Fl(l, param_key1, param_key2, ksteps_Cl, &cond_num, &Pk_index,\
                    &Tb_index, &q_index);
            ss2 << "fl with l = " << l << " is: " << fl << "\n";
            cout << ss2.str();
            Fl_file << l << " " << fl << " " << cond_num << endl;
            sum += (2*l + 1) * fl;
        }
    }
    return sum;
}

double Fisher::F_fixed_kstepsize(string param_key1, string param_key2)
{
    int Pk_index, Tb_index, q_index;
    double kstepsize = 0.0178;
    if (param_key1 == param_key2)
        initializer(param_key1, &Pk_index, &Tb_index, &q_index);
    else {
        initializer(param_key1, &Pk_index, &Tb_index, &q_index);
        initializer(param_key2, &Pk_index, &Tb_index, &q_index);
    }
    int l0 = 1000;
    int lmax = 2000;
    double sum = 0;
    // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!

    // The following line parallelizes the code
    // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
    // to define how many threads should be used.

#pragma omp parallel num_threads(7) private(Pk_index, Tb_index, q_index) 
    {
#pragma omp for reduction (+:sum)
        for (int l = l0; l < lmax; ++l) {
            stringstream ss, ss2, res;
            double cond_num = 0;
            ss << "Computation of Fl starts for l = " << l << "\n";
            cout << ss.str();
            double fl = this->compute_Fl(l, param_key1, param_key2, kstepsize, &cond_num, &Pk_index,\
                    &Tb_index, &q_index);
            ss2 << "fl with l = " << l << " is: " << fl << "\n";
            cout << ss2.str();
            res << l << " " << fl << " " << cond_num << "\n";
            Fl_file << res.str() << endl;
            sum += (2*l + 1) * fl;
        }
    }
    return sum;
}

vector<double> Fisher::give_kmodes(int l, double k_max, int steps)
{
    double k_min; 
    
    /*if (l < 100) {
        k_min = (double)l/20000.0;
        if (k_min < 0.0001)
            k_min = 0.0001;
    } else if (l < 500) {
        k_min = (double)l/15000.0;
    } else {
        k_min = (double)l/10000.0;
    }
    */
    k_min = (double)l / CALC->r_interp(fiducial_params["zmax"]);
    double stepsize = (k_max - k_min)/(double)steps;
    vector<double> range;
    for (int i = 0; i <= steps; ++i)
    {
        range.push_back(k_min + i * stepsize); 
    }
    stringstream ss;
    ss << "The k_range is [" << k_min << "," << k_max << "] in " << steps+1 <<\
        " steps for l = " << l << ".\n";
    cout << ss.str();
    return range;
}

vector<double> Fisher::give_kmodes(int l, double k_max, double kstepsize)
{
    double k_min = (double)l / CALC->r_interp(fiducial_params["zmax"]);
    int steps = (k_max - k_min)/kstepsize + 1;
    vector<double> range;
    double new_max = 0;
    for (int i = 0; i <= steps; ++i)
    {
        new_max = k_min + i * kstepsize;
        range.push_back(new_max); 
    }
    stringstream ss;
    ss << "The k_range is [" << k_min << "," << new_max << "] in " << steps+1 <<\
        " steps for l = " << l << ".\n";
    cout << ss.str();
    return range;
}

void Fisher::Fl_varying_ksteps(int l, string param_key1, string param_key2, int min_ksteps_Cl,\
        int max_ksteps_Cl, int ksteps_spacing)
{
    // TODO: write a function that initializes all the interpolants, so that I can immediately call
    //       everything in parallel.
    int Pk_index, Tb_index, q_index;
    Fl_file << "# Contains Fl's for l = " << l << ", with number of kmodes sampled ranging from " <<\
        min_ksteps_Cl + 1 << " to " << max_ksteps_Cl + 1 <<\
        ". Also, the condition number of the covariance matrix is included in the third column." << endl;
    if (param_key1 == param_key2)
        initializer(param_key1, &Pk_index, &Tb_index, &q_index);
    else {
        initializer(param_key1, &Pk_index, &Tb_index, &q_index);
        initializer(param_key2, &Pk_index, &Tb_index, &q_index);
    }

    // The following line parallelizes the code
    // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
    // to define how many threads should be used.

#pragma omp parallel num_threads(7) private(Pk_index, Tb_index, q_index) 
    {
#pragma omp for 
        for (int ksteps_Cl = min_ksteps_Cl; ksteps_Cl <= max_ksteps_Cl; ksteps_Cl += ksteps_spacing) {
            double cond_num = 0;
            stringstream ss, ss2, res;
            ss << "Computation of Fl starts for # of kmodes = " << ksteps_Cl + 1 << "\n";
            cout << ss.str();
            double fl = this->compute_Fl(l, param_key1, param_key2, ksteps_Cl, &cond_num, &Pk_index,\
                    &Tb_index, &q_index);
            ss2 << "fl with number of k modes sampled = " << ksteps_Cl + 1 << " is: " << fl << "\n";
            cout << ss2.str();
            Fl_file << ksteps_Cl + 1 << " " << fl << " " << cond_num << endl; 
            //cout << res.str();
            //Fl_file << res.str();
        }
    }
}

void Fisher::Fl_varying_ksteps_smart(int l, string param_key1, string param_key2, int min_ksteps_Cl,\
        int max_ksteps_Cl)
{
    // TODO: write a function that initializes all the interpolants, so that I can immediately call
    //       everything in parallel.
    int Pk_index, Tb_index, q_index;
    Fl_file << "# Contains Fl's for l = " << l << ", with number of kmodes sampled ranging from " <<\
        min_ksteps_Cl + 1 << " to " << max_ksteps_Cl + 1 << " sampled in a smart way" <<\
        ". Also, the condition number of the covariance matrix is included in the third column." << endl;
    if (param_key1 == param_key2)
        initializer(param_key1, &Pk_index, &Tb_index, &q_index);
    else {
        initializer(param_key1, &Pk_index, &Tb_index, &q_index);
        initializer(param_key2, &Pk_index, &Tb_index, &q_index);
    }
    vector<int> num_modes;
    num_modes.push_back(2);
    int i = 0;
    while(num_modes[i] < max_ksteps_Cl)
    {
        num_modes.push_back(2*num_modes[i] - 1);
        i++;
    }
    // The following line parallelizes the code
    // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
    // to define how many threads should be used.

#pragma omp parallel num_threads(7) private(Pk_index, Tb_index, q_index) 
    {
#pragma omp for 
        for (int i = 0; i < num_modes.size(); i++) {
            double cond_num = 0;
            stringstream ss, ss2, res;
            int ksteps_Cl = num_modes[i] - 1;
            ss << "Computation of Fl starts for # of kmodes = " << ksteps_Cl + 1 << "\n";
            cout << ss.str();
            double fl = this->compute_Fl(l, param_key1, param_key2, ksteps_Cl, &cond_num, &Pk_index,\
                    &Tb_index, &q_index);
            ss2 << "fl with number of k modes sampled = " << ksteps_Cl + 1 << " is: " << fl << "\n";
            cout << ss2.str();
            res << ksteps_Cl + 1 << " " << fl << " " << cond_num << "\n";      
            Fl_file << res.str();
        }
    }
}
bool Fisher::check_file(string filename)
{
    //returns true if the file already exists.
    ifstream file(filename);
    bool res;
    if (file.good())
        res = true;
    else 
        res = false; 
    file.close();
    return res;
}

void Fisher::write_matrix(mat matrix, string filename)
{
    //Writing matrix to binary file.
    ofstream outfile;
    outfile.open(filename, ios::binary | ios::out);
    for (int i = 0; i < matrix.n_rows; i++)
    {
        for (int j = 0; j < matrix.n_cols; j++)
        {
            double x = matrix(i,j);
            outfile.write(reinterpret_cast<char*>(&x), sizeof(double));
        }
    }
    outfile.close();
}

mat Fisher::read_matrix(string filename, int n_rows, int n_cols)
{ 
    //To read matrix file written by Fisher::write_matrix(...).
    mat result;
    result = randu<mat>(n_rows, n_cols);
    ifstream infile;
    infile.open(filename, ios::binary | ios::in);
    double value = 0;
    for (int i = 0; i < n_rows; i++)
    {
        for (int j = 0; j < n_cols; j++)
        {
            infile.read(reinterpret_cast<char*>(&value), sizeof(double));
            result(i,j) = value;
        }
    }
    infile.close();
    return result;
}

string Fisher::update_runinfo(int lmin, int lmax,\
        int lstepsize, double kstepsize)
{
    string noise_incl = "N";
    string rsd_incl = "N";
    if (noise)
        noise_incl = "Y";
    if (rsd)
        rsd_incl = "Y";

    int run_number = 0;
    stringstream filename;
    filename << "output/Fisher/";
    //set_runnumber();
    ifstream run_info_file_in;
    run_info_file_in.open("output/Fisher/RUN_INFO.dat");
    vector<string> run_info;
    string first_line;
    getline(run_info_file_in, first_line);
    char last_ch = first_line.back();
    char slast_ch = first_line.at(first_line.length() - 2);
    if (slast_ch == '0'){
        stringstream buff;
        buff << last_ch;
        istringstream(buff.str()) >> run_number;
    } else {
        stringstream buff;

        buff << slast_ch << last_ch;
        istringstream(buff.str()) >> run_number;
    }
    run_number++;
    stringstream buffss;
    if (run_number < 10)
    {
        buffss << "HIGHEST RUN NUMBER COMPLETED: 0" << run_number;
        filename << "0" << run_number << "_Fisher_";
    } else {
        buffss << "HIGHEST RUN NUMBER COMPLETED: " << run_number;
        filename << run_number << "_Fisher_";
    }

    run_info.push_back(buffss.str());
    buffss.str("");
    while (run_info_file_in.good())
    {
        string buffs;
        getline(run_info_file_in, buffs);
        run_info.push_back(buffs);
    }
    run_info_file_in.close();
    run_info.pop_back(); 
    // append new information to run_info

    buffss.str("");
    buffss << " ";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << "### run number " << run_number << " ###";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << "# numerical parameter values used #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " zmin           = " << fiducial_params["zmin"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " zmax           = " << fiducial_params["zmax"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " zsteps         = " << fiducial_params["zsteps"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " Pk_steps       = " << fiducial_params["Pk_steps"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " kmin           = " << fiducial_params["kmin"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " kmax           = " << fiducial_params["kmax"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " k_stepsize     = " << fiducial_params["k_stepsize"];
    run_info.push_back(buffss.str());
    
    buffss.str("");
    buffss << "# physical parameter values used #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " ombh2          = " << fiducial_params["ombh2"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " omch2          = " << fiducial_params["omch2"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " omnuh2         = " << fiducial_params["omnuh2"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " omk            = " << fiducial_params["omk"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " hubble         = " << fiducial_params["hubble"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " T_CMB          = " << fiducial_params["T_CMB"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " sigma8         = " << fiducial_params["sigma8"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " 100*theta_s    = " << fiducial_params["100*theta_s"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " tau_reio       = " << fiducial_params["tau_reio"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " k_pivot        = " << fiducial_params["k_pivot"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " YHe            = " << fiducial_params["YHe"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " fstar          = " << fiducial_params["fstar"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " fesc           = " << fiducial_params["fesc"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " nion           = " << fiducial_params["nion"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " flya           = " << fiducial_params["flya"];
    run_info.push_back(buffss.str());

    buffss.str("");
    buffss << "# g21 flags #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " popflag        = " << fiducial_params["popflag"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " xrayflag       = " << fiducial_params["xrayflag"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lyaxrayflag    = " << fiducial_params["lyaxrayflag"];
    run_info.push_back(buffss.str());
    
    buffss.str("");
    buffss << "# Noise #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " noise included = " << noise_incl;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " Ae             = " << fiducial_params["Ae"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " df             = " << fiducial_params["df"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " Tsys           = " << fiducial_params["Tsys"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " fcover         = " << fiducial_params["fcover"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lmax_noise     = " << fiducial_params["lmax_noise"];
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " tau_noise      = " << fiducial_params["tau_noise"];
    run_info.push_back(buffss.str());

    buffss.str("");
    buffss << "# Other #";
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " rsd included   = " << rsd_incl;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lmin           = " << lmin;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lmax           = " << lmax;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " lstepsize      = " << lstepsize;
    run_info.push_back(buffss.str());
    buffss.str("");
    buffss << " k_mode stepsize= " << kstepsize;
    run_info.push_back(buffss.str());
    
    //then write it all to file.
    ofstream run_info_file_out;
    run_info_file_out.open("output/Fisher/RUN_INFO.dat");
    for (int i = 0; i < run_info.size(); i++){
        run_info_file_out << run_info[i] << endl;
    }
    run_info_file_out.close();

    return filename.str();
}

double Fisher::F_fixed_kstepsize(int lmin, int lmax, int n_points_per_thread, int n_threads)
{
    int lsteps = n_points_per_thread * n_threads;
    int lstepsize = ((double)(lmax-lmin))/(double)lsteps;
    double kstepsize = 0.0178;
    string filename_prefix = update_runinfo(lmin, lmax, lstepsize, kstepsize);
    stringstream filename;
    filename << filename_prefix;

    // now compute F_ab's (symmetric hence = F_ba's)
    for (int i = 0; i < model_params_keys.size(); i++) {
        for (int j = i; j < model_params_keys.size(); j++) {
            filename.str("");
            string param_key1 = model_params_keys[i];
            string param_key2 = model_params_keys[j];
            filename << filename_prefix << param_key1 << "_" << param_key2 << ".dat";
            ofstream outfile;
            outfile.open(filename.str());
            int Pk_index, Tb_index, q_index;
            double kstepsize = 0.0178;
            if (param_key1 == param_key2)
                initializer(param_key1, &Pk_index, &Tb_index, &q_index);
            else {
                initializer(param_key1, &Pk_index, &Tb_index, &q_index);
                initializer(param_key2, &Pk_index, &Tb_index, &q_index);
            }
            double sum = 0;
            // IMPORTANT! l has to start at 1 since Nl_bar has j_(l-1) in it!

            // The following line parallelizes the code
            // use #pragma omp parallel num_threads(4) private(Pk_index, Tb_index, q_index) 
            // to define how many threads should be used.

            #pragma omp parallel num_threads(n_threads) private(Pk_index, Tb_index, q_index) 
            {
                #pragma omp for reduction (+:sum)
                for (int k = 1; k <= lsteps; ++k) {
                    int m;
                    if (k == lsteps)
                        m = lsteps;
                    else
                        m = ((k-1)*n_threads) % (lsteps - 1) + 1;
                    int l = lmin + m * lstepsize;
                    stringstream ss, ss2, res;
                    double cond_num = 0;
                    ss << "Computation of Fl starts for l = " << l << "\n";
                    cout << ss.str();
                    double fl = this->compute_Fl(l, param_key1, param_key2, kstepsize,\
                            &cond_num, &Pk_index, &Tb_index, &q_index);
                    ss2 << "fl with l = " << l << " is: " << fl << "\n";
                    cout << ss2.str();
                    res << l << " " << fl << " " << cond_num << "\n";
                    outfile << res.str() << endl;
                    sum += (2*l + 1) * fl;
                }
            }
            outfile.close();
        }
    }

    return 0;
}


