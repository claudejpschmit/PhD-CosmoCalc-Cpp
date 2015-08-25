#include "FisherClass.hpp"
#include <fstream>
#include <time.h>

Fisher::Fisher(map<string, double> params, string Fl_filename)
{
    cout << "... Beginning to build FisherClass ..." << endl;
    int Pk_index, Tb_index, q_index;
    CALC = new CosmoCalc(params, &Pk_index, &Tb_index, &q_index);
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

void Fisher::update_Model(map<string, double> new_params, int *Pk_index, int *Tb_index, int *q_index)
{
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
    string suffix = "nrn";
    matrix_filename << "output/matrices/Cl_" << l << "_"<<\
        krange[krange.size()-1] << "_"<< krange.size() << "_"<<\
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
mat Fisher::compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange,\
        spline1dinterpolant bessels)
{
    mat Cl = randu<mat>(krange.size(),krange.size());
    for (unsigned int i = 0; i < krange.size(); ++i) {
        double k1 = krange[i];
        for (unsigned int j = i; j < krange.size(); ++j) {
            double k2 = krange[j];
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, Pk_index, Tb_index,\
                    q_index, bessels);
            Cl(i,j) = res;
            Cl(j,i) = res;
        }
    }

    return Cl;
}
mat Fisher::compute_Cl(int l, int Pk_index, int Tb_index, int q_index, vector<double> krange,\
        spline1dinterpolant bessels, spline1dinterpolant bessels_lminus1)
{
    mat Cl = randu<mat>(krange.size(),krange.size());
    for (unsigned int i = 0; i < krange.size(); ++i) {
        double k1 = krange[i];
        for (unsigned int j = i; j < krange.size(); ++j) {
            double k2 = krange[j];
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, Pk_index, Tb_index,\
                    q_index, bessels, bessels_lminus1);
            Cl(i,j) = res;
            Cl(j,i) = res;
        }
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
    string suffix = "nrn";
    matrix_filename << "output/matrices/Cla_" << param_key << "_"<< l << "_"<<\
        krange[krange.size()-1] << "_"<< krange.size() << "_"<<\
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

mat Fisher::Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
        int *Tb_index, int *q_index, vector<double> krange, spline1dinterpolant bessels)
{
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
    mat res = randu<mat>(krange.size(),krange.size());
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
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index,\
                    *Tb_index, *q_index, bessels);

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
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index,\
                    *Tb_index, *q_index, bessels);

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
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index,\
                    *Tb_index, *q_index, bessels);

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
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index,\
                    *Tb_index, *q_index, bessels);

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

    return res;
}

mat Fisher::Cl_derivative_matrix(int l, string param_key, int *Pk_index,\
        int *Tb_index, int *q_index, vector<double> krange, spline1dinterpolant bessels,\
        spline1dinterpolant bessels_lminus1)
{
    map<string,double> working_params = fiducial_params;
    double h = this->var_params[param_key];
    double x = working_params[param_key];
    mat res = randu<mat>(krange.size(),krange.size());
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
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index,\
                    *Tb_index, *q_index, bessels, bessels_lminus1);

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
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index,\
                    *Tb_index, *q_index, bessels, bessels_lminus1);

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
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index,\
                    *Tb_index, *q_index, bessels, bessels_lminus1);

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
            double res = this->CALC->Cl(l, k1, k2, this->kmin, this->kmax, *Pk_index,\
                    *Tb_index, *q_index, bessels, bessels_lminus1);

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
    Cl_inv = Cl.i();

    cout << "-> Cl & Cl_inv are done for l = " << l << endl;

    mat product = Cl_a * Cl_inv;
    product = product * Cl_b;
    product = product * Cl_inv;

    return 0.5 * trace(product);
}

double Fisher::compute_Fl(int l, string param_key1, string param_key2, int ksteps_Cl, int *Pk_index,\
        int *Tb_index, int *q_index, spline1dinterpolant bessels)
{
    vector<double> krange = give_kmodes(l, this->fiducial_params["kmax"], ksteps_Cl); 

    mat Cl = randu<mat>(krange.size(),krange.size());
    mat Cl_inv = Cl;

    cout << "... derivative matrix calulation started" << endl;
    mat Cl_a = this->Cl_derivative_matrix(l, param_key1, Pk_index, Tb_index, q_index, krange, bessels);
    mat Cl_b = randu<mat>(krange.size(),krange.size());
    if (param_key1 == param_key2)
        Cl_b = Cl_a;
    else
        Cl_b = this->Cl_derivative_matrix(l, param_key2, Pk_index, Tb_index, q_index, krange, bessels);

    cout << "-> The derivative matrices are done for l = " << l << endl;
    cout << "... The Cl and Cl_inv matrices will be calculated for l = " << l << endl;

    Cl = compute_Cl(l, *Pk_index, *Tb_index, *q_index, krange, bessels);
    Cl_inv = Cl.i();

    cout << "-> Cl & Cl_inv are done for l = " << l << endl;

    mat product = Cl_a * Cl_inv;
    product = product * Cl_b;
    product = product * Cl_inv;

    return 0.5 * trace(product);
}

double Fisher::compute_Fl(int l, string param_key1, string param_key2, int ksteps_Cl, int *Pk_index,\
        int *Tb_index, int *q_index, spline1dinterpolant bessels, spline1dinterpolant bessels_lminus1)
{
    vector<double> krange = give_kmodes(l, this->fiducial_params["kmax"], ksteps_Cl); 

    mat Cl = randu<mat>(krange.size(),krange.size());
    mat Cl_inv = Cl;

    cout << "... derivative matrix calulation started" << endl;
    mat Cl_a = this->Cl_derivative_matrix(l, param_key1, Pk_index, Tb_index,\
            q_index, krange, bessels, bessels_lminus1);
    mat Cl_b = randu<mat>(krange.size(),krange.size());
    if (param_key1 == param_key2)
        Cl_b = Cl_a;
    else
        Cl_b = this->Cl_derivative_matrix(l, param_key2, Pk_index, Tb_index,\
                q_index, krange, bessels, bessels_lminus1);

    cout << "-> The derivative matrices are done for l = " << l << endl;
    cout << "... The Cl and Cl_inv matrices will be calculated for l = " << l << endl;

    Cl = compute_Cl(l, *Pk_index, *Tb_index, *q_index, krange, bessels, bessels_lminus1);
    Cl_inv = Cl.i();

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
            double cond_num = 0;
            cout << "Computation of Fl starts for l = " << l << endl;
            double fl = this->compute_Fl(l, param_key1, param_key2, ksteps_Cl, &cond_num, &Pk_index,\
                    &Tb_index, &q_index);
            cout << "fl with l = " << l << " is: " << fl << endl;
            Fl_file << l << " " << fl << " " << cond_num << endl;
            sum += (2*l + 1) * fl;
        }
    }
    return sum;
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
            cout << "Computation of Fl starts for # of kmodes = " << ksteps_Cl + 1 << endl;
            double fl = this->compute_Fl(l, param_key1, param_key2, ksteps_Cl, &cond_num, &Pk_index,\
                    &Tb_index, &q_index);
            cout << "fl with number of k modes sampled = " << ksteps_Cl + 1 << " is: " << fl << endl;
            Fl_file << ksteps_Cl + 1 << " " << fl << " " << cond_num << endl;            
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
