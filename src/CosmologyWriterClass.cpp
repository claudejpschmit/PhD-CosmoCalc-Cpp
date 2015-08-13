#include "CosmologyWriterClass.hpp"
#include <fstream>

CosmoWrite::CosmoWrite(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
    :
        CosmoCalc(params,Pk_index,Tb_index,q_index)
{
    cout << "... Beginning to build CosmoWriter ..." << endl;

    //this->create_bessel_interpolant_ALGLIB(this->fiducial_params["l_min"], this->fiducial_params["l_max"]);

    cout << "... CosmoWriter built ..." << endl;
}

CosmoWrite::~CosmoWrite()
{}

void CosmoWrite::calculate_Ml(int l, double k_fixed, double k2_low, double k2_high,\
        int step, double stepsize)
{
    if (stepsize == 0)
        stepsize = (k2_high - k2_low) / (double)step;
    else
        step = (int)((k2_high - k2_low)/stepsize);
    string filename;
    stringstream ss;
    ss << "output/M_" << l << "(" << k_fixed << ",[" << k2_low << "," << k2_high <<\
        "])." << step << "steps.dat";
    filename = ss.str();

    double k2, res;
    ofstream file;
    file.open(filename);
    for (int n = 0; n < step; ++n) {
        k2 = k2_low + n * stepsize;
        res = this->M(l, k_fixed, k2,0,0,0);
        file << k2 << " " << res << endl;
    }

    file.close();
}

void CosmoWrite::calculate_Nl(int l, double k_fixed, double k2_low, double k2_high,\
        int step, double stepsize)
{
    if (stepsize == 0)
        stepsize = (k2_high - k2_low) / (double)step;
    else
        step = (int)((k2_high - k2_low)/stepsize);
    string filename;
    stringstream ss;
    ss << "output/N_" << l << "(" << k_fixed << ",[" << k2_low << "," << k2_high <<\
        "])." << step << "steps.dat";
    filename = ss.str();

    double k2, res;
    ofstream file;
    file.open(filename);
    for (int n = 0; n < step; ++n) {
        k2 = k2_low + n * stepsize;
        res = this->N_bar(l, k_fixed, k2,0,0,0);
        file << k2 << " " << res << endl;
    }

    file.close();
}
void CosmoWrite::calculate_comoving(double zmax)
{
    ofstream file;
    string filename = "tabulated_values/comoving_r.dat";
    file.open(filename);
    file << "# This file contains data for comoving distances vs redshift." << endl;
    file << "# Column 2: Line of sight comoving distance D_C" << endl;

    double y, x;
    int steps = 10 * zmax;
    for (int i = 0; i < steps; i++) {
        x = (double)i/10.0;
        y = this->D_C(x); 

        file << x << " " << y << endl;
    }
    file.close();
}
void CosmoWrite::calculate_inverse_comoving(double rmax)
{
    ofstream file;
    string filename = "output/inverse_comoving_r.dat";
    file.open(filename);
    file << "# This file contains data for the inverse comoving distances, ie z vs r." << endl;
    file << "# Column 2: redshift" << endl;

    double y, x;
    int steps = 10 * rmax;
    for (int i = 0; i < steps; i++) {
        x = (double)i/10.0;
        y = this->r_inverse(x); 

        file << x << " " << y << endl;
    }
    file.close();
}

void CosmoWrite::calculate_distances(double zmax)
{
    ofstream file;
    string filename = "output/distances.dat";
    file.open(filename);
    file << "# This file contains data for 4 distances vs redshift. " << endl;
    file << "# Column 1: Redshift z" << endl;
    file << "# Column 2: Angular diameter distance D_A" << endl;
    file << "# Column 3: Luminosity distance D_L" << endl;
    file << "# Column 4: Line of sight comovong distance D_C" << endl;
    file << "# Column 5: Distance based on light travel time D_ltt" << endl;

    double norm = this->H_0/this->c * 1000.0;
    double y1, y2, y3, y4, x;
    int steps = 100 * zmax;
    for (int i = 0; i < steps; i++) {
        x = (double)i/100.0;
        y1 = this->D_A(x) * norm;
        y2 = this->D_L(x) * norm;
        y3 = this->D_now(x); //* norm;
        y4 = this->D_ltt(x) * norm;

        file << x << " " << y1 << " " << y2 << " " << y3 << " " << y4 << endl;
    }
    file.close();
}
void CosmoWrite::calculate_densities_rho(double zmax)
{
    ofstream file;
    string filename = "output/densities_rho.dat";
    file.open(filename);
    file << "# This file contains data for the densities rho vs redshift. " <<\
        endl;
    file << "# Column 1: Redshift z" << endl;
    file << "# Column 2: rho_M" << endl;
    file << "# Column 3: rho_R" << endl;
    file << "# Column 4: rho_V" << endl;

    double y1, y2, y3, x;
    int steps = 100 * zmax;
    for (int i = 0; i < steps; i++) {
        x = (double)i/100.0;
        y1 = this->rho_M(x);
        y2 = this->rho_R(x);
        y3 = this->rho_V(x);

        file << x << " " << y1 << " " << y2 << " " << y3 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_densities_Omega(double zmax)
{
    ofstream file;
    string filename = "output/densities_Omega.dat";
    file.open(filename);
    file << "# This file contains data for the densities Omega vs redshift. " <<\
        endl;
    file << "# Column 1: Redshift z" << endl;
    file << "# Column 2: Omega_M" << endl;
    file << "# Column 3: Omega_R" << endl;
    file << "# Column 4: Omega_V" << endl;
    file << "# Column 5: Sum" << endl;

    double y1, y2, y3, y4, x;
    int steps = 100 * zmax;
    for (int i = 0; i < steps; i++) {
        x = (double)i/100.0;
        y1 = this->Omega_M(x);
        y2 = this->Omega_R(x);
        y3 = this->Omega_V(x);
        y4 = y1 + y2 + y3;

        file << x << " " << y1 << " " << y2 << " " << y3 << " " << y4 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_H(double zmax)
{
    ofstream file;
    string filename = "output/hvsZ.dat";
    file.open(filename);
    file << "# This file contains data for the Hubble parameter h vs redshift. "\
        << endl;
    file << "# Column 1: Redshift z" << endl;
    file << "# Column 2: h" << endl;

    double y1, x;
    int steps = 100 * zmax;
    for (int i = 0; i < steps; i++) {
        x = (double)i/100.0;
        y1 = this->H(x)/100.0;

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_P(double kmin, double kmax, int step,\
        string unit_k, string unit_P)
{
    double stepsize = (kmax - kmin) / (double)step;

    ofstream file;
    string filename = "output/P_analytic.dat";
    file.open(filename);
    file << "# This file contains data for the analytic power spectrum\
        vs k at redshift 0." << endl;
    file << "# Column 1: scales k in units: " << unit_k << endl;
    file << "# Column 2: power spectrum in units: " << unit_P << endl;

    double y1, x;
    for (int i = 0; i < step; i++) {
        x = kmin + (double)i * stepsize ;
        y1 = this->P_delta(x, unit_k, unit_P);

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_P_CLASS(double kmin, double kmax, double z, int step)
{
    double stepsize = (kmax - kmin) / (double)step;

    ofstream file;
    string filename = "output/P_CLASS.dat";
    file.open(filename);
    file << "# This file contains data for the CLASS power spectrum\
        vs k at redshift z." << endl;
    file << "# Column 1: scales k in units: h/Mpc " << endl;
    file << "# Column 2: power spectrum in units: h/Mpc^3 " << endl;

    double y1, x;
    for (int i = 0; i < step; i++) {
        x = kmin + (double)i * stepsize ;
        y1 = this->Pk_interp(x, z,0);

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_dTb(double zmin, double zmax, int step)
{
    double stepsize = (zmax - zmin) / (double)step;

    ofstream file;
    string filename = "output/dTb.dat";
    file.open(filename);
    file << "# This file contains data for delta T_b vs redshift." << endl;
    file << "# Column 1: redshift z " << endl;
    file << "# Column 2: delta T_b" << endl; 
    file << "# Column 3: delta T_b from G21 code" << endl;
    file << "# Column 4: delta T_b from analytic interpolation" << endl;

    int index_g21 = 0;
    int index_analytic = 0;
    update_G21(fiducial_params, &index_g21);
    update_Tb_analytic(fiducial_params, &index_analytic);

    double y1, y2, y3, z;
    for (int i = 0; i < step; i++) {
        z = zmin + (double)i * stepsize ;
        y1 = delta_Tb_bar(z);
        y2 = Tb_interp(z, index_g21);
        y3 = Tb_analytic_interp(z, index_analytic);
        file << z << " " << y1 << " " << y2 << " " << y3 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_xHI(double zmin, double zmax, int step)
{
    double stepsize = (zmax - zmin) / (double)step;

    ofstream file;
    string filename = "output/xHI.dat";
    file.open(filename);
    file << "# This file contains data for x_HI vs redshift." << endl;
    file << "# Column 1: redshift z " << endl;
    file << "# Column 2: delta x_HI" << endl;

    double y1, x;
    for (int i = 0; i < step; i++) {
        x = zmin + (double)i * stepsize ;
        y1 = this->x_HI(x);

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_Ts(double zmin, double zmax, int step)
{
    double stepsize = (zmax - zmin) / (double)step;

    ofstream file;
    string filename = "output/Ts.dat";
    file.open(filename);
    file << "# This file contains data for T_S vs redshift." << endl;
    file << "# Column 1: redshift z " << endl;
    file << "# Column 2: T_S" << endl;

    double y1, x;
    for (int i = 0; i < step; i++) {
        x = zmin + (double)i * stepsize ;
        y1 = this->T_S(x);

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_Tk(double zmin, double zmax, int step)
{
    double stepsize = (zmax - zmin) / (double)step;

    ofstream file;
    string filename = "output/Tk.dat";
    file.open(filename);
    file << "# This file contains data T_K vs redshift." << endl;
    file << "# Column 1: redshift z " << endl;
    file << "# Column 2: T_K" << endl;

    double y1, x;
    for (int i = 0; i < step; i++) {
        x = zmin + (double)i * stepsize ;
        y1 = this->T_K(x);

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_P_compare(double k_low, double k_high, int kstep, double z_low, double z_high, int zstep)
{
    double stepsize_k = (k_high - k_low) / (double)kstep;
    double stepsize_z = (z_high - z_low) / (double)zstep;
    double k, z, res;
    ofstream file;
    string filename = "output/Pk_Cpp.dat";
    file.open(filename);
    file << "# The parameters in this calculation are as follows:" << endl;
    for (auto& x: this->current_params)
        file << "# " << x.first << " = " << x.second << endl;
    file << "# ";
    for (int j = 0; j < zstep; j++) {
        z = z_low + (double)j * stepsize_z;
        file << z << " ";
    }
    file << endl;

    for (int i = 0; i < kstep; i++) {
        k = k_low + (double)i * stepsize_k;
        file << k << " ";
        for (int j = 0; j < zstep; j++) {
            z = z_low + (double)j * stepsize_z;
            res = this->Pk_interp(k,z,0);

            file << res << " ";
        }
        file << endl;
    }
    file.close();
}

void CosmoWrite::calculate_bessels(int l) 
{
    ofstream file;
    string filename = "output/bessel.dat";
    file.open(filename);
    file << "# This file contains the interpolated bessel function j_" << l << " from 1 to 100000." << endl;
    file << "# Column 1: x " << endl;
    file << "# Column 2: j_"<< l << "(x)" << endl;

    double y1, x;
    for (int i = 1; i < 1000000; i++) {
        x = i/10.0;
        y1 = this->bessel_j_interp(l, x);

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_bessels_basic(int l) 
{
    ofstream file;
    string filename = "output/bessel_basic.dat";
    file.open(filename);
    file << "# This file contains the interpolated bessel function j_" << l << " from 1 to 100000." << endl;
    file << "# Column 1: x " << endl;
    file << "# Column 2: j_"<< l << "(x)" << endl;

    double y1, x;
    for (int i = 1; i < 1000000; i++) {
        x = i/10.0;
        y1 = this->bessel_j_interp_basic(l, x);

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_bessels_cubic(int l) 
{
    ofstream file;
    string filename = "output/bessel_cubic.dat";
    file.open(filename);
    file << "# This file contains the interpolated bessel function j_" << l << " from 1 to 100000." << endl;
    file << "# Column 1: x " << endl;
    file << "# Column 2: j_"<< l << "(x)" << endl;

    double y1, x;
    for (int i = 1; i < 100000; i++) {
        x = i/10.0;
        y1 = this->bessel_j_interp_cubic(l, x);

        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_bessels_exact(int l) 
{
    ofstream file;
    string filename = "output/bessel_exact.dat";
    file.open(filename);
    file << "# This file contains the bessel function j_" << l <<\
        " from 1 to 100000 from camb." << endl;
    file << "# Column 1: x " << endl;
    file << "# Column 2: j_"<< l << "(x) from camb." << endl;

    double x, y2;
    for (int i = 1; i < 100000; i++) {
        x = i/10.0;
        y2 = this->sph_bessel_camb(l,x);
        file << x << " " << y2 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_integrandMM(int l, double k1, double k2, int step) 
{
    ofstream file;
    string filename = "output/integrandMM_"+to_string(k1)+"_"+to_string(k2)+".dat";
    file.open(filename);
    file << "# This file contains the integrand k^2 MM  from 0.0001 to 2 for k1 = " << k1 << " and k2 = "<< k2 << endl;
    file << "# Column 1: k (the one that is integrated over)" << endl;
    file << "# Column 2: k^2 MM" << endl;

    double stepsize = (2-0.0001)/(double)step;
    double y1, x;
    for (int i = 0; i < step; i++) {
        x = 0.0001 + i*stepsize;
        y1 = this->integrandMM(l,k1,k2,x,0,0,0);
        file << x << " " << y1 << endl;
    }
    file.close();
}
void CosmoWrite::calculate_integrandMN(int l, double k1, double k2, int step)
{
    ofstream file;
    string filename = "output/integrandMN_"+to_string(k1)+"_"+to_string(k2)+".dat";
    file.open(filename);
    file << "# This file contains the integrand k^2 MN  from 0.001 to 5 for k1 = " << k1 << " and k2 = "<< k2 << endl;
    file << "# Column 1: k (the one that is integrated over)" << endl;
    file << "# Column 2: k MN" << endl;

    double stepsize = (5-0.001)/(double)step;
    double y1, x;
    for (int i = 0; i < step; i++) {
        x = 0.001 + i*stepsize;
        y1 = this->integrandMN(l,k1,k2,x,0,0,0);
        file << x << " " << y1 << endl;
    }
    file.close();
}
void CosmoWrite::calculate_integrandNN(int l, double k1, double k2, int step)
{
    ofstream file;
    string filename = "output/integrandNN_"+to_string(k1)+"_"+to_string(k2)+".dat";
    file.open(filename);
    file << "# This file contains the integrand NN  from 0.001 to 5 for k1 = " << k1 << " and k2 = "<< k2 << endl;
    file << "# Column 1: k (the one that is integrated over)" << endl;
    file << "# Column 2: NN" << endl;

    double stepsize = (5-0.001)/(double)step;
    double y1, x;
    for (int i = 0; i < step; i++) {
        x = 0.001 + i*stepsize;
        y1 = this->integrandNN(l,k1,k2,x,0,0,0);
        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_integrand_kappa_norsd(int l, double k1, double k2, int step) 
{
    ofstream file;
    stringstream filename;
    filename << "output/integrand_kappa_norsd_l" << l << "_" << k1 << "_" << k2 << ".dat";
    file.open(filename.str());
    file << "# This file contains the integrand kappa^2 MM for l = " << l <<\
        " from 0.0001 to 1 for k1 = " << k1 << " and k2 = "<< k2 << endl;
    file << "# Column 1: kappa" << endl;
    file << "# Column 2: kappa^2 MM_l=" << l << endl;

    double stepsize = (1-0.0001)/(double)step;
    double y1, x;
    for (int i = 0; i < step; i++) {
        x = 0.0001 + i*stepsize;
        y1 = this->integrand_kappa_norsd(l,k1,k2,x,0,0,0);
        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_integrand_kappa_rsd(int l, double k1, double k2, int step)
{
    ofstream file;
    stringstream filename;
    filename << "output/integrand_kappa_rsd_l" << l << "_" << k1 << "_" << k2 << ".dat";
    file.open(filename.str());
    file << "# This file contains the integrand (kappa^2 MM + bbeta (MN + NM) + bbeta^2 NN) for l = " <<\
        l << " from 0.0001 to 1 for k1 = " << k1 << " and k2 = "<< k2 << endl;
    file << "# Column 1: kappa" << endl;
    file << "# Column 2: kappa^2 MM_l=" << l << endl;

    double stepsize = (1-0.0001)/(double)step;
    double y1, x;
    for (int i = 0; i < step; i++) {
        x = 0.0001 + i*stepsize;
        y1 = this->integrand_kappa_rsd(l,k1,k2,x,0,0,0);
        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_qdot()
{
    update_q_prime();
    ofstream file;

    string filename = "output/qdot.dat";
    file.open(filename);
    file << "# This file contains the derivative of q" << endl;
    file << "# Column 1: n value (redshift measure)" << endl;
    file << "# Column 2: q_dot" << endl;

    double y1, x;
    for (unsigned int i = 0; i < q_p_Ml.size(); i++) {
        x = i;
        y1 = q_p_Ml[i];
        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_q()
{
    update_q_prime();
    ofstream file;

    string filename = "output/q.dat";
    file.open(filename);
    file << "# This file contains the derivative of q" << endl;
    file << "# Column 1: n value (redshift measure)" << endl;
    file << "# Column 2: q" << endl;

    double y1, x;
    for (unsigned int i = 0; i < q_Ml.size(); i++) {
        x = i;
        y1 = q_Ml[i];
        file << x << " " << y1 << endl;
    }
    file.close();
}
void CosmoWrite::calculate_integrandlong(int l, double k1, double k2, int step)
{
    ofstream file;
    string filename = "output/integrandlong_"+to_string(k1)+"_"+to_string(k2)+".dat";
    file.open(filename);
    file << "# This file contains the long integrand from z = 7 to 9 for k1 = " << k1 << " and k2 = "<< k2 << endl;
    file << "# Column 1: z" << endl;
    file << "# Column 2: long integrand" << endl;

    double stepsize = (9-7)/(double)step;
    double y1, x;
    for (int i = 0; i < step; i++) {
        x = 7 + i*stepsize;
        y1 = this->integrandlong(l,k1,k2,x,0,0,0);
        file << x << " " << y1 << endl;
    }
    file.close();
}
void CosmoWrite::calculate_integrandsimple(int l, double k1, double k2, int step)
{
    ofstream file;
    string filename = "output/integrandsimple_"+to_string(k1)+"_"+to_string(k2)+".dat";
    file.open(filename);
    file << "# This file contains the integrand simple  from z = 7 to 9 for k1 = " << k1 << " and k2 = "<< k2 << endl;
    file << "# Column 1: z" << endl;
    file << "# Column 2: simplified integrand" << endl;

    double stepsize = (9-7)/(double)step;
    double y1, x;
    for (int i = 0; i < step; i++) {
        x = 7 + i*stepsize;
        y1 = this->integrandsimple(l,k1,k2,x,0,0,0);
        file << x << " " << y1 << endl;
    }
    file.close();
}

void CosmoWrite::generate_movie(int l)
{
    double x, y1, y2;
    double stepsize = 2.0/100.0;
    int num = 0;
    for (int k = 0; k < 10; k++) {
        for (int i = 0; i < 10; i++) {
            ofstream file;
            file.open("movie/frame.dat");
            double k1, k2;
            k1 = 0.01 + k * 0.01;
            k2 = 0.01 + i * 0.01;
            for (int j = 0; j < 100; j++) {
                x = 7 + j*stepsize;
                y1 = this->integrandlong(l,k1,k2,x,0,0,0);
                y2 = this->integrandsimple(l,k1,k2,x,0,0,0);

                file << x << " " << y1 << " " << y2 << endl;
            }
            file.close();
            stringstream command;
            command << "python plot.py " << k1 << " " << k2 << " " << num;
            system(command.str().c_str());
            cout << "frame "<< num << " done " << endl;
            num++;
        }
    }
}

void CosmoWrite::calculate_Cl_simple(int l, double k, double k_min, double k_max, double k_stepsize)
{
    double x, y;
    ofstream file;
    string filename = "output/Cl_simple_"+to_string(l)+"_"+to_string(k)+".dat";
    file.open(filename);
    file << "# This file contains the simplified Cl(k = "<< k <<", k') with k' from "<< k_min <<\
        " to " << k_max << "." << endl;
    file << "# Column 1: k'" << endl;
    file << "# Column 2: Cl(k,k')" << endl;

    int step = (k_max - k_min)/k_stepsize;
    for (int i = 0; i < step; i++) {
        x = k_min + i*k_stepsize;
        y = this->Cl_simplified(l,k,x,0,0,0);
        file << x << " " << y << endl;
    }
    file.close();

}

void CosmoWrite::calculate_Cl_full(int l, double k, double k_min, double k_max, double k_stepsize)
{
    double x, y, y2;
    ofstream file;
    string filename = "output/Cl_full_"+to_string(l)+"_"+to_string(k)+"_new.dat";
    file.open(filename);
    file << "# This file contains the full Cl(k = "<< k <<", k') with k' from "<< k_min <<\
        " to " << k_max << "." << endl;
    file << "# Column 1: k'" << endl;
    file << "# Column 2: Cl(k,k')" << endl;

    int step = (k_max - k_min)/k_stepsize;
    for (int i = 0; i < step; i++) {
        x = k_min + i*k_stepsize;
        y2 = this->corr_Tb(l,k,x,0.0001,1.0,0,0,0);
        file << x << " " << y2 << endl;
    }
    file.close();
}

void CosmoWrite::calculate_Cl_simple_rsd(int l, double k, double k_min, double k_max, double k_stepsize)
{
    double x, y;
    ofstream file;
    string filename = "output/Cl_simple_rsd_"+to_string(l)+"_"+to_string(k)+".dat";
    file.open(filename);
    file << "# This file contains the simplified Cl(k = "<< k <<", k') with k' from "<< k_min <<\
        " to " << k_max << "." << endl;
    file << "# Column 1: k'" << endl;
    file << "# Column 2: Cl(k,k')" << endl;

    int step = (k_max - k_min)/k_stepsize;
    for (int i = 0; i < step; i++) {
        x = k_min + i*k_stepsize;
        y = this->Cl_simplified_rsd(l,k,x,0,0,0);
        file << x << " " << y << endl;
    }
    file.close();

}

void CosmoWrite::calculate_Cl_full_rsd(int l, double k, double k_min, double k_max, double k_stepsize)
{
    double x, y;
    ofstream file;
    string filename = "output/Cl_full_rsd"+to_string(l)+"_"+to_string(k)+".dat";
    file.open(filename);
    file << "# This file contains the full Cl(k = "<< k <<", k') with k' from "<< k_min <<\
        " to " << k_max << "." << endl;
    file << "# Column 1: k'" << endl;
    file << "# Column 2: Cl(k,k')" << endl;

    int step = (k_max - k_min)/k_stepsize;
    for (int i = 0; i < step; i++) {
        x = k_min + i*k_stepsize;
        y = this->corr_Tb_rsd(l,k,x,0.0001,1.0,0,0,0);
        file << x << " " << y << endl;
    }
    file.close();
}

void CosmoWrite::generate_movie_Cl(int l_min, int l_max, double k, double k_min,\
        double k_max, double k_stepsize)
{
    double x, y;
    int num = 0;
    int step = (k_max - k_min)/k_stepsize;
    int l_steps = l_max - l_min;
    for (int i = 0; i < l_steps; i++) {
        ofstream file;
        file.open("movie_Cl/frame.dat");
        int l = l_min + i;
        for (int j = 0; j < step; j++) {
            x = k_min + j*k_stepsize;
            y = this->Cl_simplified(l,k,x,0,0,0);

            file << x << " " << y << endl;
        }
        file.close();
        stringstream command;
        command << "python plot_Cl.py " << l << " " << num;
        system(command.str().c_str());
        cout << "frame "<< num << " done " << endl;
        num++;
    }

}
void CosmoWrite::calculate_Ml(int lmin, int lmax, double k, double kappa)
{
    double x, y;
    ofstream file;
    string filename = "output/Ml_"+to_string(k)+"_"+to_string(kappa)+".dat";
    file.open(filename);
    file << "# This file contains Ml(k = "<< k <<", kappa = " << kappa << ") with l from "<< lmin <<\
        " to " << lmax << "." << endl;
    file << "# Column 1: l'" << endl;
    file << "# Column 2: Ml(k,kappa)" << endl;

    int step = (lmax - lmin);
    for (int i = 0; i < step; i++) {
        x = lmin + i;
        y = this->M(x, k, kappa,0,0,0);
        file << x << " " << y << endl;
    }
    file.close();

}
void CosmoWrite::calculate_Nl(int lmin, int lmax, double k, double kappa)
{
    double x, y;
    ofstream file;
    string filename = "output/Nl_"+to_string(k)+"_"+to_string(kappa)+".dat";
    file.open(filename);
    file << "# This file contains Nl(k = "<< k <<", kappa = " << kappa << ") with l from "<< lmin <<\
        " to " << lmax << "." << endl;
    file << "# Column 1: l'" << endl;
    file << "# Column 2: Nl(k,kappa)" << endl;

    int step = (lmax - lmin);
    for (int i = 0; i < step; i++) {
        x = lmin + i;
        y = this->N_bar(x, k, kappa,0,0,0);
        file << x << " " << y << endl;
    }
    file.close();

}
