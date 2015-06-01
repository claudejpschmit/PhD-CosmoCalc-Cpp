#include "CosmologyWriterClass.hpp"
#include <fstream>

CosmoWrite::CosmoWrite(map<string, double> params)
    :
        CosmoCalc(params)
{

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
        res = this->M(l, k_fixed, k2);
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
        res = this->N_bar(l, k_fixed, k2);
        file << k2 << " " << res << endl;
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
        y3 = this->D_now(x) * norm;
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
        y1 = this->Pk_interp(x, z);

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

    double y1, x;
    for (int i = 0; i < step; i++) {
        x = zmin + (double)i * stepsize ;
        y1 = this->delta_Tb_bar(x);

        file << x << " " << y1 << endl;
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
            res = this->Pk_interp(k,z);

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

void CosmoWrite::calculate_bessels_exact(int l) 
{
    ofstream file;
    string filename = "output/bessel_boost.dat";
    file.open(filename);
    file << "# This file contains the bessel function j_" << l <<\
        " from 1 to 100000 from boost and camb." << endl;
    file << "# Column 1: x " << endl;
    file << "# Column 2: j_"<< l << "(x) from boost." << endl;
    file << "# Column 3: j_"<< l << "(x) from camb." << endl;

    double y1, x, y2;
    for (int i = 1; i < 1000000; i++) {
        x = i/10.0;
        y1 = this->sph_bessel(l, x);
        y2 = this->sph_bessel_camb(l,x);
        file << x << " " << y1 << " " << y2 << endl;
    }
    file.close();
}
