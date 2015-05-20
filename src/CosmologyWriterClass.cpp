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
void CosmoWrite::calculate_densities_rho(){}

void CosmoWrite::calculate_densities_Omega(){}

void CosmoWrite::calculate_H(){}

void CosmoWrite::calculate_P(){}

void CosmoWrite::calculate_P_CLASS(){}

void CosmoWrite::calculate_dTb(){}

void CosmoWrite::calculate_xHI(){}

void CosmoWrite::calculate_Ts(){}

void CosmoWrite::calculate_Tk(){}


