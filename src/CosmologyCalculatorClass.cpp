#include "CosmologyCalculatorClass.hpp"
#include <string>
#include <iostream>
#include "Integrator.hpp"
CosmoCalc::CosmoCalc(map<string, double> params)
    :
        CosmoBasis(params)
{
    this->k_steps = this->fiducial_params["k_steps"];
    this->prefactor_Ml = 2*this->b_bias*this->c/pi;
    this->zmin_Ml = this->fiducial_params["zmin"];
    this->zmax_Ml = this->fiducial_params["zmax"];
    this->zsteps_Ml = this->fiducial_params["zsteps"];
    this->stepsize_Ml = (this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    this->Pk_steps = this->fiducial_params["Pk_steps"];
    this->k_steps = this->fiducial_params["k_steps"];
    //TODO: Now we need to predefine the lists q_Ml & r_Ml & H_f....
}

CosmoCalc::~CosmoCalc()
{
}
void CosmoCalc::show_cosmo_calcs()
{
    cout << hubble_time() << endl;
    cout << hubble_dist() << endl;
    cout << comoving_radial_dist(10) << endl;
    cout << "O_m = " << O_M << ", O_V = " << O_V << "." << endl;
    cout << "Age in Gigayears "<< age_of_universe(0) * pow(10,10) * 3.08568 / (365.25 * 24 * 3600) << endl;

}

double CosmoCalc::hubble_time()
{
    return this->t_H;
}
double CosmoCalc::hubble_dist()
{
    return this->D_H;
}
double CosmoCalc::comoving_radial_dist(double z)
{
    return this->hubble_dist() * this->Z(z);
}
double CosmoCalc::D_C(double z)
{
    return this->comoving_radial_dist(z);
}

double CosmoCalc::D_now(double z)
{
    return this->comoving_radial_dist(z);
}

double CosmoCalc::comoving_dist_transverse(double z)
{
    const double Ok = 1 - this->O_M - this->O_V;
    const double arg = sqrt(abs(1-Ok))*this->D_C(z)/this->D_H;

    return this->D_H / sqrt(abs(1-Ok)) * this->S_k(arg);
}

double CosmoCalc::D_M(double z)
{
    return this->comoving_dist_transverse(z);
}

double CosmoCalc::angular_diam_dist(double z, double z2 = -1)
{
    const double Ok = 1 - this->O_M - this->O_V;
    const double root = sqrt(abs(1-this->O_tot));
    double result;
    if (z2 < 0.0)
    {
        result = this->D_H * this->S_k(root * this-> Z(z)) / ((1+z) * root);
    }
    else if (this->O_V + this->O_M <= 1.0)
    {
        const double dm = this->D_M(z);
        const double dm2 = this->D_M(z2);
        result = 1.0 / (1.0 + z2) *\
                 ( dm2 * sqrt(1 + Ok * pow(dm,2) / pow(this->D_H, 2)) -\
                   dm * sqrt(1 + Ok * pow(dm2,2) / pow(this->D_H,2)) );
    }
    else
    {
        cout << "Error: D_A12 formula invalid for O_tot > 1.0" << endl;
        result = 1.0;
    }

    return result;
}

double CosmoCalc::D_A(double z, double z2 = -1)
{
    return this->angular_diam_dist(z,z2);
}

double CosmoCalc::luminosity_dist(double z)
{
    return this->D_A(z) * pow(1+z,2);
}

double CosmoCalc::D_L(double z)
{
    return this->luminosity_dist(z);
}

double CosmoCalc::comoving_volume(double z)
{
    const double vol = 4*pi*this->D_H;
    auto integrand = [&](double x){return pow(1+x,2) * pow(this->D_A(x),2)/this->E(x);};
    const double integral = integrate(integrand, 0.0, z, 1000, simpson()); 
    return vol * integral;
}

double CosmoCalc::V_C(double z)
{
    return comoving_volume(z);
}

// infinite integrals = not fun!
// This seems to be large enough for the upper bound.

double CosmoCalc::age_of_universe(double z)
{
    auto integrand = [&](double x){return 1.0/((1+x)*this->E(x));};
    return this->t_H * integrate(integrand, z, pow(10,4), 100000, simpson());
}

double CosmoCalc::light_travel_time(double z)
{
    return this->age_of_universe(0) - this->age_of_universe(z);
}

double CosmoCalc::D_ltt(double z)
{
    return this->light_travel_time(z) * this->c / 1000.0;
}

double CosmoCalc::H(double z)
{
    return this->H_0 * this->E(z);
}

double CosmoCalc::H_SI(double z)
{
    return this->H(z) * 1000.0 / (3.08567758 * pow(10,16));
}

double CosmoCalc::rho_crit(double z)
{
    return 3.0 * pow(this->H_SI(z),2) / (8.0 * pi * this->G);
}

double CosmoCalc::rho_M(double z)
{
    const double rho_M_0 = this->O_M * this->rho_crit(0);
    return rho_M_0 * pow((1+z),3);
}

double CosmoCalc::rho_R(double z)
{
    const double rho_R_0 = this->O_R * this->rho_crit(0);
    return rho_R_0 * pow((1+z),4);
}

double CosmoCalc::rho_V(double z)
{
    (void) z;
    return this->O_V * this->rho_crit(0);
}

double CosmoCalc::Omega_M(double z)
{
    return this->O_M * pow(1+z,3) / pow(this->E(z),2);
}

double CosmoCalc::Omega_R(double z)
{
    return this->O_R * pow(1+z,4) / pow(this->E(z),2);
}

double CosmoCalc::Omega_V(double z)
{
    return this->O_V / pow(this->E(z),2);
}

double CosmoCalc::num_baryons()
{
    const double num = 4.0/3.0 * pi * pow(this->c/this->H_SI(0),3) * this->rho_crit(0);
    return num / this->m_b;
}

double CosmoCalc::T(double z)
{
    return this->T_CMB * (1+z);
}

double CosmoCalc::n_H_tot(double z)
{
    return 1.6 * pow(1+z,3);
}

double CosmoCalc::n_b(double z)
{
    return this->O_b * this->rho_crit(z) / this->m_b;
}

double CosmoCalc::n_H(double z)
{
    const double coeff = 1.0;
    return coeff * this->n_b(z);
}

double CosmoCalc::n_p(double z)
{
    return this->x_HI(z);
}

double CosmoCalc::n_e(double z)
{
    return this->n_p(z);
}

void CosmoCalc::update_q()
{
}

void CosmoCalc::Pk_update_interpolator(map<string, double> params)
{
}

double CosmoCalc::Pkz_calc(double k, double z)
{
    return this->P_growth(z) * this->P_delta(k);
}

double CosmoCalc::P_growth(double z)
{
    const double res = this->D1(z) / this->D1(0);
    return pow(res,2);
}

double CosmoCalc::D1(double z)
{
    const double prefactor = 5 * this->O_M / 2 * this->E(z);
    auto integrand = [&](double x){return (1+x)/pow(this->E(x),3);};
    return prefactor * integrate(integrand, z, pow(10, 5), 1000000, simpson());
}

double CosmoCalc::P_delta(double k, string units_k, string units_P)
{
    double keq, k_factor, delta_H;
    if (units_P == "default")
    {
        if (units_k == "default")
        {
            keq = 0.073 * this->O_M * this->h;
            k_factor = k;
        }
        else if (units_k == "Mpc-1" or units_k == "mpc-1")
        {
            keq = 0.073 * this->O_M * pow(this->h,2);
            k_factor = k / this->h;

        }
    }
    else if (units_P == "Mpc3" or units_P == "mpc3")
    {
        if (units_k == "default")
        {
            keq = 0.073 * this->O_M * this->h;
            k_factor = k / pow(this->h,3);

        }
        else if (units_k == "Mpc-1" or units_k == "mpc-1")
        {
            keq = 0.073 * this->O_M * pow(this->h,2);
            k_factor = k / pow(this->h,4);
        }
    }
    if (this->O_M == 1)
    {
        delta_H = 1.9 * pow(10,-5);
    }
    else
    {
        delta_H = 4.5 * pow(10,-5);
    }
    const double A = 2*pow(pi,2) * pow(this->c/1000.0,4) * pow(delta_H,2)/pow(10,8);
    const double x = k/keq;
    const double transfer_function_sq = pow(this->transfer(x),2);
    return A * k_factor * transfer_function_sq;
}

double CosmoCalc::transfer(double x)
{
    const double res = log(1 + 0.171 * x) / (0.171 * x);
    const double bracket = 1 + 0.284 * x + pow(1.18 * x, 2) + pow(0.399 * x, 3) + pow(0.49 * x, 4);
    return res * pow(bracket, -0.25);
}


double CosmoCalc::corr_Tb(int l, double k1, double k2, double k_low, double k_high)
{
    if (k1 == k2)
    {
        auto integrand = [&](double k){return pow(k,2) * pow(this->M(l,k1,k),2);};
        return integrate(integrand, k_low, k_high, this->k_steps, simpson());     
    }
    else
    {
        auto integrand = [&](double k){return pow(k,2) * this->M(l,k1,k) * this->M(l,k2,k);};
        return integrate(integrand, k_low, k_high, this->k_steps, simpson());
    }
}

double CosmoCalc::corr_Tb_rsd(int l, double k1, double k2, double k_low, double k_high)
{
    auto integrand = [&](double k)
    {
        double m1,n1,m2,n2;
        m1 = this->M(l,k1,k);
        n1 = this->N_bar(l,k1,k);
        if (k1 == k2)
        {
            m2 = m1;
            n2 = n1;
        }
        else
        {
            m2 = this->M(l,k2,k);
            n2 = this->N_bar(l,k2,k);
        }
        const double bb = this->b_bias * this->beta;
        const double bb2 = pow(bb,2);
        return pow(k,2) * m1 * m2 + bb * k * (m1*n2 + n1*m2) + bb2 * n1 * n2;
    };

    return integrate(integrand, k_low, k_high, this->k_steps, simpson());

}

double CosmoCalc::M(int l, double k1, double k2)
{
    auto integrand = [&](double z)
    {
        const double n_old = (z - this->zmin_Ml)/this->stepsize_Ml;
        int n;
        int n_old_int = (int)n_old;
        if (abs(n_old - (double)n_old_int) > 0.5)
            n = n_old_int + 1;
        else
            n = n_old_int;
        double r,q;
        r = 1;
        q = 1;
        
        return 0;

        //return pow(r,2) * this->delta_Tb_bar(z) * this->sph_bessel(l,k1*r) *\
                this->sph_bessel(l,k2*q) * sqrt(this->Pk_interp(k2*this->h,z)/pow(this->h,3)/\
                (this->H_f[n]*1000.0);
    };

    return 0;
}

double CosmoCalc::delta_Tb_bar(double z)
{
    return 0;
}

double CosmoCalc::T_S(double z)
{
    return 0;
}

double CosmoCalc::x_HI(double z)
{
    return 0;
}

double CosmoCalc::T_K(double z)
{
    return 0;
}

double CosmoCalc::I(int l1, int l2, double k1, double k2, double z, double r)
{
    return 0;
}

double CosmoCalc::N_bar(int l, double k1, double k2)
{
    return 0;
}

