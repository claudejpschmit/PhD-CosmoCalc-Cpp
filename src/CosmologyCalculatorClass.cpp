#include "CosmologyCalculatorClass.hpp"
#include <string>

CosmoCalc::CosmoCalc(map<string, double> params)
    :
        CosmoBasis(params)
{
}

CosmoCalc::~CosmoCalc()
{
}
void CosmoCalc::show_cosmo_calcs()
{
}

double CosmoCalc::hubble_time()
{
    return 0;
}
double CosmoCalc::hubble_dist()
{
    return 0;
}
double CosmoCalc::comoving_radial_dist(double z)
{
    return 0;
}
double CosmoCalc::D_C(double z)
{
    return 0;
}

double CosmoCalc::D_now(double z)
{
    return 0;
}

double CosmoCalc::comoving_dist_transverse(double z)
{
    return 0;
}

double CosmoCalc::D_M(double z)
{
    return 0;
}

double CosmoCalc::angular_diam_dist(double z, double z2 = -1)
{
    return 0;
}

double CosmoCalc::D_A(double z, double z2 = -1)
{
    return 0;
}

double CosmoCalc::luminosity_dist(double z)
{
    return 0;
}

double CosmoCalc::D_L(double z)
{
    return 0;
}

double CosmoCalc::comoving_volume(double z)
{
    return 0;
}

double CosmoCalc::V_C(double z)
{
    return 0;
}

double CosmoCalc::age_of_universe(double z)
{
    return 0;
}

double CosmoCalc::light_travel_time(double z)
{
    return 0;
}

double CosmoCalc::D_ltt(double z)
{
    return 0;
}

double CosmoCalc::H(double z)
{
    return 0;
}

double CosmoCalc::H_SI(double z)
{
    return 0;
}

double CosmoCalc::rho_crit(double z)
{
    return 0;
}

double CosmoCalc::rho_M(double z)
{
    return 0;
}

double CosmoCalc::rho_R(double z)
{
    return 0;
}

double CosmoCalc::rho_V(double z)
{
    return 0;
}

double CosmoCalc::Omega_M(double z)
{
    return 0;
}

double CosmoCalc::Omega_R(double z)
{
    return 0;
}

double CosmoCalc::Omega_V(double z)
{
    return 0;
}

double CosmoCalc::num_baryons()
{
    return 0;
}

double CosmoCalc::T(double z)
{
    return 0;
}

double CosmoCalc::n_H_tot(double z)
{
    return 0;
}

double CosmoCalc::n_b(double z)
{
    return 0;
}

double CosmoCalc::n_H(double z)
{
    return 0;
}

double CosmoCalc::n_p(double z)
{
    return 0;
}

double CosmoCalc::n_e(double z)
{
    return 0;
}


void CosmoCalc::Pk_update_interpolator(map<string, double> params, double zi, double zf, int nsteps)
{
}

double CosmoCalc::Pkz_calc(double k, double z)
{
    return 0;
}

double CosmoCalc::P_growth(double z)
{
    return 0;
}

double CosmoCalc::D1(double z)
{
    return 0;
}

double CosmoCalc::P_delta(double k, string units_k = "default", string units_P = "default")
{
    return 0;
}

double CosmoCalc::transfer(double x)
{
    return 0;
}


double CosmoCalc::corr_Tb(int l, double k1, double k2, double k_low, double k_high)
{
    return 0;
}

double CosmoCalc::corr_Tb_rsd(int l, double k1, double k2, double k_low, double k_high)
{
    return 0;
}

double CosmoCalc::M(int l, double k1, double k2)
{
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

double CosmoCalc::N_bar(int l, double k1, double k2, double z_low, double z_high)
{
    return 0;
}

