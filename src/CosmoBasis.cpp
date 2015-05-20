#include "CosmoBasis.hpp"
#include <iostream>
#include <string>
#include "Integrator.hpp"
#include <boost/math/special_functions/bessel.hpp>
using namespace std;

CosmoBasis::CosmoBasis(map<string,double> params)
{
    // Initializing parameters
    fiducial_params = params;
    this->check_params();
    current_params = fiducial_params;
    this->generate_params(fiducial_params);
    
    T_star = h_planck * c / (k_b *0.21);
    b_bias = pow(O_M,0.6) / beta;
    k_eq = 0.073 * O_M * pow(h,2);

}

CosmoBasis::~CosmoBasis()
{
}

void CosmoBasis::show_params()
{
    cout << c << endl;
    cout << k_b << endl;
    cout << m_b << endl;
    cout << m_e << endl;
    cout << h_planck << endl;
    cout << e << endl;
    cout << G << endl;
    cout << T_star << endl;
    
    cout << "ombh2 = "<< fiducial_params["ombh2"] << endl;
    cout << "hubble = "<< fiducial_params["hubble"] << endl;
    cout << T_CMB << endl;
}

void CosmoBasis::check_params()
{
    // insert automatically checks if the key is present and only adds if it isn't.
    // There isn't really a convenient way in C++ to store pairs with different
    // types in the same map, so I keep all as double so one might need to cast some
    // to integers when used later on.
    this->fiducial_params.insert(pair<string,double>("ombh2",0.0226));
    this->fiducial_params.insert(pair<string,double>("omch2",0.112));
    this->fiducial_params.insert(pair<string,double>("omnuh2",0.00064));
    this->fiducial_params.insert(pair<string,double>("omk",0.0));
    this->fiducial_params.insert(pair<string,double>("hubble",70.0));
    this->fiducial_params.insert(pair<string,double>("T_CMB",2.7255));
    this->fiducial_params.insert(pair<string,double>("zmin",7.0));
    this->fiducial_params.insert(pair<string,double>("zmax",9.0));
    this->fiducial_params.insert(pair<string,double>("zsteps",100));
    this->fiducial_params.insert(pair<string,double>("Pk_steps",3));
    this->fiducial_params.insert(pair<string,double>("k_steps",40000));
    
    this->fiducial_params.insert(pair<string,double>("100*theta_s",1.04));
    this->fiducial_params.insert(pair<string,double>("A_s",2.42e-9));
    this->fiducial_params.insert(pair<string,double>("n_s",0.96));
    this->fiducial_params.insert(pair<string,double>("tau_reio",0.09));
    this->fiducial_params.insert(pair<string,double>("k_pivot",0.05));
    this->fiducial_params.insert(pair<string,double>("YHe",0.25));
    this->fiducial_params.insert(pair<string,double>("z_pk",7.0));
}

void CosmoBasis::generate_params(map<string,double> params)
{
    T_CMB = params["T_CMB"];
    T_gamma = T_CMB;
    H_0 = params["hubble"];
    h = H_0 / 100.0;
    O_b = params["ombh2"] / pow(h,2);
    O_cdm = params["omch2"] / pow(h,2);
    O_nu = params["omnuh2"] / pow(h,2);
    O_gamma = pow(pi,2) * pow(T_CMB/11605.0,4) / (15.0*8.098*pow(10,-11)*pow(h,2));
    O_nu_rel = O_gamma * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
    O_R = O_gamma + O_nu_rel;
    O_k = params["omk"];
    O_M = O_b + O_cdm + O_nu;
    O_tot = 1.0 - O_k;
    O_V = O_tot - O_M - O_R;
    D_H = c / (1000.0 * H_0);
    t_H = 1.0 / H_0;
         
}

double CosmoBasis::sph_bessel(unsigned int l, double x)
{
    return boost::math::sph_bessel(l,x);
}


double CosmoBasis::E(double z)
{
    return sqrt(O_V + O_R * pow(1+z,4) + O_M * pow(1+z,3) + O_k * pow(1+z,2));
}

double CosmoBasis::Z(double z)
{
    auto integrand = [&](double x){return 1/E(x);};
    return integrate(integrand, 0.0, z, 1000, simpson());
}

double CosmoBasis::S_k(double x)
{
    if (O_tot < 1.0)
        return sinh(x);
    else if (O_tot == 1.0)
        return x;
    else
        return sin(x);
}

double CosmoBasis::mpc_to_m(double x)
{
    const double conv_factor = 3.0856776 * pow(10,22);
    return x*conv_factor;
}

double CosmoBasis::m_to_mpc(double x)
{
    const double conv_factor = 3.0856776 * pow(10,22);
    return x*conv_factor;
}
void CosmoBasis::params_to_planck15(map<string, double> params)
{
    (void) params;
}
