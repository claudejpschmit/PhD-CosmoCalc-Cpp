#include "CosmologyCalculatorClass.hpp"
#include <string>
#include <iostream>
#include "Integrator.hpp"
#include <fstream>

#define ARES_MACRO false

CosmoCalc::CosmoCalc(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
    :
        CosmoBasis(params)
{
    cout << "... Beginning to build CosmoCalc ..." << endl;
    this->prefactor_Ml = 2*this->b_bias*this->c/pi;
    this->zmin_Ml = this->fiducial_params["zmin"];
    this->zmax_Ml = this->fiducial_params["zmax"];
    this->zsteps_Ml = this->fiducial_params["zsteps"];
    this->stepsize_Ml = abs(this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    this->Pk_steps = this->fiducial_params["Pk_steps"];
    this->k_stepsize = this->fiducial_params["k_stepsize"];
    //this->zmax_interp = this->fiducial_params["zmax_interp"];
    //generate object that is the CAMB interface.
    CAMB = new CAMB_CALLER;

    cout << "... precalculating Ml dependencies ..." << endl;
    //this->update_q_full();
    //this->update_q_prime_full();
    //this->update_Hf();

    this->update_q(fiducial_params, q_index);
    //this->update_q_prime();
    
    //this->prefactor_Ml = 2*this->b_bias * this->c / this->pi;
    cout << "... Dependencies calculated ..." << endl;
    cout << "... precalculating inverse r ..." << endl;
    //this->update_r_inverse();
    cout << "... r inverse is calculated ..." << endl;
    cout << "... Initializing Pk interpolator ..." << endl;
    //this->update_Pk_interpolator(this->fiducial_params);
    
    this->update_Pk_interpolator_direct(this->fiducial_params, Pk_index);
    
    //this->update_Pk_interpolator_full(this->fiducial_params);
    cout << "... Pks calculated ..." << endl;

    cout << "... Creating Bessels ..." << endl;
    //this->create_bessel_interpolant_ALGLIB(0, this->fiducial_params["l_max"]);
    //this->create_bessel_interpolant_OWN(this->fiducial_params["l_min"],this->fiducial_params["l_max"]);
    cout << "... Bessels built ..." << endl;
    
    //set to true if ARES should be used,
    //set to false if G21 should be used.
    use_ARES = ARES_MACRO;

    cout << "... generating 21cm interface ..." << endl;
    if (!use_ARES) {
        G21 = new Global21cmInterface();
        this->update_G21(fiducial_params, Tb_index);
    } 
    else {
        ARES = new AresInterface();
        this->update_ARES(fiducial_params, Tb_index);
    }
    //this->update_G21_full(fiducial_params);
    cout << "... 21cm interface built ..." << endl;

    cout << "... generate analytic Tb ..." << endl;
    //this->update_Tb_analytic(fiducial_params, Tb_index);
    cout << "... Tb analytic done..." << endl;

    cout << "... CosmoCalc built ..." << endl;
}

bool CosmoCalc::get_useAres()
{
    return use_ARES;
}

double CosmoCalc::Cl(int l, double k1, double k2, double k_low, double k_high, int Pk_index, int Tb_index, int q_index)
{
    bool rsd = false;
    if (fiducial_params["rsd"] == 1.0)
        rsd = true;
    bool limber = false;
    if (fiducial_params["limber"] == 1.0)
        limber = true;

    if (rsd && !limber)
        return this->corr_Tb_rsd(l, k1, k2, k_low, k_high, Pk_index, Tb_index, q_index);
    else if (!rsd && !limber) 
        return this->corr_Tb(l, k1, k2, k_low, k_high, Pk_index, Tb_index, q_index);
    else if (rsd && limber)
        return this->Cl_limber_rsd(l,  k1, k2, Pk_index, Tb_index, q_index);
    else
        return this->Cl_limber(l, k1, k2, Pk_index, Tb_index, q_index);
}

double CosmoCalc::Cl_noise(int l, double k1, double k2)
{
    //TODO: integrand needs to be corrected.
    auto integrand = [&](double z)
    {
        double r;
        r = r_interp(z);
        double jl = sph_bessel_camb(l,k1*r);
        double hub = Hf_interp(z)*1000.0;
        return r*r*jl/hub; 
    };

    if (k1==k2) {
        // in mK
        double Tsys = fiducial_params["Tsys"];
        double fcover = fiducial_params["fcover"];
        double lmax = fiducial_params["lmax_noise"];
        // in seconds
        double tau = fiducial_params["tau_noise"];
        double prefactor = 2.0 *pi*c*c * Tsys*Tsys/(fcover*fcover * fiducial_params["df"] *\
                lmax * lmax * tau);
        double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
                this->zsteps_Ml);
        return prefactor * integral * integral;
    } else {
        return 0.0;
    }
}

CosmoCalc::~CosmoCalc()
{
    delete CAMB;
    delete G21;
    delete ARES;
}

void CosmoCalc::show_cosmo_calcs()
{
    cout << hubble_time() << endl;
    cout << hubble_dist() << endl;
    cout << comoving_radial_dist(10) << endl;
    cout << "O_m = " << O_M << ", O_V = " << O_V << "." << endl;
    cout << "Age in Gigayears "<< age_of_universe(0) * pow(10,10) *\
        3.08568 / (365.25 * 24 * 3600) << endl;
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

double CosmoCalc::angular_diam_dist(double z, double z2)
{
    const double Ok = 1 - this->O_M - this->O_V;
    const double root = sqrt(abs(1-this->O_tot));
    double result;
    if (z2 < 0.0) {
        if (this->O_tot == 1.0)
            result = this->D_H * this->Z(z)/(1+z);
        else
            result = this->D_H * this->S_k(root * this-> Z(z)) / ((1+z) * root);
    } else if (this->O_V + this->O_M <= 1.0) {
        const double dm = this->D_M(z);
        const double dm2 = this->D_M(z2);
        result = 1.0 / (1.0 + z2) *\
                 ( dm2 * sqrt(1 + Ok * pow(dm,2) / pow(this->D_H, 2)) -\
                   dm * sqrt(1 + Ok * pow(dm2,2) / pow(this->D_H,2)) );
    } else {
        cout << "Error: D_A12 formula invalid for O_tot > 1.0" << endl;
        result = 1.0;
    }
    return result;
}

double CosmoCalc::D_A(double z, double z2)
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
    auto integrand = [&](double x) 
    {
        return pow(1+x,2) * pow(this->D_A(x),2)/this->E(x);
    };
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
    auto integrand = [&](double x)
    {
        return 1.0/((1+x)*this->E(x));
    };
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
    return this->H(z) * 1000.0 / (3.08567758 * pow(10, 16));
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
    const double num = 4.0/3.0 * pi * pow(this->c/this->H_SI(0),3) *\
                       this->rho_crit(0);
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

void CosmoCalc::update_Hf()
{
    double z;
    real_1d_array xs, ys;
    xs.setlength(10*zmax_interp + 1);
    ys.setlength(10*zmax_interp + 1);

    for (int i = 0; i <= 10*zmax_interp; ++i) {
        z = i * 0.1;
        xs[i] = z;
        ys[i] = this->H(z);
    }
    spline1dbuildlinear(xs,ys,this->H_f_interp_full);
}
void CosmoCalc::update_r_inverse()
{
    this->r_inv.clear();
    real_1d_array xs, ys;
    xs.setlength(10*zmax_interp + 1);
    ys.setlength(10*zmax_interp + 1);

    double z;
    for (int n = 0; n <= 10*zmax_interp; ++n) {
        z =  n * 0.1;
        this->r_inv.push_back(this->D_C(z));
        xs[n] = z;
        ys[n] = r_inv[n];
    }
    spline1dbuildlinear(ys,xs,rinv_interp);
}

double CosmoCalc::r_inverse(double r)
{
    return spline1dcalc(rinv_interp, r);
}

void CosmoCalc::update_q_full()
{
    real_1d_array xs, ys;
    xs.setlength(10*zmax_interp + 1);
    ys.setlength(10*zmax_interp + 1);

    double z;
    for (int n = 0; n <= 10*zmax_interp; ++n) {
        z =  n * 0.1;
        xs[n] = z;
        ys[n] = D_C(z);
    }
    spline1dbuildlinear(xs,ys,q_interp_full);
}

void CosmoCalc::update_q_prime_full()
{
    double z;
    double h = 10e-4;
    real_1d_array xs1, ys1;
    xs1.setlength(10*zmax_interp + 1);
    ys1.setlength(10*zmax_interp + 1);

    for (int n = 0; n <= 10*zmax_interp; ++n) {
        z = 0.1 + n * 0.1;
        double res = 0;
        res = - D_C(z+2*h) + 8 * D_C(z+h) - 8 * D_C(z-h) + D_C(z-2*h);
        res = abs(res);
        xs1[n] = z;
        ys1[n] = res/(12*h);
    }
    spline1dbuildlinear(xs1,ys1,q_p_interp_full);
}

void CosmoCalc::update_q(map<string,double> params, int *q_index)
{
    bool limber = false;
    if (params["limber"] == 1.0)
        limber = true;

    // We first update q and then q' if the limber approximation is being used..
    bool do_calc = true;
    for (unsigned int i = 0; i < qs.size(); ++i) {
        if (params["ombh2"] == qs[i].ombh2 && params["omnuh2"] == qs[i].omnuh2 &&\
                params["omch2"] == qs[i].omch2 && params["omk"] == qs[i].omk &&\
                params["hubble"] == qs[i].hubble && params["T_CMB"] == qs[i].t_cmb &&\
                params["w_DE"] == qs[i].w_DE) {

            do_calc = false;
            *q_index = i;
            break;
        }
    }

    if (do_calc) {
        q_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.t_cmb = params["T_CMB"];
        interp.w_DE = params["w_DE"];
        // TODO: Do this in a way that works with parallelism....
        // UPDATE D_C to use the above parameters.

        double T_CMB2 = params["T_CMB"];
        double H_02 = params["hubble"];
        double h2 = H_02 / 100.0;
        double O_b2 = params["ombh2"] / pow(h2,2);
        double O_cdm2 = params["omch2"] / pow(h2,2);
        double O_nu2 = params["omnuh2"] / pow(h2,2);
        double O_gamma2 = pow(pi,2) * pow(T_CMB2/11605.0,4) /\
                          (15.0*8.098*pow(10,-11)*pow(h2,2));
        double O_nu_rel2 = O_gamma2 * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
        double O_R2 = O_gamma2 + O_nu_rel2;
        double O_k2 = params["omk"];
        double O_M2 = O_b2 + O_cdm2 + O_nu2;
        double O_tot2 = 1.0 - O_k2;
        double O_V2 = O_tot2 - O_M2 - O_R2;
        double D_H2 = c / (1000.0 * H_02);
        double w2 = params["w_DE"];

        real_1d_array xs, ys, qps, hs;
        xs.setlength(this->zsteps_Ml+1);
        ys.setlength(this->zsteps_Ml+1);
        qps.setlength(this->zsteps_Ml+1);
        hs.setlength(this->zsteps_Ml+1);
        double h = 10e-4;
        double z;
        for (int n = 0; n <= this->zsteps_Ml; ++n) {
            z = this->zmin_Ml + n * this->stepsize_Ml;
            xs[n] = z;

            auto integrand = [&](double zp)
            {
                return 1/sqrt(O_V2 * pow(1+zp,3*(1+w2)) + O_R2 * pow(1+zp,4) +\
                        O_M2 * pow(1+zp,3) + O_k2 * pow(1+zp,2));
            };
            double Z = integrate(integrand, 0.0, z, 1000, simpson());
            
            if (limber) {
                double dc1 = integrate(integrand, 0.0, z+2*h, 1000, simpson());
                double dc2 = integrate(integrand, 0.0, z+h, 1000, simpson());
                double dc3 = integrate(integrand, 0.0, z-h, 1000, simpson());
                double dc4 = integrate(integrand, 0.0, z-2*h, 1000, simpson());
                qps[n] =  abs(D_H2 * (-dc1 + 8 * dc2 - 8 * dc3 + dc4));
            }

            ys[n] = D_H2 * Z;
            hs[n] = H_02 * sqrt(O_V2 * pow(1+z,3*(1+w2)) + O_R2 * pow(1+z,4) +\
                    O_M2 * pow(1+z,3) + O_k2 * pow(1+z,2));
        }
        spline1dinterpolant interpolator, interpolator_Hf, interpolator_qp;
        spline1dbuildlinear(xs,ys,interpolator);
        spline1dbuildlinear(xs,hs,interpolator_Hf);
        spline1dbuildlinear(xs,qps,interpolator_qp);

        // If limber == false, the qp_interpolator will just be empty but that
        // is fine because it won't be used in that case.
        interp.h = h2;
        interp.interpolator = interpolator;
        interp.interpolator_Hf = interpolator_Hf;
        interp.interpolator_qp = interpolator_qp;

        qs.push_back(interp);
        *q_index = qs.size() - 1;
    }    
}

double CosmoCalc::Hf_interp(double z)
{
    return spline1dcalc(qs[0].interpolator_Hf,z);
}

double CosmoCalc::q_interp(double z, int q_index)
{
    return spline1dcalc(qs[q_index].interpolator,z);
}

double CosmoCalc::qp_interp(double z, int q_index)
{
    return spline1dcalc(qs[q_index].interpolator_qp,z);
}

double CosmoCalc::r_interp(double z)
{
    return spline1dcalc(qs[0].interpolator,z);
}

void CosmoCalc::update_q_prime()
{
    this->q_p_Ml.clear();
    double z;
    double h = 10e-4;

    real_1d_array xs, ys;
    xs.setlength(this->zsteps_Ml+1);
    ys.setlength(this->zsteps_Ml+1);

    for (int n = 0; n <= this->zsteps_Ml; ++n) {
        z = this->zmin_Ml + n * this->stepsize_Ml;
        double res = 0;
        res = - D_C(z+2*h) + 8 * D_C(z+h) - 8 * D_C(z-h) + D_C(z-2*h);
        res = abs(res);

        this->q_p_Ml.push_back(res/(12*h));

        xs[n] = z;
        ys[n] = q_p_Ml[n];
    }
    spline1dbuildlinear(xs,ys,q_p_interp);
}

void CosmoCalc::create_bessel_interpolant_ALGLIB(int lmin, int lmax)
{
    double xmin, xmax;
    this->lmin_bess = lmin;
    // TODO: these values are kind of arbitrary...
    // a little under (0.5)
    xmin = 0.5 * 0.001 * r_interp(this->fiducial_params["zmin"]);
    // a little over (2)
    xmax = 2 * 5 * r_interp(this->fiducial_params["zmax"]);
    //cout << "xmin for bessels is: " << xmin << ". xmax is: " << xmax << endl;

    for (int l = lmin; l <= lmax; ++l) {

        double stepsize = 1.0;

        int nsteps = (int)((xmax - xmin)/stepsize);
        //cout << "nsteps is: " << nsteps << endl; 
        real_1d_array ys, xs;
        xs.setlength(nsteps);
        ys.setlength(nsteps);
        for (int j = 0; j < nsteps; ++j) {
            double x = xmin + j * stepsize;
            xs[j] = x;
            ys[j] = this->sph_bessel_camb(l,x);
        }
        spline1dinterpolant interpolator;
        spline1dbuildcubic(xs, ys, interpolator);

        bessel_interp_list.push_back(interpolator);
    }
}

void CosmoCalc::create_bessel_interpolant_OWN(int lmin, int lmax)
{
    double xmax = 2 * 5 * r_interp(this->fiducial_params["zmax"]);
    this->lmin_bess = lmin;
    for (int l = lmin; l <= lmax; ++l) {
        vector<double> row;
        //finer grid now. (stepsize = 0.5)
        for (int j = 0; j < 2*(int)xmax; ++j) {
            row.push_back(this->sph_bessel_camb(l,(double)j / 2.0));
        }
        bessel_values.push_back(row);  
    }
}

double CosmoCalc::bessel_j_interp(int l, double x)
{
    return spline1dcalc(bessel_interp_list[l - this->lmin_bess], x);
}

double CosmoCalc::bessel_j_interp_cubic(int l, double x)
{
    int l_used = l - this->lmin_bess;
    if ((unsigned int)(2*x) + 2 >= bessel_values[l_used].size()) {
        return bessel_values[l_used][bessel_values[l_used].size()-1];
    } else if (int(2*x) - 1 < 0) {
        return bessel_values[l_used][0];
    } else {
        if (x-(int)x == 0 || x-(int)x == 0.5) {       
            return bessel_values[l_used][(int)(2 * x)];
        } else {
            double y0, y1, y2, y3, a0, a1, a2, a3, mu, mu2, x0;

            x0 = (int)(2*x);
            x0 = x0/2;
            mu = (x - x0)/0.5;
            mu2 = mu*mu;

            y0 = bessel_values[l_used][(int)(2 * x) - 1];
            y1 = bessel_values[l_used][(int)(2 * x)];
            y2 = bessel_values[l_used][(int)(2 * x) + 1];
            y3 = bessel_values[l_used][(int)(2 * x) + 2];

            a0 = 0.5*y3 - 1.5 * y2 - 0.5 * y0 + 1.5 * y1;
            a1 = y0 - 2.5*y1 +2*y2 - 0.5*y3;
            a2 = 0.5*y2 - 0.5*y0;
            a3 = y1;

            return a0*mu*mu2 + a1*mu2 + a2*mu + a3;
        }
    }
}

double CosmoCalc::bessel_j_interp_basic(int l, double x)
{
    if ((unsigned int)(2*x) + 1 >= bessel_values[l].size()) {
        return bessel_values[l][bessel_values[l].size()-1];
    } else {
        if (x-(int)x == 0 || x-(int)x == 0.5) {       
            return bessel_values[l][(int)(2 * x)];
        } else {
            double two_x0 = (int)(2*x);
            return bessel_values[l][(int)(2 * x)] +\
                (bessel_values[l][(int)(2 * x) + 1] -\
                 bessel_values[l][(int)(2 * x)]) * (x - (two_x0)/2) / 0.5;
        }
    }
}

void CosmoCalc::update_G21_full(map<string,double> params, int *Tb_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Tbs.size(); ++i) {
        if (params["ombh2"] == Tbs_full[i].ombh2 && params["omnuh2"] == Tbs_full[i].omnuh2 &&\
                params["omch2"] == Tbs_full[i].omch2 && params["omk"] == Tbs_full[i].omk &&\
                params["hubble"] == Tbs_full[i].hubble && params["sigma8"] == Tbs_full[i].s8 &&\
                params["T_CMB"] == Tbs_full[i].T_CMB && params["n_s"] == Tbs_full[i].n_s &&\
                params["fstar"] == Tbs_full[i].fstar && params["fesc"] == Tbs_full[i].fesc &&\
                params["nion"] == Tbs_full[i].nion && params["fx"] == Tbs_full[i].fx &&\
                params["flya"] == Tbs_full[i].flya && params["w_DE"] == Tbs_full[i].w_DE) {

            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {
        Tb_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.s8 = params["sigma8"];
        interp.T_CMB = params["T_CMB"];
        interp.n_s = params["n_s"];
        interp.fstar = params["fstar"];
        interp.fesc = params["fesc"];
        interp.nion = params["nion"];
        interp.fx = params["fx"];
        interp.flya = params["flya"];
        interp.w_DE = params["w_DE"];

        cout << "G21 is being updated" << endl;
        G21->updateGlobal21cm_full(params);
        vector<double> vz, vTb;
        G21->getTb(&vz, &vTb);

        real_1d_array g21_z, g21_Tb;
        g21_z.setlength(vz.size());
        g21_Tb.setlength(vTb.size());

        for (unsigned int i = 0; i < vz.size(); i++){
            g21_z[i] = vz[i];
        }
        for (unsigned int i = 0; i < vTb.size(); i++){
            g21_Tb[i] = vTb[i];
        }

        spline1dinterpolant interpolator;
        spline1dbuildcubic(g21_z, g21_Tb, interpolator);
        interp.interpolator = interpolator;

        Tbs_full.push_back(interp);
        *Tb_index = Tbs_full.size() - 1;

    }
}

double CosmoCalc::Tb_interp_full(double z, int Tb_index)
{
    // The * 1000.0 is so we get the result in mK
    return spline1dcalc(Tbs_full[Tb_index].interpolator,z) * 1000.0;
}

void CosmoCalc::update_G21(map<string,double> params, int *Tb_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Tbs.size(); ++i) {
        if (params["ombh2"] == Tbs[i].ombh2 && params["omnuh2"] == Tbs[i].omnuh2 &&\
                params["omch2"] == Tbs[i].omch2 && params["omk"] == Tbs[i].omk &&\
                params["hubble"] == Tbs[i].hubble && params["sigma8"] == Tbs[i].s8 &&\
                params["T_CMB"] == Tbs[i].T_CMB && params["n_s"] == Tbs[i].n_s &&\
                params["fstar"] == Tbs[i].fstar && params["fesc"] == Tbs[i].fesc &&\
                params["nion"] == Tbs[i].nion && params["fx"] == Tbs[i].fx &&\
                params["flya"] == Tbs[i].flya && params["w_DE"] == Tbs[i].w_DE) {
            cout << "found precalculated G21" << endl;
            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {
        Tb_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.s8 = params["sigma8"];
        interp.T_CMB = params["T_CMB"];
        interp.n_s = params["n_s"];
        interp.fstar = params["fstar"];
        interp.fesc = params["fesc"];
        interp.nion = params["nion"];
        interp.fx = params["fx"];
        interp.flya = params["flya"];
        interp.w_DE = params["w_DE"];

        cout << "G21 is being updated" << endl;
        G21->updateGlobal21cm(params);
        vector<double> vz, vTb;
        G21->getTb(&vz, &vTb);

        real_1d_array g21_z, g21_Tb;
        g21_z.setlength(vz.size());
        g21_Tb.setlength(vTb.size());

        for (unsigned int i = 0; i < vz.size(); i++){
            g21_z[i] = vz[i];
        }
        for (unsigned int i = 0; i < vTb.size(); i++){
            g21_Tb[i] = vTb[i];
        }

        spline1dinterpolant interpolator;
        spline1dbuildcubic(g21_z, g21_Tb, interpolator);
        interp.interpolator = interpolator;

        Tbs.push_back(interp);
        *Tb_index = Tbs.size() - 1;

        cout << "G21 update done" << endl;
    }
}

double CosmoCalc::Tb_interp(double z, int Tb_index)
{
    if (!use_ARES) {
        // The * 1000.0 is so we get the result in mK
        return spline1dcalc(Tbs[Tb_index].interpolator,z) * 1000.0;
    }
    else {
        return spline1dcalc(Tbs_ares[Tb_index].interpolator,z);
    }
}

void CosmoCalc::update_ARES(map<string,double> params, int *Tb_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Tbs_ares.size(); ++i) {
        if (params["ombh2"] == Tbs_ares[i].ombh2 &&\
            params["omnuh2"] == Tbs_ares[i].omnuh2 &&\
            params["omch2"] == Tbs_ares[i].omch2 &&\ 
            params["omk"] == Tbs_ares[i].omk &&\
            params["hubble"] == Tbs_ares[i].hubble &&\
            params["sigma8"] == Tbs_ares[i].s8 &&\
            params["T_CMB"] == Tbs_ares[i].T_CMB &&\
            params["n_s"] == Tbs_ares[i].n_s &&\
            params["fstar"] == Tbs_ares[i].fstar &&\
            params["fesc"] == Tbs_ares[i].fesc &&\
            params["nion"] == Tbs_ares[i].nion &&\
            params["fx"] == Tbs_ares[i].fX ) {
                
            /* **** These parameters aren't part of the fiducial parameter
             * set, or, as is the case for w_DE, aren't used by ARES.
                params["Tmin"] == Tbs_ares[i].Tmin &&\
                params["w_DE"] == Tbs_ares[i].w_DE &&\
                params["Nlw"] == Tbs_ares[i].Nlw &&\
                params["cX"] == Tbs_ares[i].cX &&\
                params["HeByMass"] == Tbs_ares[i].HeByMass
            */
            cout << "found precalculated Ares" << endl;
            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {
        Tb_interpolator_ares interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.s8 = params["sigma8"];
        interp.T_CMB = params["T_CMB"];
        interp.n_s = params["n_s"];
        interp.fstar = params["fstar"];
        interp.fesc = params["fesc"];
        interp.nion = params["nion"];
        interp.fX = params["fx"];

        
        interp.w_DE = -1; //params["w_DE"];
        interp.Tmin = -1; //params["Tmin"];
        interp.Nlw = -1; //params["Nlw"];
        interp.cX = -1; //params["cX"];
        interp.HeByMass = -1; //params["HeByMass"];
        

        cout << "Ares is being updated" << endl;
        ARES->updateAres(params);
        vector<double> vz, vTb;
        ARES->getTb(&vz, &vTb);

        real_1d_array Ares_z, Ares_Tb;
        Ares_z.setlength(vz.size());
        Ares_Tb.setlength(vTb.size());

        for (unsigned int i = 0; i < vz.size(); i++){
            Ares_z[i] = vz[i];
        }
        for (unsigned int i = 0; i < vTb.size(); i++){
            Ares_Tb[i] = vTb[i];
        }

        spline1dinterpolant interpolator;
        spline1dbuildcubic(Ares_z, Ares_Tb, interpolator);
        interp.interpolator = interpolator;

        Tbs_ares.push_back(interp);
        *Tb_index = Tbs.size() - 1;

        cout << "Ares dTb update done" << endl;
    }
}

double CosmoCalc::Tb_interp_ARES(double z, int Tb_index)
{
    return spline1dcalc(Tbs_ares[Tb_index].interpolator,z);
}

void CosmoCalc::update_Pk_interpolator_full(map<string, double> params, int *Pk_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Pks_full.size(); ++i) {
        if (params["ombh2"] == Pks_full[i].ombh2 && params["omnuh2"] == Pks_full[i].omnuh2 &&\
                params["omch2"] == Pks_full[i].omch2 && params["omk"] == Pks_full[i].omk &&\
                params["hubble"] == Pks_full[i].hubble && params["w_DE"] == Pks_full[i].w_DE) {

            do_calc = false;
            *Pk_index = i;
            break;
        }
    }


    if (do_calc) {
        Pk_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.w_DE = params["w_DE"];

        CAMB->call_full(params);    
        vector<double> vk = CAMB->get_k_values();
        vector<vector<double>> Pz = CAMB->get_Pz_values();

        double z_stepsize = 1;
        vector<double> vz, vP;
        for (unsigned int i = 0; i < Pz.size(); ++i) {
            vz.push_back(i * z_stepsize);
            vP.insert(vP.end(), Pz[i].begin(), Pz[i].end());
        }

        real_1d_array matterpowerspectrum_k, matterpowerspectrum_z, matterpowerspectrum_P;
        matterpowerspectrum_k.setlength(vk.size());
        matterpowerspectrum_z.setlength(vz.size());
        matterpowerspectrum_P.setlength(vP.size());
        for (unsigned int i = 0; i < vk.size(); i++){
            matterpowerspectrum_k[i] = vk[i];
        }
        for (unsigned int i = 0; i < vP.size(); i++){
            matterpowerspectrum_P[i] = vP[i];
        }
        for (unsigned int i = 0; i < vz.size(); i++){
            matterpowerspectrum_z[i] = vz[i];
        }

        spline2dinterpolant interpolator;
        spline2dbuildbilinearv(matterpowerspectrum_k, vk.size(),matterpowerspectrum_z, vz.size(),\
                matterpowerspectrum_P, 1, interpolator);
        interp.interpolator = interpolator;

        Pks_full.push_back(interp);
        *Pk_index = Pks_full.size() - 1;
    }
}
double CosmoCalc::Pk_interp_full(double k, double z, int Pk_index)
{
    return spline2dcalc(Pks_full[Pk_index].interpolator, k, z);
}


void CosmoCalc::update_Pk_interpolator_direct(map<string, double> params, int *Pk_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Pks.size(); ++i) {
        if (params["ombh2"] == Pks[i].ombh2 && params["omnuh2"] == Pks[i].omnuh2 &&\
                params["omch2"] == Pks[i].omch2 && params["omk"] == Pks[i].omk &&\
                params["hubble"] == Pks[i].hubble && params["T_CMB"] == Pks[i].tcmb &&\
                params["w_DE"] == Pks[i].w_DE && params["n_s"] == Pks[i].n_s &&\
                params["A_s"] == Pks[i].A_s){

            do_calc = false;
            *Pk_index = i;
            break;
        }
    }


    if (do_calc) {
        Pk_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.omk = params["omk"];
        interp.hubble = params["hubble"];
        interp.tcmb = params["T_CMB"];
        interp.w_DE = params["w_DE"];
        interp.n_s = params["n_s"];
        interp.A_s = params["A_s"];


        CAMB->call(params);    
        vector<double> vk = CAMB->get_k_values();
        vector<vector<double>> Pz = CAMB->get_Pz_values();

        double z_stepsize = (params["zmax"] - params["zmin"])/(params["Pk_steps"] - 1);
        vector<double> vz, vP;
        for (unsigned int i = 0; i < Pz.size(); ++i) {
            vz.push_back(params["zmin"] + i * z_stepsize);
            vP.insert(vP.end(), Pz[i].begin(), Pz[i].end());
        }

        real_1d_array matterpowerspectrum_k, matterpowerspectrum_z, matterpowerspectrum_P;
        matterpowerspectrum_k.setlength(vk.size());
        matterpowerspectrum_z.setlength(vz.size());
        matterpowerspectrum_P.setlength(vP.size());
        for (unsigned int i = 0; i < vk.size(); i++){
            matterpowerspectrum_k[i] = vk[i];
        }
        for (unsigned int i = 0; i < vP.size(); i++){
            matterpowerspectrum_P[i] = vP[i];
        }
        for (unsigned int i = 0; i < vz.size(); i++){
            matterpowerspectrum_z[i] = vz[i];
        }

        spline2dinterpolant interpolator;
        spline2dbuildbilinearv(matterpowerspectrum_k, vk.size(),matterpowerspectrum_z, vz.size(),\
                matterpowerspectrum_P, 1, interpolator);
        interp.interpolator = interpolator;

        Pks.push_back(interp);
        *Pk_index = Pks.size() - 1;
    }
}

double CosmoCalc::Pk_interp(double k, double z, int Pk_index)
{
    return spline2dcalc(Pks[Pk_index].interpolator, k, z);
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
    auto integrand = [&](double x)
    {
        return (1+x)/pow(this->E(x),3);
    };
    return prefactor * integrate(integrand, z, pow(10, 5), 1000000, simpson());
}

double CosmoCalc::P_delta(double k, string units_k, string units_P)
{
    double keq, k_factor, delta_H;
    keq = 0;
    k_factor = 0;
    if (units_P == "default") {
        if (units_k == "default") {
            keq = 0.073 * this->O_M * this->h;
            k_factor = k;
        } else if (units_k == "Mpc-1" or units_k == "mpc-1") {
            keq = 0.073 * this->O_M * pow(this->h,2);
            k_factor = k / this->h;

        }
    } else if (units_P == "Mpc3" or units_P == "mpc3") {
        if (units_k == "default") {
            keq = 0.073 * this->O_M * this->h;
            k_factor = k / pow(this->h,3);

        } else if (units_k == "Mpc-1" or units_k == "mpc-1") {
            keq = 0.073 * this->O_M * pow(this->h,2);
            k_factor = k / pow(this->h,4);
        }
    } 

    if (this->O_M == 1) {
        delta_H = 1.9 * pow(10,-5);
    } else {
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
    const double bracket = 1 + 0.284 * x + pow(1.18 * x, 2) + pow(0.399 * x, 3) +\
                           pow(0.49 * x, 4);
    return res * pow(bracket, -0.25);
}

double CosmoCalc::Cl_new(int l, double k1, double k2, double k_low,\
        double k_high, int n_levin, int Pk_index, int Tb_index, int q_index)
{
    double a;
    double low;
    double hhh = pow(qs[q_index].h,3);
    if (l < 50){
        low = k_low;
        a = (double)(l+1000)/(1.5*10000);

    } else if (l < 1000){
        low = (double)l/(1.2*10000);
        a = (double)(l+1000)/(1.5*10000);
    } else {
        low = (double)l/(10000);
        a = (double)(l+1000)/(1.5*10000);
    }
    double lower_kappa_bound;
    if (low > k_low)
        lower_kappa_bound = low;
    else
        lower_kappa_bound = k_low;
    
   
    if (a < lower_kappa_bound)
        a += 0.1;
    if (a > k_high)
        a = k_high;
    int steps = (int)(abs(a - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;
    //cout << lower_kappa_bound << " " << a << endl;

    auto integrand1 = [&](double z)
    {
        double r,q;
        r = r_interp(z);
        q = q_interp(z, q_index);

        auto integrand2 = [&](double zp)
        {
            double rp,qp;
            rp = r_interp(zp);
            qp = q_interp(zp, q_index);

            auto integrand3 = [&](double kappa)
            {
                double sP = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
                double sPp = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
                return kappa*kappa * sP * sPp * this->sph_bessel_camb(l,kappa*q) *\
                    this->sph_bessel_camb(l,kappa*qp);
            };
            double integral = integrate_simps(integrand3, lower_kappa_bound, a, steps);

            if (a < k_high) {
                Levin LEVIN(a, k_high);

                auto foo = [&](double kappa)
                {
                    double sP1 = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
                    double sPp1 = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
                    return kappa*kappa * sP1 * sPp1;     
                };

                double integral3;
                if (z == zp)
                    integral3 = LEVIN.integrate_2sphj_1r(foo,q,l,n_levin);
                else
                    integral3 = LEVIN.integrate_2sphj_2r(foo,q,qp,l,n_levin);


                integral += integral3;
            }
            return rp*rp / (Hf_interp(zp)*1000.0) * this->Tb_interp(zp, Tb_index) *\
                this->sph_bessel_camb(l,k2*rp) * integral;
        };
        double integral2 = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
        return r*r / (Hf_interp(z)*1000.0) * this->Tb_interp(z, Tb_index) *\
            this->sph_bessel_camb(l,k1*r) * integral2;
    };
    double integral1 = integrate_simps(integrand1,this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return pow(this->prefactor_Ml,2) * integral1;
}

double CosmoCalc::Cl_new_analyticTb(int l, double k1, double k2, double k_low,\
        double k_high, int n_levin, int Pk_index, int Tb_index, int q_index)
{
    double a;
    double low;
    double hhh = pow(qs[q_index].h,3);

    auto integrand1 = [&](double z)
    {
        double r,q;
        r = r_interp(z);
        q = q_interp(z, q_index);

        if (l < 1000){
            low = (double)l/(1.2*q);
            a = (double)(l+1000)/(1.5*q);
        } else {
            low = (double)l/(q);
            a = (double)(l+1000)/(1.5*q);
        }
        double lower_kappa_bound;
        if (low > k_low)
            lower_kappa_bound = low;
        else
            lower_kappa_bound = k_low;

        int steps = (int)((a - lower_kappa_bound)/this->k_stepsize);
        if (steps % 2 == 1)
            ++steps;

        auto integrand2 = [&](double zp)
        {
            double rp,qp;
            rp = r_interp(zp);
            qp = q_interp(zp, q_index);

            auto integrand3 = [&](double kappa)
            {
                double sP = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
                double sPp = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
                return kappa*kappa * sP * sPp * this->sph_bessel_camb(l,kappa*q) *\
                    this->sph_bessel_camb(l,kappa*qp);
            };
            double integral = integrate_simps(integrand3, lower_kappa_bound, a, steps);
            double integral3;

            Levin LEVIN(a, k_high);

            auto foo = [&](double kappa)
            {
                double sP1 = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
                double sPp1 = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
                return kappa*kappa * sP1 * sPp1;     
            };


            if (z == zp)
                integral3 = LEVIN.integrate_2sphj_1r(foo,q,l,n_levin);
            else
                integral3 = LEVIN.integrate_2sphj_2r(foo,q,qp,l,n_levin);


            integral3 += integral;
            return rp*rp / (Hf_interp(zp)*1000.0) * this->Tb_analytic_interp(zp, Tb_index) *\
                this->sph_bessel_camb(l,k2*rp) * integral3;
        };
        double integral2 = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
        return r*r / (Hf_interp(z)*1000.0) * this->Tb_analytic_interp(z, Tb_index) *\
            this->sph_bessel_camb(l,k1*r) * integral2;
    };
    double integral1 = integrate_simps(integrand1,this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return pow(this->prefactor_Ml,2) * integral1;
}

double CosmoCalc::Cl_limber_rsd(int l, double k1, double k2, int Pk_index,\
        int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r,rr,q,qq,qp,k1r,k2r;
        r = r_interp(z);
        rr = r*r;
        q = q_interp(z, q_index);
        qq = q*q;
        qp = qp_interp(z, q_index);
        k1r = k1 * r;
        k2r = k2 * r;
        double hh = pow(Hf_interp(z)*1000.0, 2);
        double j1,j2,j3,j4;
        j1 = this->sph_bessel_camb(l,k1r);
        j2 = this->sph_bessel_camb(l,k2r);
        j3 = this->sph_bessel_camb(l-1,k1r);
        j4 = this->sph_bessel_camb(l-1,k2r);

        double Jl1 =  j1 * j2;
        vector<double> Jl2, Jl3, Jl4, L2, L3, L4;
        Jl2.push_back(Jl1);
        Jl2.push_back(j1 * j4);
        Jl3.push_back(Jl1);
        Jl3.push_back(j3 * j2);
        Jl4.push_back(Jl1);
        Jl4.push_back(Jl2[1]);
        Jl4.push_back(Jl3[1]);
        Jl4.push_back(j3*j4);
        double LL1 = -1.0/(double)(l*(2*l+1)*(2*l+1));
        //LL1 =- 1.0/(double)(300*300*300);
        L2.push_back(LL1 * (l+1));
        L2.push_back(-LL1 * k2r);
        L3.push_back(L2[0]);
        L3.push_back(-LL1 * k1r);
        double LL2 = 2.0*(-32.0*pow(l,6) - 24.0*pow(l,5) + 48.0*pow(l,4) +\
                46.0*pow(l,3) - 5.0*l - 1.0);
        LL2 = LL2 / (double)(pow(l,3)*pow(2*l-1,2)*pow(2*l+1,4));
        //LL2 = -1.0/(double)(300*300*300);
        L4.push_back(LL2 * pow(l+1,2));
        L4.push_back(-LL2 * k2r * (l+1));
        L4.push_back(-LL2 * k1r * (l+1));       
        L4.push_back(LL2 * k1r * k2r);
        double A = rr * pow(this->Tb_interp(z, Tb_index),2) *\
                   this->Pk_interp((double)l/q * qs[q_index].h,z,Pk_index)/\
                   (pow(qs[q_index].h,3)*hh*qp);   
        double JL2 = 0;
        double JL3 = 0;
        for (int i = 0; i < 2; i++) {
            JL2 += Jl2[i] * L2[i];
            JL3 += Jl3[i] * L3[i];
        }
        double JL4 = 0;
        for (int i = 0; i < 4; i++) {
            JL4 += Jl4[i] * L4[i];
        }
        double bb = this->b_bias * this->beta;
        double bb2 = bb*bb;
        double bracket = rr/qq * Jl1 + bb*r/((1+z)*q) * (JL2+JL3) +\
                         bb2/pow(1+z,2) * JL4; 

        return A * bracket;

    };
    double prefact = pow(this->prefactor_Ml,2) * this->pi / 2.0;
    return prefact * integrate_simps(integrand, this->zmin_Ml,\
            this->zmax_Ml, this->zsteps_Ml);
}

double CosmoCalc::Cl_limber(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r,q,qp,rr;
        r = r_interp(z);
        q = q_interp(z, q_index);
        qp = qp_interp(z, q_index);
        rr = r*r;
        double hh = pow(Hf_interp(z)*1000.0, 2);
        double A = rr * this->Pk_interp(((double)l + 0.5)/q * qs[q_index].h,z,Pk_index)/(pow(qs[q_index].h,3)*hh*qp) *\
                   pow(this->Tb_interp(z, Tb_index),2);

        //TODO: check whether we need to multiply py h.
        return A * rr / (q*q) * this->sph_bessel_camb(l,k1*r) * this->sph_bessel_camb(l, k2*r);
    };

    double pre = pow(this->prefactor_Ml,2) * this->pi / 2.0;
    return  pre * integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
            this->zsteps_Ml);
}

double CosmoCalc::Cl_simplified_levin(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r, q, qp;
        r = r_interp(z);
        q = q_interp(z,q_index);
        qp = spline1dcalc(q_p_interp,z);
        double rr = r*r;
        double hh = pow(Hf_interp(z)*1000.0, 2);
        double A = rr * this->Pk_interp(((double)l + 0.5)/q * qs[q_index].h,z, Pk_index)/(pow(qs[q_index].h,3)*hh*qp) *\
                   pow(this->Tb_interp(z, Tb_index),2);

        //TODO: check whether we need to multiply py h.
        // here: changed it back to non - interpolation, because it isn't necessary for the simplified case.
        return A * rr / (q*q) * this->sph_bessel_camb(l,k1*r) * this->sph_bessel_camb(l, k2*r);
    };

    double pre = pow(this->prefactor_Ml,2) * this->pi / 2.0;
    return pre * qromb(integrand,this->zmin_Ml,this->zmax_Ml,1.0E-5);
}


double CosmoCalc::corr_Tb(int l, double k1, double k2, double k_low,\
        double k_high, int Pk_index, int Tb_index, int q_index)
{
    //This determines the lower bound of the kappa integral
    double low;
    if (l < 50){
        low = k_low;
    } else if (l < 1000){
        low = (double)l/(1.2*10000);
    } else {
        low = (double)l/(10000);
    }
    double lower_kappa_bound;// = k_low;
    if (low > k_low)
        lower_kappa_bound = low;
    else
        lower_kappa_bound = k_low;
    
    //This determines the upper bound of the kappa integral
    double higher_kappa_bound = max(k1,k2) + 0.1;
    
    int steps = (int)(abs(higher_kappa_bound - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;

    //The boundaries and thus the #steps is determined by the properties of Bessel functions.
    if (k1 == k2)
    {
        auto integrand = [&](double kappa)
        {
            return pow(kappa,2) * pow(this->M(l,k1,kappa,Pk_index,Tb_index,q_index),2);
        };

        //return integrate(integrand, k_low, k_high, this->k_steps, simpson());
        return integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);
    } else {
        auto integrand = [&](double kappa)
        {
            return pow(kappa,2) * this->M(l,k1,kappa,Pk_index,Tb_index,q_index) *\
                this->M(l,k2,kappa,Pk_index,Tb_index,q_index);
        };

        //return integrate(integrand, k_low, k_high, this->k_steps, simpson());
        return integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);
    }
}

double CosmoCalc::Cl_simplified2(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index)
{
    double res2 = 2 * pow(this->b_bias,2) * pow(this->c,2)/this->pi;
    double zstar1 = r_inverse((l+0.5)/k1);
    double zstar2 = r_inverse((l+0.5)/k2);

    double H1 = spline1dcalc(H_f_interp_full, zstar1)*1000.0;
    double H2 = spline1dcalc(H_f_interp_full, zstar2)*1000.0;
    double Tb1 = Tb_interp_full(zstar1, Tb_index);
    double Tb2 = Tb_interp_full(zstar2, Tb_index);
    double q1 = spline1dcalc(q_interp_full, zstar1);
    double q2 = spline1dcalc(q_interp_full, zstar2);
    double qp1 = spline1dcalc(q_p_interp_full, zstar1);
    double qp2 = spline1dcalc(q_p_interp_full, zstar2);

    double hhh = pow(qs[q_index].h,3);
    auto integrand = [&](double kappa)
    {
        double sP1 = sqrt(this->Pk_interp_full(kappa*qs[q_index].h,zstar1, Pk_index)/hhh);
        double sP2 = sqrt(this->Pk_interp_full(kappa*qs[q_index].h,zstar2, Pk_index)/hhh);

        double jl1 = sph_bessel_camb(l,q1*kappa);
        double jl2 = sph_bessel_camb(l,q2*kappa);

        return pow(kappa,2) * sP1*sP2*jl1*jl2;
    };
    double integral = integrate_simps(integrand, 0.1,1,100000);
    res2 = res2*pow(l+0.5,3)/(k1*k1*k1*k2*k2*k2* H1 * H2 * qp1 * qp2) * Tb1 * Tb2 * integral;
    return res2;

}

double CosmoCalc::Cl_simplified3(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index)
{   
    double hhh = pow(qs[q_index].h,3);
    auto integrand1 = [&] (double zp)
    {
        double rp,qp;
        rp = r_interp(zp);
        qp = q_interp(zp, q_index);


        auto integrand2 = [&](double z)
        {
            double r,q;
            r = r_interp(z);
            q = q_interp(z, q_index);

            double sP = sqrt(this->Pk_interp((l+0.5)/q*qs[q_index].h,z, Pk_index)/hhh);
            double sPp = sqrt(this->Pk_interp((l+0.5)/q*qs[q_index].h,zp, Pk_index)/hhh);
            double jqq = this->bessel_j_interp_cubic(l,(l+0.5)*qp/q);
            return r*r / (q*q*q * Hf_interp(z)*1000.0) * this->Tb_interp(z, Tb_index) *\
                this->bessel_j_interp_cubic(l,k1*r) * jqq * sP * sPp;
        };
        double integral = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
        return rp*rp / (Hf_interp(zp)*1000.0) * this->Tb_interp(zp, Tb_index) * this->bessel_j_interp_cubic(l,k2*rp) * integral;
    };
    double integral2 = integrate_simps(integrand1, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return pow(this->prefactor_Ml,2) * sqrt(this->pi/(2*(l+0.5))) * pow((l+0.5),2) * integral2;
}

double CosmoCalc::corr_Tb_rsd(int l, double k1, double k2, double k_low,\
        double k_high, int Pk_index, int Tb_index, int q_index)
{       
    double low;
    if (l < 50){
        low = k_low;
    } else if (l < 1000){
        low = (double)l/(1.2*10000);
    } else {
        low = (double)l/(10000);
    }
    double lower_kappa_bound;
    if (low > k_low)
        lower_kappa_bound = low;
    else
        lower_kappa_bound = k_low;

    double higher_kappa_bound = max(k1,k2) + 0.1;

    int steps = (int)(abs(higher_kappa_bound - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;
    

    auto integrand = [&](double k)
    {
        double m1,n1,m2,n2;
        m1 = this->M(l,k1,k,Pk_index,Tb_index,q_index);
        n1 = this->N_bar(l,k1,k,Pk_index,Tb_index,q_index);
        if (k1 == k2) {
            m2 = m1;
            n2 = n1;
        } else {
            m2 = this->M(l,k2,k,Pk_index,Tb_index,q_index);
            n2 = this->N_bar(l,k2,k,Pk_index,Tb_index,q_index);
        }
        const double bb = this->b_bias * this->beta;
        const double bb2 = pow(bb,2);

        return pow(k,2) * m1 * m2 + bb * k * (m1*n2 + n1*m2) + bb2 * n1 * n2;
    };
    //return integrate(integrand, k_low, k_high, this->k_steps, simpson());

    return integrate_simps(integrand, lower_kappa_bound, higher_kappa_bound, steps);
}

// This is the same as corr_Tb just reordered the integrals.
double CosmoCalc::corr_Tb_new(int l, double k1, double k2, double k_low,\
        double k_high, int Pk_index, int Tb_index, int q_index)
{
    double low;
    double hhh = pow(qs[q_index].h,3);
    if (l < 50){
        low = k_low;
    } else if (l < 1000){
        low = (double)l/(1.2*10000);
    } else {
        low = (double)l/(10000);
    }
    double lower_kappa_bound;
    if (low > k_low)
        lower_kappa_bound = low;
    else
        lower_kappa_bound = k_low;

    int steps = (int)(abs(k_high - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;

    auto integrand1 = [&](double z)
    {
        double r,q;
        r = r_interp(z);
        q = q_interp(z, q_index);


        auto integrand2 = [&](double zp)
        {
            double rp,qp;
            rp = r_interp(zp);
            qp = q_interp(zp, q_index);


            auto integrand3 = [&](double kappa)
            {
                double hhh = pow(qs[q_index].h,3);
                double sP = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
                double sPp = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
                return kappa*kappa * sP * sPp * this->sph_bessel_camb(l,kappa*q) *\
                    this->sph_bessel_camb(l,kappa*qp);
            };
            double integral3 = integrate_simps(integrand3, k_low, k_high, steps);
            return rp*rp / (Hf_interp(zp)*1000.0) * this->Tb_interp(zp, Tb_index) *\
                this->sph_bessel_camb(l,k2*rp) * integral3;
        };
        double integral2 = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
        return r*r / (Hf_interp(z)*1000.0) * this->Tb_interp(z, Tb_index) *\
            this->sph_bessel_camb(l,k1*r) * integral2;
    };
    double integral1 = integrate_simps(integrand1,this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return pow(this->prefactor_Ml,2) * integral1;
}

//same as corr_Tb_new but using qromb as integration method
double CosmoCalc::corr_Tb_new2(int l, double k1, double k2, double k_low,\
        double k_high, int Pk_index, int Tb_index, int q_index)
{
    int steps = (int)((k_high - k_low)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;

    auto integrand1 = [&](double z)
    {
        double hh = Hf_interp(z)*1000.0;
        double r,q;
        r = r_interp(z);
        q = q_interp(z,q_index);

        auto integrand2 = [&](double zp)
        {
            double hhp = Hf_interp(zp)*1000.0;
            double rp,qp;
            rp = r_interp(zp);
            qp = q_interp(zp,q_index);

            auto integrand3 = [&](double kappa)
            {
                double hhh = pow(qs[q_index].h,3);
                double sP = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
                double sPp = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
                return kappa*kappa * sP * sPp * this->bessel_j_interp_cubic(l,kappa*q) *\
                    this->bessel_j_interp_cubic(l,kappa*qp);
            };
            double integral3 = qromb(integrand3, k_low, k_high, 1.0E-6);
            return rp*rp / hhp * this->Tb_interp(zp, Tb_index) *\
                this->bessel_j_interp_cubic(l,k2*rp) * integral3;
        };
        double integral2 = qromb(integrand2, this->zmin_Ml, this->zmax_Ml, 1.0E-2);
        return r*r / hh * this->Tb_interp(z, Tb_index) *\
            this->bessel_j_interp_cubic(l,k1*r) * integral2;
    };
    double integral1 = qromb(integrand1,this->zmin_Ml, this->zmax_Ml, 1.0E-1);
    return pow(this->prefactor_Ml,2) * integral1;
}

void CosmoCalc::compare(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index)
{
    double k_high = 0.6;
    double k_low = 0.4;
    int steps = (int)((k_high - k_low)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;
    double z = 8.0;

    double r,q;
    r = r_interp(z);
    q = q_interp(z, q_index);


    auto integrand2 = [&](double zp)
    {
        double rp,qp;
        rp = r_interp(zp);
        qp = q_interp(zp, q_index);


        auto integrand3 = [&](double kappa)
        {
            double hhh = pow(qs[q_index].h,3);
            double sP = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
            double sPp = sqrt(this->Pk_interp(kappa*qs[q_index].h,zp, Pk_index)/hhh);
            return kappa*kappa * sP * sPp * this->bessel_j_interp_cubic(l,kappa*q) *\
                this->bessel_j_interp_cubic(l,kappa*qp);
        };
        double integral3 = integrate_simps(integrand3, k_low, k_high, steps);
        return rp*rp / (Hf_interp(zp)*1000.0) * this->Tb_interp(zp, Tb_index) *\
            this->bessel_j_interp_cubic(l,k2*rp) * integral3;
    };
    int zstep = this->zsteps_Ml;
    double integral2 = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml, zstep);
    double res1 = pow(this->prefactor_Ml,2) * r*r / (Hf_interp(z)*1000.0) * this->Tb_interp(z, Tb_index) *\
                  this->bessel_j_interp_cubic(l,k1*r) * integral2;
    const double n_old = (z - this->zmin_Ml)/this->stepsize_Ml;
    int n;
    int n_old_int = (int)n_old;
    if (abs(n_old - (double)n_old_int) > 0.5)
        n = n_old_int + 1;
    else
        n = n_old_int;
    double qp, rr; 
    qp = this->q_p_Ml[n];
    rr = r*r;
    double hh = pow(Hf_interp(z)*1000.0, 2);
    double A = rr * this->Pk_interp(((double)l + 0.5)/q * qs[q_index].h,z, Pk_index)/(pow(qs[q_index].h,3)*hh*qp) *\
               pow(this->Tb_interp(z, Tb_index),2);

    double pre = 2*this->b_bias*this->b_bias*this->c*this->c/this->pi;
    double res2 = pre * A * rr / (q*q) * this->sph_bessel_camb(l,k1*r) * this->sph_bessel_camb(l, k2*r);

    cout << res1 << endl;
    cout << res2 << endl;
    cout << res2 / (14.6*res1) << endl;
}

double CosmoCalc::M(int l, double k1, double kappa, int Pk_index, int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r,q;
        r = r_interp(z);
        q = q_interp(z, q_index);


        //TODO: check whether we need to multiply py h.
        //return pow(r,2) * this->Tb_interp_full(z, Tb_index) * this->bessel_j_interp_cubic(l,k1*r) *\
        //   this->bessel_j_interp_cubic(l,k2*q) * sqrt(this->Pk_interp_full(k2*qs[q_index].h,z,Pk_index)/\
        //            pow(qs[q_index].h,3)) / (Hf_interp(z)*1000.0);

        return pow(r,2) * this->Tb_interp(z, Tb_index) * this->sph_bessel_camb(l,k1*r) *\
            this->sph_bessel_camb(l,kappa*q) * sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/\
                pow(qs[q_index].h,3)) / (Hf_interp(z)*1000.0);
    };

    //double integral = integrate(integrand, this->zmin_Ml, this->zmax_Ml,
    //this->zsteps_Ml, simpson());
    //
    //TODO: The problem here is that q_ml etc are not calculated necessarily for those values which
    //is why it all goes to shit.
    //int zstep = give_optimal_zstep(k1,k2);
    int zstep = this->zsteps_Ml;
    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml, zstep);

    return this->prefactor_Ml * integral;
}

double CosmoCalc::N_bar(int l, double k1, double k2, int Pk_index, int Tb_index, int q_index)
{
    auto integrand = [&](double z)
    {
        double r,q;
        r = r_interp(z);
        q = q_interp(z, q_index);


        double pref = r / (Hf_interp(z)*1000.0*(1+z));
        double pk = sqrt(this->Pk_interp(k2*qs[q_index].h, z, Pk_index)/pow(qs[q_index].h, 3));
        double dtb = this->Tb_interp(z, Tb_index);
        double pkdtb = pk * dtb;
        //double jl1r = this->bessel_j_interp_cubic(l - 1, k1 * r);
        //double jl2r = this->bessel_j_interp_cubic(l, k1 * r);
        //double jl1q = this->bessel_j_interp_cubic(l - 1, k2 * q);
        //double jl2q = this->bessel_j_interp_cubic(l, k2 * q);

        double jl1r = this->sph_bessel_camb(l - 1, k1 * r);
        double jl2r = this->sph_bessel_camb(l, k1 * r);
        double jl1q = this->sph_bessel_camb(l - 1, k2 * q);
        double jl2q = this->sph_bessel_camb(l, k2 * q);

        double sums = k1 * r * jl1r * jl1q -\
                      k1 * r * ((double)l+1.0) / (k2 * q) * jl1r * jl2q -\
                      ((double)l+1.0) * jl2r * jl1q +\
                      pow((double)l+1.0,2) / (k2 * q) * jl2r * jl2q;
        return pref * pkdtb * sums;

    };

    //double integral = integrate(integrand, this->zmin_Ml, this->zmax_Ml,
    //this->zsteps_Ml, simpson());

    int zstep = this->zsteps_Ml;//10000;//give_optimal_zstep(k1,k2);
    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml, zstep);

    return integral * this->prefactor_Ml;
}

void CosmoCalc::update_Tb_analytic(map<string, double> params, int *Tb_index)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Tbas.size(); ++i) {
        if (params["ombh2"] == Tbas[i].ombh2 && params["omnuh2"] == Tbas[i].omnuh2 &&\
                params["omch2"] == Tbas[i].omch2 && params["hubble"] == Tbas[i].hubble &&\
                params["T_CMB"] == Tbas[i].t_cmb) {

            do_calc = false;
            *Tb_index = i;
            break;
        }
    }
    if (do_calc) {
        Tb_analytic_interpolator interp;
        interp.ombh2 = params["ombh2"];
        interp.omnuh2 = params["omnuh2"];
        interp.omch2 = params["omch2"];
        interp.hubble = params["hubble"];
        interp.t_cmb = params["T_CMB"];

        cout << "Tb_analytic is being updated" << endl;

        double T_CMB2 = params["T_CMB"];
        double H_02 = params["hubble"];
        double h2 = H_02 / 100.0;
        double O_b2 = params["ombh2"] / pow(h2,2);
        double O_cdm2 = params["omch2"] / pow(h2,2);
        double O_nu2 = params["omnuh2"] / pow(h2,2);
        double O_M2 = O_b2 + O_cdm2 + O_nu2;

        real_1d_array xs, ys, hs;
        xs.setlength(this->zsteps_Ml+1);
        ys.setlength(this->zsteps_Ml+1);
        double z;
        for (int n = 0; n <= this->zsteps_Ml; ++n) {
            z = this->zmin_Ml + n * this->stepsize_Ml;
            xs[n] = z;

            double constant_A = 27.0 * O_b2 * pow(h2, 2) / 0.023 *\
                                sqrt(0.015 / (O_M2 * pow(h2,2)));
            double T_S = this->T_S(z);
            double Tz = T_CMB2 * (1+z);
            ys[n] = constant_A * this->x_HI(z) * (T_S - Tz)/T_S * sqrt(1+z);

        }
        spline1dinterpolant interpolator;
        spline1dbuildlinear(xs,ys,interpolator);

        interp.interpolator = interpolator;

        Tbas.push_back(interp);
        *Tb_index = Tbas.size() - 1;

        cout << "Tb analytic update done" << endl;
    }
}

double CosmoCalc::Tb_analytic_interp(double z, int Tb_index)
{
    return spline1dcalc(Tbas[Tb_index].interpolator,z);
}

double CosmoCalc::delta_Tb_bar(double z)
{
    double constant_A = 27.0 * this->O_b * pow(this->h, 2) / 0.023 *\
                        sqrt(0.015 / (this->O_M * pow(this->h,2)));
    double T_S = this->T_S(z);
    return constant_A * this->x_HI(z) * (T_S - this->T(z))/T_S * sqrt(1+z);
}
double CosmoCalc::delta_Tb_bar_G21(double z)
{
    return G21->getTb_interp_cubic(z);
}
double CosmoCalc::T_S(double z)
{
    double Ti = 10.0;
    double Tf = 500.0;
    double rate = 2.0;

    double zi = this->z_rei;
    double zf = zi - this->delta_z_rei;
    double deltaT = abs(Ti -Tf);

    return deltaT * (1.0/this->pi * atan(rate * (-(z - (zi + zf)/2.0))) +\
            0.5) + Ti;
}

double CosmoCalc::x_HI(double z)
{
    double rate = 2.0;
    double zi = this->z_rei;
    double zf = zi - this->delta_z_rei;

    return 1.0/this->pi * atan(rate * (z - (zi + zf)/2.0)) + 0.5;
}

double CosmoCalc::T_K(double z)
{
    double zd = 200.0;
    double res;
    if (z > zd)
        res = this->T(z);
    else {
        double Td = this->T(zd);
        double Tf = 1000.0;
        double rate = 1.0;
        double z_on = this->z_rei;
        double z0 = (2 * z_on - 4)/2.0;
        double tanh_term = (0.5 * (tanh(rate * (-(z-z0)))+1)) * Tf;
        res = Td * pow(1+z,2) / pow(1+zd,2) + tanh_term;
    }
    return res;
}

double CosmoCalc::integrandMM(int l, double k1, double k2, double k, int Pk_index, int Tb_index, int q_index)
{
    if (k1 == k2)
        return pow(k,2) * pow(this->M(l,k1,k,Pk_index,Tb_index,q_index),2);
    else 
        return pow(k,2) * this->M(l,k1,k,Pk_index,Tb_index,q_index) * this->M(l,k2,k,Pk_index,Tb_index,q_index);
}
double CosmoCalc::integrandMN(int l, double k1, double k2, double k, int Pk_index, int Tb_index, int q_index)
{
    return k * this->M(l,k1,k,Pk_index,Tb_index,q_index) * this->N_bar(l,k2,k,Pk_index,Tb_index,q_index);
}

double CosmoCalc::integrandNN(int l, double k1, double k2, double k, int Pk_index, int Tb_index, int q_index)
{
    if (k1 == k2)
        return pow(this->N_bar(l,k1,k,Pk_index,Tb_index,q_index),2);
    else 
        return this->N_bar(l,k1,k,Pk_index,Tb_index,q_index) * this->N_bar(l,k2,k,Pk_index,Tb_index,q_index);
}

double CosmoCalc::integrandsimple(int l, double k1, double k2, double z, int Pk_index, int Tb_index, int q_index)
{
    const double n_old = (z - this->zmin_Ml)/this->stepsize_Ml;
    int n;
    int n_old_int = (int)n_old;
    if (abs(n_old - (double)n_old_int) > 0.5)
        n = n_old_int + 1;
    else
        n = n_old_int;
    double r,q,qp;
    r = r_interp(z);
    q = q_interp(z, q_index);

    qp = this->q_p_Ml[n];
    double hh = pow(Hf_interp(z)*1000.0, 2);

    return pow(this->prefactor_Ml,2)* pow(r,4) / abs(qp) * pow(this->Tb_interp(z, Tb_index),2) *\
        this->pi / (2*pow(q,2)) *\
        this->bessel_j_interp_cubic(l,k1*r) *\
        this->bessel_j_interp_cubic(l,k2*r) *\
        this->Pk_interp(((double)l + 0.5)/q * qs[q_index].h,z,Pk_index)/(pow(qs[q_index].h,3)*hh);
}

double CosmoCalc::help_long(int l, double kp, double kappa, int Pk_index, int Tb_index, int q_index)
{
    auto integrand1 = [&](double zz)
    {

        double r,q,hh,Tb;
        r = r_interp(zz);
        q = q_interp(zz, q_index);

        hh = Hf_interp(zz)*1000.0;
        Tb = this->Tb_interp(zz, Tb_index);

        return pow(r,2)/hh * Tb * sqrt(this->Pk_interp(kappa*qs[q_index].h,zz, Pk_index)/pow(qs[q_index].h,3)) *\
            this->bessel_j_interp_cubic(l,kp*r) * this->bessel_j_interp_cubic(l,kappa*q);
    };

    double integral = integrate_simps(integrand1, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return integral;
}
double CosmoCalc::integrandlong(int l, double k1, double k2, double z, int Pk_index, int Tb_index, int q_index)
{
    double r,hh, Tb;
    r = r_interp(z);
    hh = Hf_interp(z)*1000.0;
    Tb = this->Tb_interp(z, Tb_index);

    int steps = (int)((1 - 0.001)/0.0001);
    if (steps % 2 == 1)
        ++steps;


    double res = pow(this->prefactor_Ml,2) * pow(r,2)/hh * Tb * this->bessel_j_interp_cubic(l,k1*r);


    auto integrand1 = [&](double kappa)
    {
        double q;
        q = q_interp(z, q_index);


        return pow(kappa,2) * sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/pow(qs[q_index].h,3)) *\
            bessel_j_interp_cubic(l, kappa*q) * help_long(l,k2,kappa, Pk_index, Tb_index, q_index);
    };

    double integral = integrate_simps(integrand1, 0.001, 1, steps);
    return res * integral;
}

double CosmoCalc::limber(int l, double r)
{
    auto integrand2 = [&](double q)
    {
        auto integrand1 = [&](double k)
        {
            //return boost::math::sph_bessel(l,k*r)*boost::math::sph_bessel(l,k*q);
            return sph_bessel_camb(l, k*r) * sph_bessel_camb(l, k*q);
        };

        //return integrate_simps(integrand1,0,1000,2000);
        return integrate_levin(integrand1,0,100);
    };

    //double res = integrate_simps(integrand2,0.9,1.1, 1000);
    double res = integrate_levin(integrand2,0.9,1.1);
    double nu = pow(l+0.5,2); 
    return res * 2 /this->pi * nu;
}

double CosmoCalc::limber2(int l, double z)
{
    auto integrand2 = [&](double zp)
    {
        auto integrand1 = [&](double k)
        {
            return pow(k,3) * sph_bessel_camb(l, k/z) * sph_bessel_camb(l, k/zp);
        };

        return integrate_simps(integrand1,10000,30000,30000);
    };

    double range = 0.2;
    double res = integrate_simps(integrand2,z-range,z+range, 1000);
    //double nu = pow(l+0.5,2); 
    return res * 2 /this->pi / (l+0.5);
}

double CosmoCalc::integrand_kappa_norsd(int l, double k1, double k2, double kappa, int Pk_index,\
        int Tb_index, int q_index)
{
    if (k1 == k2)
        return pow(kappa,2) * pow(this->M(l,k1,kappa,Pk_index,Tb_index,q_index),2);
    else 
        return pow(kappa,2) * this->M(l,k1,kappa,Pk_index,Tb_index,q_index) *\
            this->M(l,k2,kappa,Pk_index,Tb_index,q_index);

}

double CosmoCalc::integrand_kappa_rsd(int l, double k1, double k2, double kappa, int Pk_index,\
        int Tb_index, int q_index)
{
    double m1,n1,m2,n2;
    m1 = this->M(l,k1,kappa,Pk_index,Tb_index,q_index);
    n1 = this->N_bar(l,k1,kappa,Pk_index,Tb_index,q_index);
    if (k1 == k2) {
        m2 = m1;
        n2 = n1;
    } else {
        m2 = this->M(l,k2,kappa,Pk_index,Tb_index,q_index);
        n2 = this->N_bar(l,k2,kappa,Pk_index,Tb_index,q_index);
    }
    const double bb = this->b_bias * this->beta;
    const double bb2 = pow(bb,2);

    return pow(kappa,2) * m1 * m2 + bb * kappa * (m1*n2 + n1*m2) + bb2 * n1 * n2;
}
