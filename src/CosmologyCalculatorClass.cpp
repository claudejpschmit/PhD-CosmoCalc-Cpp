#include "CosmologyCalculatorClass.hpp"
#include <string>
#include <iostream>
#include "Integrator.hpp"
#include <fstream>

CosmoCalc::CosmoCalc(map<string, double> params)
    :
        CosmoBasis(params)
{
    cout << "... Beginning to build CosmoCalc ..." << endl;
    this->prefactor_Ml = 2*this->b_bias*this->c/pi;
    this->zmin_Ml = this->fiducial_params["zmin"];
    this->zmax_Ml = this->fiducial_params["zmax"];
    this->zsteps_Ml = this->fiducial_params["zsteps"];
    this->stepsize_Ml = (this->zmax_Ml - this->zmin_Ml)/(double)this->zsteps_Ml;
    this->Pk_steps = this->fiducial_params["Pk_steps"];
    this->k_steps = this->fiducial_params["k_steps"];

    pars.add("100*theta_s",0);
    pars.add("omega_b",0);
    pars.add("omega_cdm",0);
    pars.add("A_s",0);
    pars.add("n_s",0);
    pars.add("tau_reio",0);
    pars.add("k_pivot",0);
    pars.add("YHe",0);
    pars.add("z_pk", 0);
    //This doesn't work currently
    //pars.add("h", 0);
    pars.add("Omega_k", 0);
    //pars.add("Omega_Lambda", 0);
    pars.add("T_cmb", 0);
    //pars.add("bias", 0);

    //Unchanging
    pars.add("output","mPk"); //pol +clphi
    pars.add("P_k_max_h/Mpc", 100);
    cout << "... Initializing Class ..." << endl;
    updateClass(this->fiducial_params);
    cout << "... Class initialized ..." << endl;
    
    cout << "... precalculating Ml dependencies ..." << endl;
    this->update_q();
    this->r_Ml = this->q_Ml;
    
    double z;
    for (int i = 0; i <= this->zsteps_Ml; ++i) {
        z = this->zmin_Ml + i * this->stepsize_Ml;
        this->H_f.push_back(this->H(z));
    }
    this->prefactor_Ml = 2*this->b_bias * this->c / this->pi;
    cout << "... Dependencies calculated ..." << endl;

    cout << "... Initializing Pk interpolator ..." << endl;
    this->update_Pk_interpolator(this->fiducial_params);
    cout << "... Pks calculated ..." << endl;

    cout << "... Creating Bessels ..." << endl;
    this->create_bessel_interpolant(0, this->fiducial_params["l_max"]);
    cout << "... Bessels built ..." << endl; 
    
    cout << "... CosmoCalc built ..." << endl;

}

CosmoCalc::~CosmoCalc()
{
    delete CLASS;
}

void CosmoCalc::write_pks(string filename, double z)
{
    ofstream file;
    file.open(filename);
    cout.precision( 16 );
    this->CLASS->writePks(file, z);
    file.close();

    filename = "new_"+filename;
    file.open(filename);
    
    double k = 0.0001;
    double kstep = k;
    for (int n = 0; n < 10000; ++n) {
        k += kstep;
        file << k << " " << this->Pk_interp(k,z) << endl; 
    }
    file.close();

}
void CosmoCalc::show_cosmo_calcs()
{
    write_pks("outputtest_z0.dat", 0);
    write_pks("outputtest_z1.dat", 1);
    write_pks("outputtest_z2.dat", 2);
    write_pks("outputtest_z755.dat", 7.55);
    //updateClass(current_params);
    current_params["ombh2"] = 0.04;
    updateClass(current_params);
    write_pks("outputtestuuu.dat", 0);
    
    cout << hubble_time() << endl;
    cout << hubble_dist() << endl;
    cout << comoving_radial_dist(10) << endl;
    cout << "O_m = " << O_M << ", O_V = " << O_V << "." << endl;
    cout << "Age in Gigayears "<< age_of_universe(0) * pow(10,10) *\
        3.08568 / (365.25 * 24 * 3600) << endl;

}

void CosmoCalc::updateClass(map<string, double> params)
{
    pars.updateParam("omega_b", params["ombh2"]);
    pars.updateParam("omega_cdm", params["omch2"]);
    //This doesn't work currently
    //pars.updateParam("h", params["hubble"]/100.0);
    pars.updateParam("Omega_k", params["omk"]);
    //pars.updateParam("Omega_Lambda", this->O_V);
    pars.updateParam("T_cmb", params["T_CMB"]);
    pars.updateParam("A_s", params["A_s"]);
    pars.updateParam("n_s", params["n_s"]);
    pars.updateParam("tau_reio", params["tau_reio"]);
    pars.updateParam("k_pivot", params["k_pivot"]);
    pars.updateParam("YHe", params["YHe"]);
    pars.updateParam("z_pk", params["z_pk"]);
    pars.updateParam("100*theta_s", params["100*theta_s"]);
    
    //Not quite sure if this is = b_bias...
    //pars.updateParam("bias", 1);
    
    CLASS = new ClassEngine(pars);
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

void CosmoCalc::update_q()
{
    this->q_Ml.clear();
    double z;
    for (int n = 0; n <= this->zsteps_Ml; ++n) {
        z = this->zmin_Ml + n * this->stepsize_Ml;
        this->q_Ml.push_back(this->D_C(z));
    }
}

void CosmoCalc::create_bessel_interpolant(int lmin, int lmax)
{
    double xmin, xmax;
    this->lmin_bess = lmin;
    // TODO: these values are kind of arbitrary...
    // a little under (0.5)
    xmin = 0.5 * 0.001 * this->r_Ml[0];
    // a little over (2)
    xmax = 2 * 5 * this->r_Ml[this->zsteps_Ml];
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

    //Trying to interpolate myself.
    for (int l = 0; l <= lmax; ++l) {
        vector<double> row;
        for (int j = 0; j < (int)xmax; ++j) {
            row.push_back(this->sph_bessel_camb(l,(double)j));
        }
        bessel_values.push_back(row);  
    }

}

double CosmoCalc::bessel_j_interp(int l, double x)
{
    return spline1dcalc(bessel_interp_list[l - this->lmin_bess], x);
}

double CosmoCalc::bessel_j_interp_basic(int l, double x)
{

    if ((int)x + 1 >= bessel_values[l].size()) {
        return bessel_values[l][bessel_values[l].size()-1];
    } else {
        if (x-(int)x != 0) {       
            return bessel_values[l][(int)x] +\
                   (bessel_values[l][(int)x + 1] -\
                   bessel_values[l][(int)x]) * (x - (int)x);
        } else {
            return bessel_values[l][(int)x];
        }
    }
}

void CosmoCalc::update_Pk_interpolator(map<string, double> params)
{
    
    double z_stepsize = (params["zmax"] - params["zmin"])/params["Pk_steps"];
    vector<double> vk, vz, vP;

    for (int i = 0; i < (int)params["Pk_steps"]; ++i) {
        vz.push_back(params["z_pk"] + i * z_stepsize);
        stringstream command;
        command << "python Pk.py";
        command << " --H_0 " << params["hubble"];
        command << " --ombh2 " << params["ombh2"];
        command << " --omnuh2 " << params["omnuh2"];
        command << " --omch2 " << params["omch2"];
        command << " --omk " << params["omk"];
        command << " --z " << vz[i];
        system(command.str().c_str());
        ifstream file("Pks.dat");
        
        double k, P;
        vk.clear();

        while (file >> k >> P) {
            vk.push_back(k);
            vP.push_back(P);
        }
        file.close();
    }

    matterpowerspectrum_k.setlength(vk.size());
    matterpowerspectrum_z.setlength(vz.size());
    matterpowerspectrum_P.setlength(vP.size());
    for (int i = 0; i < vk.size(); i++){
        matterpowerspectrum_k[i] = vk[i];
    }
    for (int i = 0; i < vP.size(); i++){
        matterpowerspectrum_P[i] = vP[i];
    }
    for (int i = 0; i < vz.size(); i++){
        matterpowerspectrum_z[i] = vz[i];
    }


    spline2dbuildbilinearv(matterpowerspectrum_k, vk.size(),matterpowerspectrum_z, vz.size(),\
                        matterpowerspectrum_P, 1, Pk_interpolator); 

    
    
    // that should create a file containing the Pks.
/*        
        readin(pks.dat)
        interpolate(pks)

        update interpolation function.

     */
}

double CosmoCalc::Pk_interp(double k, double z)
{
    return spline2dcalc(Pk_interpolator, k, z);
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

double CosmoCalc::Cl(int l, double k1, double k2, double k_low, double k_high)
{
    //double lambda = this->current_params["ombh2"];
    //cout << lambda << endl;
    //return pow(lambda, l) * (k1+k2);
    return this->corr_Tb_rsd(l, k1, k2, k_low, k_high);
}


double CosmoCalc::corr_Tb(int l, double k1, double k2, double k_low,\
                          double k_high)
{
    if (k1 == k2)
    {
        auto integrand = [&](double k)
        {
            return pow(k,2) * pow(this->M(l,k1,k),2);
        };

        //return integrate(integrand, k_low, k_high, this->k_steps, simpson());
        return integrate_simps(integrand, k_low, k_high, this->k_steps);
    } else {
        auto integrand = [&](double k)
        {
            return pow(k,2) * this->M(l,k1,k) * this->M(l,k2,k);
        };

        //return integrate(integrand, k_low, k_high, this->k_steps, simpson());
        return integrate_simps(integrand, k_low, k_high, this->k_steps);
    }
}

double CosmoCalc::corr_Tb_rsd(int l, double k1, double k2, double k_low,\
                              double k_high)
{
    auto integrand = [&](double k)
    {
        double m1,n1,m2,n2;
        m1 = this->M(l,k1,k);
        n1 = this->N_bar(l,k1,k);
        if (k1 == k2) {
            m2 = m1;
            n2 = n1;
        } else {
            m2 = this->M(l,k2,k);
            n2 = this->N_bar(l,k2,k);
        }
        const double bb = this->b_bias * this->beta;
        const double bb2 = pow(bb,2);
        
        return pow(k,2) * m1 * m2 + bb * k * (m1*n2 + n1*m2) + bb2 * n1 * n2;
    };
    //return integrate(integrand, k_low, k_high, this->k_steps, simpson());

    return integrate_simps(integrand, k_low, k_high, this->k_steps);
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
        r = this->r_Ml[n];
        q = this->q_Ml[n];
       
        //TODO: check whether we need to multiply py h.
        return pow(r,2) * this->delta_Tb_bar(z) * this->bessel_j_interp(l,k1*r) *\
                this->bessel_j_interp(l,k2*q) * sqrt(this->Pk_interp(k2*this->h,z)/\
                pow(this->h,3)) / (this->H_f[n]*1000.0);

        //return pow(r,2) * this->delta_Tb_bar(z) * this->sph_bessel_camb(l,k1*r) *\
                this->sph_bessel_camb(l,k2*q) * sqrt(this->Pk_interp(k2*this->h,z)/\
                pow(this->h,3)) / (this->H_f[n]*1000.0);
    };
    
    //double integral = integrate(integrand, this->zmin_Ml, this->zmax_Ml,\
                                this->zsteps_Ml, simpson());

    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
                                this->zsteps_Ml);
    return this->prefactor_Ml * integral;
}

double CosmoCalc::N_bar(int l, double k1, double k2)
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
        r = this->r_Ml[n];
        q = this->q_Ml[n];
        
        double pref = 1.0 / (this->H_f[n]*1000.0*(1+z)) * this->prefactor_Ml;
        double pk = sqrt(this->Pk_interp(k2*this->h, z)/pow(this->h, 3));
        double dtb = this->delta_Tb_bar(z);
        double pkdtb = pk * dtb;
        double jl1r = this->bessel_j_interp(l - 1, k1 * r);
        double jl2r = this->bessel_j_interp(l, k1 * r);
        double jl1q = this->bessel_j_interp(l - 1, k2 * q);
        double jl2q = this->bessel_j_interp(l, k2 * q);
/*
        double jl1r = this->sph_bessel_camb(l - 1, k1 * r);
        double jl2r = this->sph_bessel_camb(l, k1 * r);
        double jl1q = this->sph_bessel_camb(l - 1, k2 * q);
        double jl2q = this->sph_bessel_camb(l, k2 * q);
*/
        double sums = k1 * r * jl1r * jl1q -\
                      k1 * r * (l+1) / (k2 * q) * jl1r * jl2q -\
                      (l+1) * jl2r * jl1q +\
                      pow(l+1,2) / (k2 * q) * jl2r * jl2q;
        return pref * r * pkdtb * sums;
  
    };
    
    //double integral = integrate(integrand, this->zmin_Ml, this->zmax_Ml,\
                                this->zsteps_Ml, simpson());
    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
                                this->zsteps_Ml);

    return integral;
}

double CosmoCalc::delta_Tb_bar(double z)
{
    double constant_A = 27 * this->O_b * pow(this->h, 2) / 0.023 *\
                        sqrt(0.015 / (this->O_M * pow(this->h,2)));
    double T_S = this->T_S(z);
    return constant_A * this->x_HI(z) * (T_S - this->T(z))/T_S * sqrt(1+z);
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
