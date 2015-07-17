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
    this->k_stepsize = this->fiducial_params["k_stepsize"];

    //generate object that is the CAMB interface.
    CAMB = new CAMB_CALLER;

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
    //updateClass(this->fiducial_params);
    cout << "... Class initialized ..." << endl;

    cout << "... precalculating Ml dependencies ..." << endl;
    this->update_q();
    this->update_q_prime();
    this->r_Ml = this->q_Ml;

    double z;
    for (int i = 0; i <= this->zsteps_Ml; ++i) {
        z = this->zmin_Ml + i * this->stepsize_Ml;
        this->H_f.push_back(this->H(z));
    }
    this->prefactor_Ml = 2*this->b_bias * this->c / this->pi;
    cout << "... Dependencies calculated ..." << endl;

    cout << "... Initializing Pk interpolator ..." << endl;
    //this->update_Pk_interpolator(this->fiducial_params);
    this->update_Pk_interpolator_direct(this->fiducial_params);
    cout << "... Pks calculated ..." << endl;

    cout << "... Creating Bessels ..." << endl;
    //this->create_bessel_interpolant_ALGLIB(0, this->fiducial_params["l_max"]);
    this->create_bessel_interpolant_OWN(this->fiducial_params["l_min"],this->fiducial_params["l_max"]);
    cout << "... Bessels built ..." << endl;

    cout << "... generating 21cm interface ..." << endl;
    G21 = new Global21cmInterface();
    this->update_G21(fiducial_params);
    cout << "... 21cm interface built ..." << endl;

    cout << "... CosmoCalc built ..." << endl;

}

double CosmoCalc::Cl(int l, double k1, double k2, double k_low, double k_high)
{
    //double lambda = this->current_params["ombh2"];
    //cout << lambda << endl;
    //return pow(lambda, l) * (k1+k2);
    //return this->corr_Tb(l, k1, k2, k_low, k_high);
    //return this->corr_Tb_rsd(l, k1, k2, k_low, k_high);
    //return this->Cl_simplified(l, k1, k2);
    //return this->Cl_simplified_rsd(l,k1,k2);
    return this->Cl_simplified(l,k1,k2) + this->Cl_noise(l,k1,k2);
}

double CosmoCalc::Cl_noise(int l, double k1, double k2)
{
    //TODO: integrand needs to be corrected.
    auto integrand = [&](double z)
    {
        const double n_old = (z - this->zmin_Ml)/this->stepsize_Ml;
        int n;
        int n_old_int = (int)n_old;
        if (abs(n_old - (double)n_old_int) > 0.5)
            n = n_old_int + 1;
        else
            n = n_old_int;
        double r;
        r = this->r_Ml[n];
        double jl = sph_bessel_camb(l,k1*r);
        double hub = this->H_f[n]*1000.0;
        return r*r*jl/hub; 
    };
    
    if (k1==k2) {
        // in mK
        double Tsys = 700000;
        double fcover = 1.0;
        double lmax = 5000;
        double tau = 365.25*24*60*60;
        double prefactor = 2.0 *pi*c*c * Tsys*Tsys/(fcover*fcover * fiducial_params["df"] * lmax * lmax * tau);
        double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
            this->zsteps_Ml);
        return prefactor * integral * integral;
    } else {
        return 0.0;
    }
}

CosmoCalc::~CosmoCalc()
{
    // delete CLASS;
    delete CAMB;
    delete G21;
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

void CosmoCalc::update_q_prime()
{
    this->q_p_Ml.clear();
    double z;
    double h = 10e-4;
    for (int n = 0; n <= this->zsteps_Ml; ++n) {
        z = this->zmin_Ml + n * this->stepsize_Ml;
        double res = 0;
        res = - D_C(z+2*h) + 8 * D_C(z+h) - 8 * D_C(z-h) + D_C(z-2*h);
        res = abs(res);

        this->q_p_Ml.push_back(res/(12*h));
    }
}

void CosmoCalc::create_bessel_interpolant_ALGLIB(int lmin, int lmax)
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
}

void CosmoCalc::create_bessel_interpolant_OWN(int lmin, int lmax)
{
    double xmax = 2 * 5 * this->r_Ml[this->zsteps_Ml];
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

void CosmoCalc::update_G21(map<string,double> params)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Tbs.size(); ++i) {
        if (params["ombh2"] == Tbs[i].ombh2 && params["omnuh2"] == Tbs[i].omnuh2 &&\
                params["omch2"] == Tbs[i].omch2 && params["omk"] == Tbs[i].omk &&\
                params["hubble"] == Tbs[i].hubble && params["sigma8"] == Tbs[i].s8 &&\
                params["T_CMB"] == Tbs[i].T_CMB && params["n_s"] == Tbs[i].n_s &&\
                params["fstar"] == Tbs[i].fstar && params["fesc"] == Tbs[i].fesc &&\
                params["nion"] == Tbs[i].nion && params["fx"] == Tbs[i].fx &&\
                params["flya"] == Tbs[i].flya) {

            do_calc = false;
            this->Tb_index = i;
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
        this->Tb_index = Tbs.size() - 1;

    }
}

double CosmoCalc::Tb_interp(double z)
{
    // The * 1000.0 is so we get the result in mK
    return spline1dcalc(Tbs[Tb_index].interpolator,z) * 1000.0;
}

void CosmoCalc::update_Pk_interpolator_direct(map<string, double> params)
{
    bool do_calc = true;
    for (unsigned int i = 0; i < Pks.size(); ++i) {
        if (params["ombh2"] == Pks[i].ombh2 && params["omnuh2"] == Pks[i].omnuh2 &&\
                params["omch2"] == Pks[i].omch2 && params["omk"] == Pks[i].omk &&\
                params["hubble"] == Pks[i].hubble) {

            do_calc = false;
            this->Pk_index = i;
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
        this->Pk_index = Pks.size() - 1;
    }
}
double CosmoCalc::Pk_interp(double k, double z)
{
    return spline2dcalc(Pks[this->Pk_index].interpolator, k, z);
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



double CosmoCalc::Cl_simplified_rsd(int l, double k1, double k2)
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
        double r,rr,q,qq,qp,k1r,k2r;
        r = this->r_Ml[n];
        rr = r*r;
        q = this->q_Ml[n];
        qq = q*q;
        qp = this->q_p_Ml[n];
        k1r = k1 * r;
        k2r = k2 * r;
        double hh = pow(this->H_f[n]*1000.0, 2);
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
        double LL2 = 2.0*(-32.0*pow(l,6) - 24.0*pow(l,5) + 48.0*pow(l,4) + 46.0*pow(l,3) - 5.0*l - 1.0);
        LL2 = LL2 / (double)(pow(l,3)*pow(2*l-1,2)*pow(2*l+1,4));
        //LL2 = -1.0/(double)(300*300*300);
        L4.push_back(LL2 * pow(l+1,2));
        L4.push_back(-LL2 * k2r * (l+1));
        L4.push_back(-LL2 * k1r * (l+1));       
        L4.push_back(LL2 * k1r * k2r);
        double A = rr * pow(this->Tb_interp(z),2) * this->Pk_interp((double)l/q * this->h,z)/\
                   (pow(this->h,3)*hh*qp);   
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
        double bracket = rr/qq * Jl1 + bb*r/((1+z)*q) * (JL2+JL3) + bb2/pow(1+z,2) * JL4; 

        return A * bracket;

    };
    double prefact = pow(this->prefactor_Ml,2) * this->pi / 2.0;
    return prefact * integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
}

double CosmoCalc::Cl_simplified(int l, double k1, double k2)
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
        double r,q,qp,rr;
        r = this->r_Ml[n];
        q = this->q_Ml[n];
        qp = this->q_p_Ml[n];
        rr = r*r;
        double hh = pow(this->H_f[n]*1000.0, 2);
        double A = rr * this->Pk_interp(((double)l + 0.5)/q * this->h,z)/(pow(this->h,3)*hh*qp) *\
                   pow(this->Tb_interp(z),2);

        //TODO: check whether we need to multiply py h.
        // here: changed it back to non - interpolation, because it isn't necessary for the simplified case.
        return A * rr / (q*q) * this->sph_bessel_camb(l,k1*r) * this->sph_bessel_camb(l, k2*r);
    };

    double pre = pow(this->prefactor_Ml,2) * this->pi / 2.0;
    return  pre * integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml,\
            this->zsteps_Ml);
}

double CosmoCalc::corr_Tb(int l, double k1, double k2, double k_low,\
        double k_high)
{
    int steps = (int)((k_high - k_low)/this->k_stepsize);

    if (k1 == k2)
    {
        auto integrand = [&](double k)
        {
            return pow(k,2) * pow(this->M(l,k1,k),2);
        };

        //return integrate(integrand, k_low, k_high, this->k_steps, simpson());
        return integrate_simps(integrand, k_low, k_high, steps);
    } else {
        auto integrand = [&](double k)
        {
            return pow(k,2) * this->M(l,k1,k) * this->M(l,k2,k);
        };

        //return integrate(integrand, k_low, k_high, this->k_steps, simpson());
        return integrate_simps(integrand, k_low, k_high, steps);
    }
}

double CosmoCalc::corr_Tb_rsd(int l, double k1, double k2, double k_low,\
        double k_high)
{       
    int steps = (int)((k_high - k_low)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;

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

    return integrate_simps(integrand, k_low, k_high, steps);
}

double CosmoCalc::corr_Tb_new(int l, double k1, double k2, double k_low,\
        double k_high)
{
    int steps = (int)((k_high - k_low)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;

    auto integrand1 = [&](double z)
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

        auto integrand2 = [&](double zp)
        {
            const double n_old2 = (zp - this->zmin_Ml)/this->stepsize_Ml;
            int n2;
            int n_old_int2 = (int)n_old2;
            if (abs(n_old2 - (double)n_old_int2) > 0.5)
                n2 = n_old_int2 + 1;
            else
                n2 = n_old_int2;
            double rp,qp;
            rp = this->r_Ml[n2];
            qp = this->q_Ml[n2];

            auto integrand3 = [&](double kappa)
            {
                double hhh = pow(this->h,3);
                double sP = sqrt(this->Pk_interp(kappa*this->h,z)/hhh);
                double sPp = sqrt(this->Pk_interp(kappa*this->h,zp)/hhh);
                return kappa*kappa * sP * sPp * this->bessel_j_interp_cubic(l,kappa*q) *\
                    this->bessel_j_interp_cubic(l,kappa*qp);
            };
            double integral3 = integrate_simps(integrand3, k_low, k_high, steps);
            return rp*rp / (this->H_f[n2]*1000.0) * this->Tb_interp(zp) *\
                this->bessel_j_interp_cubic(l,k2*rp) * integral3;
        };
        int zstep = 1000;
        double integral2 = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml, zstep);
        return r*r / (this->H_f[n]*1000.0) * this->Tb_interp(z) *\
            this->bessel_j_interp_cubic(l,k1*r) * integral2;
    };
    int zstep = 1000;
    double integral1 = integrate_simps(integrand1,this->zmin_Ml, this->zmax_Ml, zstep);
    return pow(this->prefactor_Ml,2) * integral1;
}

void CosmoCalc::compare(int l, double k1, double k2)
{
    double k_high = 0.6;
    double k_low = 0.4;
    int steps = (int)((k_high - k_low)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;
    double z = 8.0;
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

    auto integrand2 = [&](double zp)
    {
        const double n_old2 = (zp - this->zmin_Ml)/this->stepsize_Ml;
        int n2;
        int n_old_int2 = (int)n_old2;
        if (abs(n_old2 - (double)n_old_int2) > 0.5)
            n2 = n_old_int2 + 1;
        else
            n2 = n_old_int2;
        double rp,qp;
        rp = this->r_Ml[n2];
        qp = this->q_Ml[n2];

        auto integrand3 = [&](double kappa)
        {
            double hhh = pow(this->h,3);
            double sP = sqrt(this->Pk_interp(kappa*this->h,z)/hhh);
            double sPp = sqrt(this->Pk_interp(kappa*this->h,zp)/hhh);
            return kappa*kappa * sP * sPp * this->bessel_j_interp_cubic(l,kappa*q) *\
                this->bessel_j_interp_cubic(l,kappa*qp);
        };
        double integral3 = integrate_simps(integrand3, k_low, k_high, steps);
        return rp*rp / (this->H_f[n2]*1000.0) * this->Tb_interp(zp) *\
            this->bessel_j_interp_cubic(l,k2*rp) * integral3;
    };
    int zstep = 10000;
    double integral2 = integrate_simps(integrand2, 7.5, 8.5, zstep);
    double res1 = pow(this->prefactor_Ml,2) * r*r / (this->H_f[n]*1000.0) * this->Tb_interp(z) *\
                  this->bessel_j_interp_cubic(l,k1*r) * integral2;

    double qp, rr; 
    qp = this->q_p_Ml[n];
    rr = r*r;
    double hh = pow(this->H_f[n]*1000.0, 2);
    double A = rr * this->Pk_interp(((double)l + 0.5)/q * this->h,z)/(pow(this->h,3)*hh*qp) *\
               pow(this->Tb_interp(z),2);

    double pre = 2*this->b_bias*this->b_bias*this->c*this->c/this->pi;
    double res2 = pre * A * rr / (q*q) * this->sph_bessel_camb(l,k1*r) * this->sph_bessel_camb(l, k2*r);

    cout << res1 << endl;
    cout << res2 << endl;
    cout << res2 / res1 << endl;
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
        return pow(r,2) * this->Tb_interp(z) * this->bessel_j_interp_cubic(l,k1*r) *\
            this->bessel_j_interp_cubic(l,k2*q) * sqrt(this->Pk_interp(k2*this->h,z)/\
                    pow(this->h,3)) / (this->H_f[n]*1000.0);

        //return pow(r,2) * this->delta_Tb_bar(z) * this->sph_bessel_camb(l,k1*r) *
        //this->sph_bessel_camb(l,k2*q) * sqrt(this->Pk_interp(k2*this->h,z)/
        //        pow(this->h,3)) / (this->H_f[n]*1000.0);
    };

    //double integral = integrate(integrand, this->zmin_Ml, this->zmax_Ml,
    //this->zsteps_Ml, simpson());
    int zstep = give_optimal_zstep(k1,k2);
    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml, zstep);

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

        double pref = r / (this->H_f[n]*1000.0*(1+z));
        double pk = sqrt(this->Pk_interp(k2*this->h, z)/pow(this->h, 3));
        double dtb = this->Tb_interp(z);
        double pkdtb = pk * dtb;
        double jl1r = this->bessel_j_interp_cubic(l - 1, k1 * r);
        double jl2r = this->bessel_j_interp_cubic(l, k1 * r);
        double jl1q = this->bessel_j_interp_cubic(l - 1, k2 * q);
        double jl2q = this->bessel_j_interp_cubic(l, k2 * q);

        //double jl1r = this->sph_bessel_camb(l - 1, k1 * r);
        //double jl2r = this->sph_bessel_camb(l, k1 * r);
        //double jl1q = this->sph_bessel_camb(l - 1, k2 * q);
        //double jl2q = this->sph_bessel_camb(l, k2 * q);

        double sums = k1 * r * jl1r * jl1q -\
                      k1 * r * ((double)l+1.0) / (k2 * q) * jl1r * jl2q -\
                      ((double)l+1.0) * jl2r * jl1q +\
                      pow((double)l+1.0,2) / (k2 * q) * jl2r * jl2q;
        return pref * pkdtb * sums;

    };

    //double integral = integrate(integrand, this->zmin_Ml, this->zmax_Ml,
    //this->zsteps_Ml, simpson());

    int zstep = 10000;//give_optimal_zstep(k1,k2);
    double integral = integrate_simps(integrand, this->zmin_Ml, this->zmax_Ml, zstep);

    return integral * this->prefactor_Ml;
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

double CosmoCalc::integrandMM(int l, double k1, double k2, double k)
{
    if (k1 == k2)
        return pow(k,2) * pow(this->M(l,k1,k),2);
    else 
        return pow(k,2) * this->M(l,k1,k) * this->M(l,k2,k);
}
double CosmoCalc::integrandMN(int l, double k1, double k2, double k)
{
    return k * this->M(l,k1,k) * this->N_bar(l,k2,k);
}

double CosmoCalc::integrandNN(int l, double k1, double k2, double k)
{
    if (k1 == k2)
        return pow(this->N_bar(l,k1,k),2);
    else 
        return this->N_bar(l,k1,k) * this->N_bar(l,k2,k);
}

double CosmoCalc::integrandsimple(int l, double k1, double k2, double z)
{
    const double n_old = (z - this->zmin_Ml)/this->stepsize_Ml;
    int n;
    int n_old_int = (int)n_old;
    if (abs(n_old - (double)n_old_int) > 0.5)
        n = n_old_int + 1;
    else
        n = n_old_int;
    double r,q,qp;
    r = this->r_Ml[n];
    q = this->q_Ml[n];
    qp = this->q_p_Ml[n];
    double hh = pow(this->H_f[n]*1000.0, 2);

    return pow(this->prefactor_Ml,2)* pow(r,4) / abs(qp) * pow(this->Tb_interp(z),2) *\
        this->pi / (2*pow(q,2)) *\
        this->bessel_j_interp_cubic(l,k1*r) *\
        this->bessel_j_interp_cubic(l,k2*r) *\
        this->Pk_interp(((double)l + 0.5)/q * this->h,z)/(pow(this->h,3)*hh);
}




double CosmoCalc::help_long(int l, double kp, double kappa)
{
    auto integrand1 = [&](double zz)
    {
        const double n_old = (zz - this->zmin_Ml)/this->stepsize_Ml;
        int n;
        int n_old_int = (int)n_old;
        if (abs(n_old - (double)n_old_int) > 0.5)
            n = n_old_int + 1;
        else
            n = n_old_int;
        double r,q,hh,Tb;
        r = this->r_Ml[n];
        q = this->q_Ml[n];
        hh = this->H_f[n]*1000.0;
        Tb = this->Tb_interp(zz);

        return pow(r,2)/hh * Tb * sqrt(this->Pk_interp(kappa*this->h,zz)/pow(this->h,3)) *\
            this->bessel_j_interp_cubic(l,kp*r) * this->bessel_j_interp_cubic(l,kappa*q);
    };

    double integral = integrate_simps(integrand1, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return integral;
}
double CosmoCalc::integrandlong(int l, double k1, double k2, double z)
{
    const double n_old = (z - this->zmin_Ml)/this->stepsize_Ml;
    int n;
    int n_old_int = (int)n_old;
    if (abs(n_old - (double)n_old_int) > 0.5)
        n = n_old_int + 1;
    else
        n = n_old_int;
    double r,hh, Tb;
    r = this->r_Ml[n];
    hh = this->H_f[n]*1000.0;
    Tb = this->Tb_interp(z);

    int steps = (int)((1 - 0.001)/0.0001);
    if (steps % 2 == 1)
        ++steps;


    double res = pow(this->prefactor_Ml,2) * pow(r,2)/hh * Tb * this->bessel_j_interp_cubic(l,k1*r);


    auto integrand1 = [&](double kappa)
    {
        const double n_old = (z - this->zmin_Ml)/this->stepsize_Ml;
        int n;
        int n_old_int = (int)n_old;
        if (abs(n_old - (double)n_old_int) > 0.5)
            n = n_old_int + 1;
        else
            n = n_old_int;
        double q;
        q = this->q_Ml[n];

        return pow(kappa,2) * sqrt(this->Pk_interp(kappa*this->h,z)/pow(this->h,3)) *\
            bessel_j_interp_cubic(l, kappa*q) * help_long(l,k2,kappa);
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
            return sph_bessel_camb(l, k*r) * sph_bessel_camb(l, k*q);
        };

        return integrate_simps(integrand1,0,10000,20000);
    };

    double res = integrate_simps(integrand2,0.9,1.1, 1000);
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
