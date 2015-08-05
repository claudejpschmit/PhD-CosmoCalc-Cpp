#include "SanityChecker.hpp"
#include "Integrator.hpp"
#include <time.h>

SanityChecker::SanityChecker(map<string, double> params)
    :
        CosmoCalc(params)
{
    cout << "... Beginning to build SanityChecker ..." << endl;

    //this->create_bessel_interpolant_ALGLIB(this->fiducial_params["l_min"], this->fiducial_params["l_max"]);

    cout << "... SanityChecker built ..." << endl;
}

SanityChecker::~SanityChecker()
{}

void SanityChecker::kappa_integral(int l, double z, double zp, double *out1, double k_low, double k_high, int n)
{
    double q;
    q = spline1dcalc(q_interp, z);   

    double qp;
    qp = spline1dcalc(q_interp,zp);

    double hhh = pow(this->h,3);

    double a;
    //cout<< a << endl;
    double low;
    if (l < 1000){
        low = (double)l/(1.2*q);
        a = 2*(double)(l+1000)/q;
    } else {
        low = (double)l/(q);
        a = (double)(l+1000)/q;
    }
    double lower_kappa_bound;
    //cout << low << endl;
    if (low > k_low)
        lower_kappa_bound = low;
    else
        lower_kappa_bound = k_low;

    int steps = (int)((k_high - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;



    auto integrand3 = [&](double kap)
    {
        double sP1 = sqrt(this->Pk_interp(kap*this->h,z)/hhh);
        double sPp1 = sqrt(this->Pk_interp(kap*this->h,zp)/hhh);
        return kap*kap * sP1 * sPp1 * this->sph_bessel_camb(l,kap*q) *\
            this->sph_bessel_camb(l,kap*qp);
    };
    double integral2 = integrate_simps(integrand3, lower_kappa_bound, k_high, steps);

    LEVIN = new Levin(a, k_high);

    auto foo = [&](double kappa)
    {
        double sP1 = sqrt(this->Pk_interp(kappa*this->h,z)/hhh);
        double sPp1 = sqrt(this->Pk_interp(kappa*this->h,zp)/hhh);
        return kappa*kappa * sP1 * sPp1;     
    };

    steps = (int)((a - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;

    cout << lower_kappa_bound << " to " << a << " with " << steps << " steps"<< endl;
    //cout << steps << endl;
    double integral = integrate_simps(integrand3, lower_kappa_bound, a, steps);
    double integral3;
    if (z == zp)
        integral3 = LEVIN->integrate_2sphj_1r(foo,q,l,n);
    else
        integral3 = LEVIN->integrate_2sphj_2r(foo,q,qp,l,n);

    delete LEVIN;
    integral3 += integral;
    *out1 = integral2/integral3;
}

double SanityChecker::kappa_integrand(int l, double z, double zp, double kappa)
{
    double q = spline1dcalc(q_interp, z);   
    double qp = spline1dcalc(q_interp,zp);

    double hhh = pow(this->h,3);
    double sP1 = sqrt(this->Pk_interp(kappa*this->h,z)/hhh);
    double sPp1 = sqrt(this->Pk_interp(kappa*this->h,zp)/hhh);
    return kappa*kappa * sP1 * sPp1 * this->sph_bessel_camb(l,kappa*q) *\
        this->sph_bessel_camb(l,kappa*qp);
}

void SanityChecker::Compare_Cl(int l, double k1, double k2, double k_low, double k_high, int n_levin, double *ratio, double *time_r)
{
    clock_t t1,t2;
    t1 = clock();
    double Cl_levin = Cl_new(l,k1,k2,k_low,k_high,n_levin);
    t2 = clock();
    double ta = (double)t2 - (double)t1;
    cout << "C_levin = " << Cl_levin << endl;
    cout << "time = " << ta/CLOCKS_PER_SEC << endl;

    t1 = clock();
    double Cl_full = corr_Tb_new(l, k1, k2, k_low, k_high);
    t2 = clock();
    double tb = (double)t2 - (double)t1;
    cout << "C_full = " << Cl_full << endl;
    cout << "time = " << tb/CLOCKS_PER_SEC << endl;
    cout << "ratio = " << Cl_levin/Cl_full << endl;
    *ratio = Cl_levin/Cl_full;
    *time_r = tb/ta;
    cout << "time_ratio = " << tb/ta << endl;
}

double SanityChecker::Cl_new(int l, double k1, double k2, double k_low,\
        double k_high, int n_levin)
{
    double a;
    double low;
    double hhh = pow(this->h,3);

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

        if (l < 1000){
            low = (double)l/(1.2*q);
            a = (double)(l+1000)/q;
        } else {
            low = (double)l/(q);
            a = (double)(l+1000)/q;
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
                return kappa*kappa * sP * sPp * this->sph_bessel_camb(l,kappa*q) *\
                    this->sph_bessel_camb(l,kappa*qp);
            };
            double integral = integrate_simps(integrand3, lower_kappa_bound, a, steps);
            double integral3;
            
            LEVIN = new Levin(a, k_high);

            auto foo = [&](double kappa)
            {
                double sP1 = sqrt(this->Pk_interp(kappa*this->h,z)/hhh);
                double sPp1 = sqrt(this->Pk_interp(kappa*this->h,zp)/hhh);
                return kappa*kappa * sP1 * sPp1;     
            };

            
            if (z == zp)
                integral3 = LEVIN->integrate_2sphj_1r(foo,q,l,n_levin);
            else
                integral3 = LEVIN->integrate_2sphj_2r(foo,q,qp,l,n_levin);

            delete LEVIN;
            integral3 += integral;
            return rp*rp / (this->H_f[n2]*1000.0) * this->Tb_interp(zp) *\
                this->sph_bessel_camb(l,k2*rp) * integral3;
        };
        double integral2 = integrate_simps(integrand2, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
        return r*r / (this->H_f[n]*1000.0) * this->Tb_interp(z) *\
            this->sph_bessel_camb(l,k1*r) * integral2;
    };
    double integral1 = integrate_simps(integrand1,this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return pow(this->prefactor_Ml,2) * integral1;
}

