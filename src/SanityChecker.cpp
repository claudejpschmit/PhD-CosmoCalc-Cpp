#include "SanityChecker.hpp"
#include "Integrator.hpp"
#include <time.h>

SanityChecker::SanityChecker(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index)
    :
        CosmoCalc(params,Pk_index,Tb_index,q_index)
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
    q = q_interp(z,0);   

    double qp;
    qp = q_interp(zp,0);

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
        double sP1 = sqrt(this->Pk_interp(kap*this->h,z,0)/hhh);
        double sPp1 = sqrt(this->Pk_interp(kap*this->h,zp,0)/hhh);
        return kap*kap * sP1 * sPp1 * this->sph_bessel_camb(l,kap*q) *\
            this->sph_bessel_camb(l,kap*qp);
    };
    double integral2 = integrate_simps(integrand3, lower_kappa_bound, k_high, steps);

    Levin LEVIN(a, k_high);

    auto foo = [&](double kappa)
    {
        double sP1 = sqrt(this->Pk_interp(kappa*this->h,z,0)/hhh);
        double sPp1 = sqrt(this->Pk_interp(kappa*this->h,zp,0)/hhh);
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
        integral3 = LEVIN.integrate_2sphj_1r(foo,q,l,n);
    else
        integral3 = LEVIN.integrate_2sphj_2r(foo,q,qp,l,n);

    integral3 += integral;
    *out1 = integral2/integral3;
}

double SanityChecker::kappa_integrand(int l, double z, double zp, double kappa)
{
    double q = q_interp(z,0);   
    double qp = q_interp(zp,0);

    double hhh = pow(this->h,3);
    double sP1 = sqrt(this->Pk_interp(kappa*this->h,z,0)/hhh);
    double sPp1 = sqrt(this->Pk_interp(kappa*this->h,zp,0)/hhh);
    return kappa*kappa * sP1 * sPp1 * this->sph_bessel_camb(l,kappa*q) *\
        this->sph_bessel_camb(l,kappa*qp);
}

void SanityChecker::Compare_Cl(int l, double k1, double k2, double k_low, double k_high, int n_levin, double *ratio, double *time_r)
{
    clock_t t1,t2;
    t1 = clock();
    double Cl_levin = Cl_new(l,k1,k2,k_low,k_high,n_levin,0,0,0);
    t2 = clock();
    double ta = (double)t2 - (double)t1;
    cout << "C_levin = " << Cl_levin << endl;
    cout << "time = " << ta/CLOCKS_PER_SEC << endl;

    t1 = clock();
    double Cl_full = corr_Tb_new(l, k1, k2, k_low, k_high,0,0,0);
    t2 = clock();
    double tb = (double)t2 - (double)t1;
    cout << "C_full = " << Cl_full << endl;
    cout << "time = " << tb/CLOCKS_PER_SEC << endl;
    cout << "ratio = " << Cl_levin/Cl_full << endl;
    *ratio = Cl_levin/Cl_full;
    *time_r = tb/ta;
    cout << "time_ratio = " << tb/ta << endl;
}

void SanityChecker::plot_integrad_z(int l, double k1, double k2, int zsteps, string filename)
{
    int kappa_steps = (int)((0.001 - 2.0)/this->k_stepsize);
    if (kappa_steps % 2 == 1)
        ++kappa_steps;
    ofstream file;
    file.open(filename);

    double zstepsize = 2.0/(double)zsteps;
    double hhh = pow(qs[0].h,3);
    for (int i = 0; i < zsteps; i++)
    {
        double z = 7.0 + i*zstepsize;
        double r,q;
        r = r_interp(z);
        q = q_interp(z, 0);


        auto integrand2 = [&](double zp)
        {
            double rp,qp;
            rp = r_interp(zp);
            qp = q_interp(zp, 0);


            auto integrand3 = [&](double kappa)
            {
                double sP = sqrt(this->Pk_interp(kappa*qs[0].h,z, 0)/hhh);
                double sPp = sqrt(this->Pk_interp(kappa*qs[0].h,zp, 0)/hhh);
                return kappa*kappa * sP * sPp * this->sph_bessel_camb(l,kappa*q) *\
                    this->sph_bessel_camb(l,kappa*qp);
            };
            double integral3 = integrate_simps(integrand3, 0.3, 2.0, kappa_steps);
            return rp*rp / (Hf_interp(zp)*1000.0) * this->Tb_interp(zp, 0) *\
                this->sph_bessel_camb(l,k2*rp) * integral3;
        };
        double integral2 = integrate_simps(integrand2, 7.0, 9.0, zsteps);
        double result = r*r / (Hf_interp(z)*1000.0) * this->Tb_interp(z,0) *\
            this->sph_bessel_camb(l,k1*r) * integral2;
        
        file << z << " " << result << endl;
    }
    file.close();

}
