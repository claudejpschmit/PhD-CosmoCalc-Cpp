#include "SanityChecker.hpp"
#include "Integrator.hpp"

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

void SanityChecker::kappa_integral(int l, double k1, double k2, double z, double zp, double *out1, double k_low, double k_high)
{
    double q;
    q = spline1dcalc(q_interp, z);   

    double qp;
    qp = spline1dcalc(q_interp,zp);

    double hhh = pow(this->h,3);

    int steps = 4*(int)((k_high - k_low)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;
    


    auto integrand3 = [&](double kap)
    {
        double sP1 = sqrt(this->Pk_interp(kap*this->h,z)/hhh);
        double sPp1 = sqrt(this->Pk_interp(kap*this->h,zp)/hhh);
        return kap*kap * sP1 * sPp1 * this->sph_bessel_camb(l,kap*q) *\
            this->sph_bessel_camb(l,kap*qp);
    };
    double integral2 = integrate_simps(integrand3, k_low, k_high, steps);
    cout << "first done" << endl;
    LEVIN = new Levin(k_low, k_high);

    auto foo = [&](double kappa)
    {
        double sP1 = sqrt(this->Pk_interp(kappa*this->h,z)/hhh);
        double sPp1 = sqrt(this->Pk_interp(kappa*this->h,zp)/hhh);
        return kappa*kappa * sP1 * sPp1;     
    };
    int n = 10;
    double integral3 = LEVIN->integrate_2sphj_2r(foo,q,qp,l,n);

    delete LEVIN;
    *out1 = integral2/integral3;
}
