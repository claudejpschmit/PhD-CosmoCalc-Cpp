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


void SanityChecker::compare_interp(int l, double k1, double k2, double z, double *out1, double *out2)
{
    cout << "Compare using interpolation for z = " << z << ": " << endl;
    double k_high = 1.0;
    double k_low = 0.01;
    int steps = 2*(int)((k_high - k_low)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;
    double r,q;
    r = spline1dcalc(r_interp, z);
    q = spline1dcalc(q_interp, z);
    
    double hubble_z = spline1dcalc(H_f_interp,z)*1000.0;

    auto integrand2 = [&](double zp)
    {
        double rp,qp, hubble_zp;
        rp = spline1dcalc(r_interp,zp);
        qp = spline1dcalc(q_interp,zp);
        hubble_zp = spline1dcalc(H_f_interp,zp)*1000.0;
        auto integrand3 = [&](double kappa)
        {
            double hhh = pow(this->h,3);
            double sP = sqrt(this->Pk_interp(kappa*this->h,z)/hhh);
            double sPp = sqrt(this->Pk_interp(kappa*this->h,zp)/hhh);
            return kappa*kappa * sP * sPp * this->sph_bessel_camb(l,kappa*q) *\
                this->sph_bessel_camb(l,kappa*qp);
        };
        double integral3 = integrate_simps(integrand3, k_low, k_high, steps);
        return rp*rp / hubble_zp * this->Tb_interp(zp) *\
            this->sph_bessel_camb(l,k2*rp) * integral3;
    };
    int zstep = 400;
    double integral2 = integrate_simps(integrand2, z-0.2, z+0.2, zstep);
    double res1 = pow(this->prefactor_Ml,2) * r*r / hubble_z * this->Tb_interp(z) *\
                  this->sph_bessel_camb(l,k1*r) * integral2;

    double qp, rr; 
    qp = spline1dcalc(q_p_interp,z);
    rr = r*r;
    double hh = pow(hubble_z, 2);
    double A = rr * this->Pk_interp(((double)l + 0.5)/q * this->h,z)/(pow(this->h,3)*hh*qp) *\
               pow(this->Tb_interp(z),2);

    double pre = 2*this->b_bias*this->b_bias*this->c*this->c/this->pi;
    double res2 = pre * A * rr / (q*q) * this->sph_bessel_camb(l,k1*r) * this->sph_bessel_camb(l, k2*r);
    *out1 = res1;
    *out2 = res2;

    cout << res1 << endl;
    cout << res2 << endl;
    cout << res2 / res1 << endl;
}

void SanityChecker::compare_new(int l, double k1, double k2, double kappa, double *out1, double *out2)
{
    double res1 = pow(kappa,2) * this->M(l,k1,kappa) * this->M(l,k2,kappa);
    double res2 = pow(kappa,2) * 2 * pow(this->b_bias,2) * pow(this->c,2)/this->pi;
    double zstar1 = r_inverse((l+0.5)/k1);
    double zstar2 = r_inverse((l+0.5)/k2);

    cout << zstar1 << " " << zstar2 << endl;
    double H1 = spline1dcalc(H_f_interp_full, zstar1)*1000.0;
    double H2 = spline1dcalc(H_f_interp_full, zstar2)*1000.0;
    double Tb1 = Tb_interp_full(zstar1);
    double Tb2 = Tb_interp_full(zstar2);
    double q1 = spline1dcalc(q_interp_full, zstar1);
    double q2 = spline1dcalc(q_interp_full, zstar2);
    double hhh = pow(this->h,3);
    double sP1 = sqrt(this->Pk_interp_full(kappa*this->h,zstar1)/hhh);
    double sP2 = sqrt(this->Pk_interp_full(kappa*this->h,zstar2)/hhh);
    double qp1 = spline1dcalc(q_p_interp_full, zstar1);
    double qp2 = spline1dcalc(q_p_interp_full, zstar2);


    double jl1 = sph_bessel_camb(l,q1*kappa);
    double jl2 = sph_bessel_camb(l,q2*kappa);

    res2 = res2*pow(l+0.5,3)/(k1*k1*k1*k2*k2*k2* H1 * H2 * qp1 * qp2) * sP1 * sP2 * Tb1 * Tb2 * jl1 * jl2;
    *out1 = res1;
    *out2 = res2;

    cout << "Full calculation: " << res1 << endl;
    cout << "Simplified calculation: " << res2 << endl;


}

double SanityChecker::limber3(int l, double z)
{
    auto integrand2 = [&](double zp)
    {
        auto integrand1 = [&](double k)
        {
            return pow(k,3) * sph_bessel_camb(l, k*z) * sph_bessel_camb(l, k*zp);
        };

        return integrate_simps(integrand1,0,3000,100000);
    };

    double range = 0.5;
    double A = integrate_simps(integrand2,z-range,z+range, 10000);
    //double nu = pow(l+0.5,2); 
    double res = A * 2.0 /(this->pi * ((double)l+0.5));
    cout << res << endl;
    cout << 1.0/pow(z,3) << endl;
    return res;
}

void SanityChecker::kappa_integrand(int l, double k1, double k2, double z, double kappa, double zp, double *out1)
{
    double r,q;
    r = spline1dcalc(r_interp, z);
    q = spline1dcalc(q_interp, z);   
    double hubble_z = spline1dcalc(H_f_interp,z)*1000.0;
    double jz = this->sph_bessel_camb(l, k1 * r);

    double rp,qp, hubble_zp;
    rp = spline1dcalc(r_interp,zp);
    qp = spline1dcalc(q_interp,zp);
    hubble_zp = spline1dcalc(H_f_interp,zp)*1000.0;
    double jzp = this->sph_bessel_camb(l, k2 * r);

    double hhh = pow(this->h,3);
    double sP = sqrt(this->Pk_interp(kappa*this->h,z)/hhh);
    double sPp = sqrt(this->Pk_interp(kappa*this->h,zp)/hhh);
    
    double res1 = kappa*kappa * sP * sPp * this->sph_bessel_camb(l,kappa*q) *\
                  this->sph_bessel_camb(l,kappa*qp);

    double res2 = rp*rp / hubble_zp * this->Tb_interp(zp) * jzp * res1;

    double res3 = pow(this->prefactor_Ml,2) * r*r / hubble_z * this->Tb_interp(z) * jz * res2;
   
    *out1 = res1;
}

void SanityChecker::kappa_integral(int l, double k1, double k2, double z, double zp, double *out1, double k_low, double k_high)
{
    double r,q;
    r = spline1dcalc(r_interp, z);
    q = spline1dcalc(q_interp, z);   
    double hubble_z = spline1dcalc(H_f_interp,z)*1000.0;
    double jz = this->bessel_j_interp_cubic(l, k1 * r);

    double rp,qp, hubble_zp;
    rp = spline1dcalc(r_interp,zp);
    qp = spline1dcalc(q_interp,zp);
    hubble_zp = spline1dcalc(H_f_interp,zp)*1000.0;
    double jzp = this->bessel_j_interp_cubic(l, k2 * rp);

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
    //double integral2 = qromb(integrand3, k_low, k_high, 10E-10);
    //cout << integral2 << endl;
    *out1 = integral2;
}

void SanityChecker::zp_integrand(int l, double k1, double k2, double z, double zp, double *out1, double k_low, double k_high)
{
    double integral;
    kappa_integral(l, k1, k2, z, zp, &integral, k_low, k_high);
    
    double r,q;
    r = spline1dcalc(r_interp, z);
    q = spline1dcalc(q_interp, z);   
    double hubble_z = spline1dcalc(H_f_interp,z)*1000.0;
    double jz = this->bessel_j_interp_cubic(l, k1 * r);

    double rp,qp, hubble_zp;
    rp = spline1dcalc(r_interp,zp);
    qp = spline1dcalc(q_interp,zp);
    hubble_zp = spline1dcalc(H_f_interp,zp)*1000.0;
    double jzp = this->bessel_j_interp_cubic(l, k2 * rp);

    double res1 = integral * rp*rp/hubble_zp * this->Tb_interp(zp) * jzp;
    double res2 = pow(this->prefactor_Ml,2) * r*r / hubble_z * this->Tb_interp(z) * jz * res1;

    *out1 = res1;
   
}
void SanityChecker::zp_integrand2(int l, double k1, double k2, double z, double zp, double *out1, double k_low, double k_high)
{
    double integral;

    double r,q;
    r = spline1dcalc(r_interp, z);
    q = spline1dcalc(q_interp, z);   
    double jz = this->bessel_j_interp_cubic(l, k1 * r);

    double rp,qp, hubble_zp;
    rp = spline1dcalc(r_interp,zp);
    qp = spline1dcalc(q_interp,zp);
    hubble_zp = spline1dcalc(H_f_interp,zp)*1000.0;
    double jzp = this->bessel_j_interp_cubic(l, k2 * rp);

    int steps = 4*(int)((k_high - k_low)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;
    
    auto integrand3 = [&](double kap)
        {
            return this->sph_bessel_camb(l,kap*q) * this->sph_bessel_camb(l,kap*qp);
        };
    integral = integrate_simps(integrand3, k_low, k_high, steps);

    
    double res1 = integral * jzp;

    *out1 = res1;
   
}


void SanityChecker::Cl_compare(int l, double k1, double k2)
{
    double k_low = 0.1;
    double k_high = 1.0;
    double res_simple = this->Cl_simplified3(l, k1, k2);
    cout << "simple done " << res_simple << endl;
    double res_full = this->corr_Tb_new(l,k1,k2,k_low, k_high);
    cout << "full done " << res_full << endl;
    cout << "Cl_simple/Cl_full = " << res_simple/res_full << endl;
}
    
