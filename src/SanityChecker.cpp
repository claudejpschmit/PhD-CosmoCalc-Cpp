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
    double Cl_levin = Cl_gauss(l,k1,k2,k_low,k_high,0,0,0);
    t2 = clock();
    double ta = (double)t2 - (double)t1;
    cout << "C_levin = " << Cl_levin << endl;
    cout << "time = " << ta/CLOCKS_PER_SEC << endl;

    t1 = clock();
    double Cl_full = Cl_new(l,k1,k2,k_low,k_high,n_levin,0,0,0);
    t2 = clock();
    double tb = (double)t2 - (double)t1;
    cout << "C_full = " << Cl_full << endl;
    cout << "time = " << tb/CLOCKS_PER_SEC << endl;
    cout << "ratio = " << Cl_levin/Cl_full << endl;
    *ratio = Cl_levin/Cl_full;
    *time_r = tb/ta;
    cout << "time_ratio = " << tb/ta << endl;
}

void SanityChecker::plot_integrand_z(int l, double k1, double k2, int zsteps, string filename)
{
    int kappa_steps = (int)(abs(0.001 - 2.0)/this->k_stepsize);
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
            double integral3 = integrate_simps(integrand3, 0.001, 2.0, kappa_steps);
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

void SanityChecker::plot_intjj(int l, double zp, int zsteps, string filename)
{
    double low;
    double k_low = 0.0001;
    double k_high = 2.0;
    double hhh = pow(qs[0].h,3);
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
    ofstream file;
    file.open(filename);

    double rp = r_interp(zp);
    double qp = q_interp(zp, 0);
    
    auto integrand1 = [&](double kappa)
    {
        double sP = sqrt(this->Pk_interp(kappa*qs[0].h,zp, 0)/hhh);
        return kappa * kappa * sP * sP * pow(sph_bessel_camb(l,kappa*qp),2);
    };
    double integ = integrate_simps(integrand1, lower_kappa_bound, k_high, steps);

    double aa = integ; 
    double bb = zp;
    double delta = 2.26*pow((double)l,-0.677);
    cout << delta << endl;
    if (delta > 0.1)
        delta = 0.1;
    double qp_delta = q_interp(zp+delta,0);
    auto integrand2 = [&](double kappa)
    {
        double sP = sqrt(this->Pk_interp(kappa*qs[0].h,zp+delta, 0)/hhh);
        double sPp = sqrt(this->Pk_interp(kappa*qs[0].h,zp, 0)/hhh);
        return kappa*kappa*sP*sPp*this->sph_bessel_camb(l,kappa*qp_delta) *\
            this->sph_bessel_camb(l,kappa*qp);

    };
    double integI = integrate_simps(integrand2, lower_kappa_bound, k_high, steps);


    double cc2 = - delta*delta/(2.0*log(integI/aa));
    double zstepsize = 2.0*4.29*sqrt(cc2)/(double)zsteps;
    cout << zstepsize << endl;
    for (int i = 0; i < zsteps; i++)
    {
        double z = zp - 4.29*sqrt(cc2) + i*zstepsize;
        double q = q_interp(z, 0);
        double r = r_interp(z);

        auto integrand3 = [&](double kappa)
        {
            double sP = sqrt(this->Pk_interp(kappa*qs[0].h,z, 0)/hhh);
            double sPp = sqrt(this->Pk_interp(kappa*qs[0].h,zp, 0)/hhh);
            return kappa*kappa*sP*sPp*this->sph_bessel_camb(l,kappa*q) *\
                this->sph_bessel_camb(l,kappa*qp);

        };
        double integral3 = integrate_simps(integrand3, lower_kappa_bound, k_high, steps);
        
        double approx = aa * exp(-pow(z-bb,2)/(2.0*cc2));

        double hf = Hf_interp(z)*1000.0;
        double tb = Tb_interp(z,0);
        double pre = sph_bessel_camb(l, 1.5 * r) * r*r/hf *tb;
        file << z << " " << integral3 << " " << approx << endl;
    }
    file.close();
}


void SanityChecker::plot_integrand_zp(int l, double z, double k2, int zsteps, string filename)
{
    int kappa_steps = 4*(int)(abs(0.001 - 2.0)/this->k_stepsize);
    if (kappa_steps % 2 == 1)
        ++kappa_steps;
    ofstream file;
    file.open(filename);

    double zstepsize = 2.0/(double)zsteps;
    double q = q_interp(z, 0);
    double hhh = pow(qs[0].h,3);

    for (int i = 0; i < zsteps; i++)
    {
        double zp = 7.0 + i*zstepsize;
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
        cout << kappa_steps << endl;
        double integral3 = integrate(integrand3, 0.0001, 2.0, kappa_steps, simpson());
        cout << integral3 << endl;
        //double integral3 = integrate_simps(integrand3, 0.001, 2.0, kappa_steps);
        double result = rp*rp / (Hf_interp(zp)*1000.0) * this->Tb_interp(zp, 0) *\
            this->sph_bessel_camb(l,k2*rp) * integral3;
        
       
        file << zp << " " << result << endl;
    }
    file.close();
}

double SanityChecker::Cl_gauss(int l, double k1, double k2, double k_low,double k_high, int Pk_index, int Tb_index, int q_index)
{
    // This determines the lower bound of the kappa integral, basically it gets rid
    // of the kappas for which k < l/q, since the integral is close to zero there.
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

    // once the lower bound for kappa is found, we determine the best number of
    // integration steps.
    int steps = (int)(abs(k_high - lower_kappa_bound)/this->k_stepsize);
    if (steps % 2 == 1)
        ++steps;

    // delta is independent of kappa, z & zp so we can calculate it once.
    //double delta = pow(pi/(double)l,0.5);
    double delta = 2.26*pow((double)l,-0.677);
    if (delta > 0.1)
        delta = 0.1;
    // Now comes the integral.
    auto integrand_z = [&](double z)
    {
        double q = q_interp(z,q_index);
        double r = r_interp(z);
        //calculate I(z,delta) & change delta if necessary
        double q_delta = q_interp(z+delta, q_index);
        auto I = [&](double kappa)
        {
            double sP = sqrt(this->Pk_interp(kappa*qs[q_index].h,z, Pk_index)/hhh);
            double sPp = sqrt(this->Pk_interp(kappa*qs[q_index].h,z+delta, Pk_index)/hhh);
            return kappa*kappa*sP*sPp*this->sph_bessel_camb(l,kappa*q) *\
                this->sph_bessel_camb(l,kappa*q_delta);
        };
        double Izd = integrate_simps(I, lower_kappa_bound, k_high, steps);
        // in case I(z,delta) is negative we need to fix at a different delta, so that ln(I) exists.
        while (Izd <= 0) {
            cout << "Izd calculated again, from delta = " << delta << " to delta = ";
            delta -= 0.0001;
            cout << delta << endl;
            q_delta = q_interp(z+delta, q_index);
            Izd = integrate_simps(I, lower_kappa_bound, k_high, steps);  
        }
        
        // Calculate a_l(z):
        auto a_integrand = [&](double kappa)
        {
            double P = this->Pk_interp(kappa*qs[q_index].h, z, Pk_index)/hhh;
            double jlkap = sph_bessel_camb(l, kappa*q);
            return kappa*kappa * P * jlkap * jlkap;
        };
        double alz = integrate_simps(a_integrand, lower_kappa_bound, k_high, steps); 
        // calculate ln(I/a):
        double lnIA = log(Izd/alz);
        // and ln(I/a)/delta^2:
        double frac = lnIA/pow(delta,2);
        // This concludes calculating all z dependencies for zp integral, so:

        auto integrand_zp = [&](double zp)
        {
            double rp = r_interp(zp);
            double hfp = Hf_interp(zp) * 1000.0;
            double Tbp = Tb_interp(zp, Tb_index);
            double jlp = sph_bessel_camb(l, k2*rp);
            double exponential = exp(pow(z-zp,2)*frac);
            return rp*rp/hfp * Tbp * jlp * exponential;
        };
        double integral_zp = integrate_simps(integrand_zp, this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
        double hf = Hf_interp(z) * 1000.0;
        double Tb = Tb_interp(z, Tb_index);
        double jl = sph_bessel_camb(l, k1 * r);
        return r*r/hf * Tb * jl * alz * integral_zp;
    };
    double integral_z = integrate_simps(integrand_z,this->zmin_Ml, this->zmax_Ml, this->zsteps_Ml);
    return pow(this->prefactor_Ml,2) * integral_z;
}
