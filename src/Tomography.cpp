#include "Tomography.hpp"
#include "Integrator.hpp"

Tomography::Tomography(map<string, double> params)
    :
    CosmoBasis(params)
{
    cout << "Now Building Tomography" << endl;
    cout << "Building interpolator for kappa_10(TK)" << endl;
    //table contains the values tabulated in Allison & Dalgarno 1969
    real_1d_array xs, ys;
    xs.setlength(19);
    ys.setlength(19);
    xs[0] = 1;
    xs[1] = 2;
    xs[2] = 4;
    xs[3] = 6;
    xs[4] = 8;
    xs[5] = 10;
    xs[6] = 15;
    xs[7] = 20;
    xs[8] = 25;
    xs[9] = 30;
    xs[10] = 40;
    xs[11] = 50;
    xs[12] = 60;
    xs[13] = 70;
    xs[14] = 80;
    xs[15] = 90;
    xs[16] = 100;
    xs[17] = 200;
    xs[18] = 300;

    ys[0] = 2.2E-14;
    ys[1] = 4.2E-14;   
    ys[2] = 1.8E-13;
    ys[3] = 5.1E-13;
    ys[4] = 1.2E-12;   
    ys[5] = 2.3E-12;
    ys[6] = 7.4E-12;
    ys[7] = 1.5E-11;   
    ys[8] = 2.3E-11;
    ys[9] = 3.0E-11;
    ys[10] = 4.4E-11;   
    ys[11] = 5.6E-11;
    ys[12] = 6.6E-11;
    ys[13] = 7.4E-11;   
    ys[14] = 8.2E-11;
    ys[15] = 8.9E-11;
    ys[16] = 9.5E-11;   
    ys[17] = 1.4E-10;   
    ys[18] = 1.6E-10;   
      
    spline1dbuildlinear(xs,ys,kappa10_interpolator);

    cout << "Foreground initialization" << endl;

    Ai.push_back(10.0);
    Ai.push_back(0.014);
    Ai.push_back(700.0);
    Ai.push_back(0.088);
    
    betai.push_back(1.1);
    betai.push_back(1.0);
    betai.push_back(2.4);
    betai.push_back(3.0);
    
    alphai.push_back(2.07);
    alphai.push_back(2.1);
    alphai.push_back(2.8);
    alphai.push_back(2.15);
    
    xii.push_back(1.0);
    xii.push_back(35.0);
    xii.push_back(4.0);
    xii.push_back(35.0);



}
Tomography::~Tomography()
{}
double Tomography::Cl(int l, double f1, double f2)
{
    return Cl_s(l,f1,f2) + Cl_f(l,f1,f2) + Cl_n(l,f1,f2);
}
double Tomography::Cl_s(int l, double f1, double f2)
{
    auto integrand = [&](double k)
    {
        double z1 = redshift_z(f1);
        double z2 = redshift_z(f2);
        double res = k*k*P_dcdc(k,f1,f2)* (F(k,f1)*F(k,f2) * I(l,k,f1) * I(l,k,f2) +\
                f(z1)*f(z2) * J(l,k,f1) * J(l,k,f2) -\
                F(k,f1)*f(z2) * I(l,k,f1) * J(l,k,f2) -\
                F(k,f2)*f(z1) * I(l,k,f2) * J(l,k,f1) );
        return res;
    };
    double kmin = 0;
    double kmax = 100;
    double In = integrate(integrand, kmin, kmax, 1000, simpson());
    return 2.0/pi*In * T21(f1) * T21(f2);
}
double Tomography::Cl_n(int l, double f1, double f2)
{
    //TODO: Do the noise.
    if (f1 == f2)
        return 1;
    else
        return 0;
}
double Tomography::Cl_f(int l, double f1, double f2)
{
    double res = 0;
    for (int i = 0; i < 4; i++)
    {
        res = Cl_f(i,l,f1,f2);
    }
    return res;
}

// ------------ Private Functions -------------- //
double Tomography::T21(double f)
{
    double zz = redshift_z(f);
    double result = Tc(zz) * ytot(zz)/(1+ytot(zz)) * (1-T(zz)/Tk(zz)); 
    return result;
}
double Tomography::Tc(double z)
{
    return 23.0 * xHI(z)*(0.7/h)*(O_b*h*h/0.02) * sqrt(0.15/(O_M*h*h) * (1+z)/10.0);
}
double Tomography::F(double k, double f)
{
    double zz = redshift_z(f);
    double bb = b_b(k,zz);
    double gg = _g(zz);
    double RL= 1; //TODO: right value...
    double R = 1; //TODO: right value... R = characteristic size of reionization patches in MPc
    double res = bb + b_xH(k,zz)*exp(-k*k * R*R /2.0) + T(zz)/(Tk(zz) - T(zz)) *\
                 gg*bb + b_alpha(k,zz)/(1+ytot(zz))*exp(-k*k * RL*RL/2.0);
    return res;
}
double Tomography::I(int l, double k, double f)
{
    auto integrand = [&](double zp)
    {
        //currently in MHz
        double dfdz = - 1420/((1+zp)*(1+zp)); 
        double fp = 1420/(1+zp);
        double rr = r(zp);
        double res = dfdz * W(f, fp) * sph_bessel_camb(l, k*rr);
        return res;
    };
    double In = integrate(integrand, 0.0, 100.0, 1000, simpson());
    return In;
}
double Tomography::J(int l, double k, double f)
{
    auto integrand = [&](double zp)
    {
        //currently in MHz
        double dfdz = - 1420/((1+zp)*(1+zp));
        double fp = 1420/(1+zp);
        double rr = r(zp);
        double jm2 = sph_bessel_camb(l-2, k*rr);
        double jm1 = sph_bessel_camb(l-1, k*rr);
        double j = sph_bessel_camb(l, k*rr);
        double jp1 = sph_bessel_camb(l+1, k*rr);
        double jp2 = sph_bessel_camb(l+2, k*rr);
        double deriv_j = 0.25 * (jm2 + jp2) + 1.0/(2.0*k*rr) * (jp1 - jm1) +\
                         (3 - 2*(k*k*rr*rr))/(4*k*k*rr*rr) * j;
        double res = dfdz * W(f, fp) * deriv_j;
        return res;
    };
    double Jn = integrate(integrand, 0.0, 100.0, 1000, simpson());
    return Jn;
}
double Tomography::f(double z)
{
    // f(z) = dlnG(z)/dlna(z)
    // derivative (P_growth(z), lna(z))
    // deriv = d lnG(z)/dz
    auto lnG = [&](double z)
    {
        return log(P_growth(z));
    };
    double h = 0.0001;
    double deriv = (lnG(z+h)-lnG(z-h))/(2*h);
    return (-(1+z) * deriv);
}
double Tomography::W(double f0, double f)
{
    return 0;
}

double Tomography::redshift_z(double f)
{
    return (f21/f - 1);
}

double Tomography::T(double z)
{
    return this->T_CMB * (1+z);
}
double Tomography::Tk(double z)
{
    // Adiabatic cooling starts at z = 145;
    double A = 397.85/(146.0*146.0);
    return A * (1+z)*(1+z);
}
double Tomography::ytot(double z)
{
    return yc(z) + yalpha(z);
}
double Tomography::yc(double z)
{
    // collisions are only important at redshifts above 30, 
    // so set the coupling to 0 below.
    if (z >= 30)
        return 4*kappa10(Tk(z))*nH(z)*T_star/(3.0 * A_10 * T(z));
    else 
        return 0;
}
double Tomography::kappa10(double Tk)
{
    return spline1dcalc(kappa10_interpolator, Tk);
}
double Tomography::nH(double z)
{
    //TODO: still need to write this, though I 
    return 0;
}
double Tomography::yalpha(double z)
{
    //TODO: Important! Find an expression for this 
    return 0;
}
double Tomography::xHI(double z)
{
    //TODO: Important! Find an expression for this 
    return 0;
}

double Tomography::b_b(double k, double z)
{
    //TODO: Important! Find an expression for this
    return 0;
}

double Tomography::b_alpha(double k, double z)
{
    //TODO: Important! Find an expression for this
    return 0;
}

double Tomography::b_xH(double k, double z)
{
    (void)k;
    (void)z;
    return 0;
}

double Tomography::_g(double z)
{
    return 2.0/3.0;
}

double Tomography::P_dcdc(double k, double f1, double f2)
{
    double z1 = redshift_z(f1);
    double z2 = redshift_z(f2);
    return sqrt(Pkz_calc(k, z1) * Pkz_calc(k, z2));
}

double Tomography::r(double z)
{
    return D_H * Z(z);
}
double Tomography::Cl_f(int i, int l, double f1, double f2)
{
    double nu_f = 130; //MHz
    double Iii = exp(-log(f1/f2) * log(f1/f2) / (2.0 * xii[i]*xii[i]));
    double C1 = Ai[i] * pow(1000.0/(double)l, betai[i]) * pow(nu_f/f1, 2.0*alphai[i]);
    double C2 = Ai[i] * pow(1000.0/(double)l, betai[i]) * pow(nu_f/f2, 2.0*alphai[i]);
    
    return Iii * sqrt(C1*C2);
}


double Tomography::Pkz_calc(double k, double z)
{
    return this->P_growth(z) * this->P_delta(k);
}

double Tomography::P_growth(double z)
{
    const double res = this->D1(z) / this->D1(0);
    return pow(res,2);
}

double Tomography::D1(double z)
{
    const double prefactor = 5 * this->O_M / 2 * this->E(z);
    auto integrand = [&](double x)
    {
        return (1+x)/pow(this->E(x),3);
    };
    return prefactor * integrate(integrand, z, pow(10, 5), 1000000, simpson());
}

double Tomography::P_delta(double k, string units_k, string units_P)
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
        delta_H = 6.229 * pow(10,-5);
    } else {
        delta_H = 6.229 * pow(10,-5);
    }

    const double A = 2*pow(pi,2) * pow(this->c/1000.0,4) * pow(delta_H,2)/pow(10,8);
    const double x = k/keq;
    const double transfer_function_sq = pow(this->transfer(x),2);

    return A * k_factor * transfer_function_sq;
}

double Tomography::transfer(double x)
{
    const double res = log(1 + 0.171 * x) / (0.171 * x);
    const double bracket = 1 + 0.284 * x + pow(1.18 * x, 2) + pow(0.399 * x, 3) +\
                           pow(0.49 * x, 4);
    return res * pow(bracket, -0.25);
}
