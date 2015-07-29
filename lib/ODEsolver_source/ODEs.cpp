#include "ODEs.hpp"
#include <cmath>

ODE::~ODE()
{
}

double ODE::first_derivative(double y, double t)
{
    (void)y;
    (void)t;
    return 0;
}

double ODE::initial_value()
{
    return 0;
}

double ODE::exact_solution(double t)
{
    (void)t;
    return 0;
}

/* ****      Polynomial      **** */


PolynomialODE::PolynomialODE(double t0)
{
    this->t0 = t0;
}
double PolynomialODE::first_derivative(double y, double t)
{
    (void)y;
    return 3 * t * t - 6 * t + 1;
}
    
double PolynomialODE::initial_value()
{
    return exact_solution(t0);
}

double PolynomialODE::exact_solution(double t)
{
    return 1 + t - 3 * t * t + t * t * t;
}

/* ****     Exponential     **** */

ExponentialODE::ExponentialODE(double t0)
{
    this->t0 = t0;
}

double ExponentialODE::first_derivative(double y, double t)
{
    (void)y;
    return exp(t);
}
double ExponentialODE::initial_value()
{
    return exact_solution(t0);
}
double ExponentialODE::exact_solution(double t)
{
    return exp(t);
}

/* ****    Levin Example    **** */

LevinODE::LevinODE(double t0)
{
    this->t0 = t0;
}

double LevinODE::first_derivative(double y, double t)
{
    (void)y;
    return exp(t);
}
double LevinODE::initial_value()
{
    return 0;
}
double LevinODE::exact_solution(double t)
{
    return 0;
}
