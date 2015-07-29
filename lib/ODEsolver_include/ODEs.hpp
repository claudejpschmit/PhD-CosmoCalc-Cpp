#pragma once

#include <iostream>

class ODE {
        
public:
    virtual ~ODE();
    virtual double first_derivative(double y, double t);
    virtual double initial_value();
    virtual double exact_solution(double t);

protected:
    double t0;
};

class PolynomialODE : public ODE {
    
public:
    PolynomialODE(double t0);

    double first_derivative(double y, double t);
    double initial_value();
    double exact_solution(double t);

};

class ExponentialODE : public ODE {

public:
    ExponentialODE(double t0);

    double first_derivative(double y, double t);
    double initial_value();
    double exact_solution(double t);


};

class LevinODE : public ODE {

public:
    LevinODE(double t0);

    double first_derivative(double y, double t);
    double initial_value();
    double exact_solution(double t);


};
