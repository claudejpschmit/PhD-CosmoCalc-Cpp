#include "CosmoBasis.hpp"
#include <iostream>

using namespace std;

CosmoBasis::CosmoBasis(map<char*,double> params)
    :
        c(299792458.0),
        k_b(1.3806488*pow(10,-23)),
        m_b(1.674927351*pow(10,-27)),
        m_e(9.10938291*pow(10,-31)),
        e(1.60217657*pow(10,-19)),
        h_planck(6.62606957*pow(10,-34)),
        G(6.67384*pow(10,-11))
{
    T_star = h_planck * c / (k_b *0.21);
    this->params = params;
}

void CosmoBasis::show_params()
{
    cout << c << endl;
    cout << k_b << endl;
    cout << m_b << endl;
    cout << m_e << endl;
    cout << h_planck << endl;
    cout << e << endl;
    cout << G << endl;
    cout << T_star << endl;
}
