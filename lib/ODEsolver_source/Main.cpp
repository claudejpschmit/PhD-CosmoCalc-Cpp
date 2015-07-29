#include "ODEs.hpp"
#include "ODE_Solver.hpp"

#include <iostream>
#include <fstream>

#define DT 0.05
#define T0 0.0

using namespace std;

int main( ) {

    ExponentialODE E_ODE(T0);
    PolynomialODE P_ODE(T0);

    EulerMethod E_method_1(DT, &E_ODE, T0);
    MRKMethod M_method_1(DT, &E_ODE, T0);
    FORKMethod F_method_1(DT, &E_ODE, T0);
    EulerMethod E_method_2(DT, &P_ODE, T0);
    MRKMethod M_method_2(DT, &P_ODE, T0);
    FORKMethod F_method_2(DT, &P_ODE, T0);

    ofstream euler_1("euler1.txt");
    ofstream mkm_1("mkm1.txt");
    ofstream fork_1("fork1.txt");
    ofstream euler_2("euler2.txt");
    ofstream mkm_2("mkm2.txt");
    ofstream fork_2("fork2.txt");
    ofstream exact1("exact1.txt");
    ofstream exact2("exact2.txt");

    for (int i = 1; i < 100; ++i) {
        euler_1 << i * DT << " " << E_method_1.step() << endl;
        mkm_1 << i * DT << " " << M_method_1.step() << endl;
        fork_1 << i * DT << " " << F_method_1.step() << endl;
        euler_2 << i * DT << " " << E_method_2.step() << endl;
        mkm_2 << i * DT << " " << M_method_2.step() << endl;
        fork_2 << i * DT << " " << F_method_2.step() << endl;
        exact1 << i * DT << " " << E_ODE.exact_solution(i * DT) << endl;
        exact2 << i * DT << " " << P_ODE.exact_solution(i * DT) << endl;
    }

    euler_1.close();
    mkm_1.close();
    fork_1.close();
    euler_2.close();
    mkm_2.close();
    fork_2.close();
    exact1.close();
    exact2.close();

    return 0;
}
