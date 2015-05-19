#include <iostream>
#include <string>
#include <map>
#include <cmath>
#include "CosmoBasis.hpp"
#include <fstream>
#include <sys/stat.h>
#include "Integrator.hpp"
#include "CosmologyCalculatorClass.hpp"
#include "ClassEngine.hpp"

using namespace std;

double f (double x)
{
    return x*x;
}

int main(int argc, char* argv[])
{
    (void) argc;
    (void) argv;
    struct stat sb;

    const int l_max_scalars=5000;

    /*
    ClassParams pars;
    //pars.add("H0",70.3);
    pars.add("100*theta_s",1.04);
    pars.add("omega_b",0.0220);
    pars.add("omega_cdm",0.1116);
    pars.add("A_s",2.42e-9);
    pars.add("n_s",.96);
    pars.add("tau_reio",0.09);

    pars.add("k_pivot",0.05);
    pars.add("YHe",0.25);
    pars.add("output","mPk"); //pol +clphi
    pars.add("P_k_max_h/Mpc", 100);
    pars.add("z_pk", 2.0);
    //pars.add("l_max_scalars",l_max_scalars);
    //pars.add("lensing",false); //note boolean

    ClassEngine * KKK(0);
    cout << "before initialize" << endl;
    try{
        //le calculateur de spectres
        if (argc==2){
            string pre=string(argv[1]);
            KKK=new ClassEngine(pars,pre);
        }
        else{
            KKK=new ClassEngine(pars);
        }

        cout << "after initialize" << endl;
        ofstream file;
        file.open("outputtest.dat");
        cout.precision( 16 );
        KKK->writePks(file);
        file.close();
        //KKK->writeCls(cout);
    }
    catch (std::exception &e){
        cout << "GIOSH" << e.what() << endl;
    }
    */
    const char output_path[] = "output";
    if (stat(output_path, &sb) == 0 && S_ISDIR(sb.st_mode)) {
        cout << "output directory already exists." << endl;
    } else {
        mkdir(output_path, 0700);
        cout << "output directory generated!" << endl; 
    }

    auto g = [](double x) 
    {
        return 2*f(x);
    };
    cout << "integration yields "<< integrate(g, 0.0, 1.0, 100, simpson()) << endl;

    map<string,double> params;

    params["O"] = 0.02;
    params["ombh2"] = 0.02;
    cout << params["O"] << endl;
    CosmoBasis base(params);
    base.show_params(); 
    //cout << 1.34 * pow(10,2) << endl;

    ofstream output;
    const char filename[] = "bessels.dat";
    char path[50] = "";
    strcat(path, output_path);
    strcat(path, "/");
    strcat(path, filename);
    cout << path << endl;
    output.open(path);
    for (int n = 0; n < 100000; ++n) {
        output << n << " " << base.sph_bessel(10,double(n)/100.0) << endl;
    }
    output.close();

    CosmoCalc calc(params);
    calc.show_cosmo_calcs();

    return 0;
}
