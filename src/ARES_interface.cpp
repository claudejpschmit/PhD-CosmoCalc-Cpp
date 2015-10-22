#include "ARES_interface.hpp"
#include <iostream>
#include <fstream>

// What to do about Dark Energy w???

AresInterface::AresInterface()
    :
        omega_m_0(-1),
        omega_b_0(-1),
        omega_l_0(-1),
        hubble_0(-1),
        helium_by_mass(-1),
        cmb_temp_0(-1),
        sigma_8(-1),
        primordial_index(-1),
        fstar(-1),
        Tmin(-1),
        Nion(-1),
        fesc(-1),
        Nlw(-1),
        cX(-1),
        fX(-1)
{}
AresInterface::~AresInterface()
{}

void AresInterface::updateAres(map<string,double> params)
{
    // This is hubble_0
    hubble_0 = params["hubble"] / 100.0;

    // This is omega_b_0
    omega_b_0 = params["ombh2"] / (h*h);

    cmb_temp_0 = params["T_CMB"];
    double O_cdm = params["omch2"] / pow(hubble_0,2);
    double O_nu = params["omnuh2"] / pow(hubble_0,2);
    double O_gamma = pow(pi,2) * pow(cmb_temp_0/11605.0,4) / (15.0*8.098*pow(10,-11)*pow(hubble_0,2));
    double O_nu_rel = O_gamma * 3.0 * 7.0/8.0 * pow(4.0/11.0, 4.0/3.0);
    double O_R = O_gamma + O_nu_rel;
    double O_k = params["omk"];
    double O_tot = 1.0 - O_k;

    // This parameter is currently not used.
    double w = params["w_DE"];

    // This is omega_m_0
    omega_m_0 = omega_b_0 + O_cdm + O_nu;
    // This is omega_lambda_0
    omega_l_0 = O_tot - omega_m_0 - O_R;
    // This is primordial_index
    primordial_index = params["n_s"];
    // This is sigma_8
    sigma_8 = params["sigma8"];
    

    omNu = O_nu;

    fstar = params["fstar"];
    fesc = params["fesc"];
    Nion = params["nion"];
    fX = params["fx"];


    // Call run_ares.py with the necessary parameters
    stringstream command;
    command << "python run_ares.py";
    if (omega_m_0 != -1)
        command << " --omega_m_0 " << omega_m_0;
    if (omega_b_0 != -1)
        command << " --omega_b_0 " << omega_b_0; 
    if (omega_l_0 != -1)
        command << " --omega_l_0 " << omega_l_0; 
    if (hubble_0 != -1)
        command << " --hubble_0 " << hubble_0; 
    if (helium_by_mass != -1)
        command << " --helium_by_mass " << helium_by_mass; 
    if (cmb_temp_0 != -1)
        command << " --cmb_temp_0 " << cmb_temp_0; 
    if (sigma_8 != -1)
        command << " --sigma_8 " << sigma_8; 
    if (primordial_index != -1)
        command << " --primordial_index " << primordial_index; 
    if (fstar != -1)
        command << " --fstar " << fstar; 
    if (Tmin != -1)
        command << " --Tmin " << Tmin; 
    if (Nion != -1)
        command << " --Nion " << Nion; 
    if (fesc != -1)
        command << " --fesc " << fesc; 
    if (Nlw != -1)
        command << " --Nlw " << Nlw; 
    if (cX != -1)
        command << " --cX " << cX;
    if (fX != -1)
        command << " --fX " << fX; 
    
    char* command_buff = new char[command.str().length() + 1];
    strcpy(command_buff, command.str().c_str());
    int r = system(command_buff);
    (void)r;
}

void AresInterface::getTb(vector<double>* zp, vector<double>* Tbp)
{
    /*
    *zp = Tb_z;
    *Tbp = Tb;
    */
}

void AresInterface::calc_Tb(double zmin, double zmax, int zsteps)
{
    /*
    ofstream fout;
    fout.open("output/Tb_g.dat");
    Tb_z.clear();
    Tb.clear();
    double *result;
    double tk,lyaflux,xi,xe;
    double tb;

    result=dvector(1,3);

    //file="xc_history.dat";
    //fout.open(file);

    double z, stepsize;
    stepsize = (zmax-zmin)/(double)zsteps;

    for (int i = 0; i < zsteps; i++) {
        z = zmin + i * stepsize;
        a->getTIGM(z,result);
        tk=result[1];
        xi=result[2];
        xe=result[3];
        lyaflux=a->lyaFlux(z);
        tb=(1.0-xi)*tocm->tBrightGen(z,tk,xe,lyaflux);
        fout << z << " " << tb << endl;
        Tb_z.push_back(z);
        Tb.push_back(tb);
    }
    free_dvector(result,1,3);
    fout.close();
    */
}
