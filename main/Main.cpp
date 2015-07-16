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
#include "CosmologyWriterClass.hpp"
#include <time.h>
#include "FisherClass.hpp"
#include <armadillo>
#include "stdafx.h"
#include "interpolation.h"
#include "CAMB_interface.hpp"
#include "Global21cmInterface.hpp"

//driver
#include "math.h"
#include "astrophysics.h"
#include "dnumrecipes.h"
#include "dcosmology.h"
#include "twentyonecm.h"
#include "spline.h"

using namespace std;
using namespace arma;
using namespace alglib;

double f (double x)
{
    return x*x;
}

double fl (double x)
{
    return 0.551742 * x + 9778.15;
}

int main(int argc, char* argv[])
{ 


/*
    // Cosmology parameters
  double om0(0.3);
  double lam0(0.7);
  double omb(0.046);
  double h(0.7);
  double s8(0.9);
  double n(1.0);
  double omNu(0.0);

  double zin_arg(15.0);
  int tin(2);
  int i,j;
  char *file;
  ofstream fout;
  ifstream fin;
  int xin(1);
  int lyaxray_in(0);
  int popflag_in(0);
  int lyaxrayflag, popflag, xrayflag, sourceflag;
  double fstar, fx, nion, flya, fesc;
  int file_flag(0);

  // Handle arguments 
  while ((argc>1) && (argv[1][0]=='-')) {
    switch (argv[1][1]) {
    case 'z':
      zin_arg=atof(&argv[1][2]);
      break;
    case 't':
      tin=atoi(&argv[1][2]);
      break;
    case 'x':
      xin=atoi(&argv[1][2]);
      break;
    case 'l':
      lyaxray_in=atoi(&argv[1][2]);
      break;
    case 'p':
      popflag_in=atoi(&argv[1][2]);
      break;
    case 'f':
      file_flag=1;
      break;
    default:
      cerr << "Bad option" <<argv[1] <<"\n";
    }
    --argc;
    ++argv;
  }

  
  s8=0.77;
  omb=0.044;
  om0=0.26;
  lam0=0.74;
  h=0.72;
  n=0.95;

  xrayflag=-1;
  sourceflag=-1;

  if(file_flag==0){
    cout<<"Enter cosmology parameters"<<endl;
    cout<<"Enter s8:"<<endl;
    cin>>s8;
    cout<<"Enter omb:"<<endl;
    cin>>omb;
    cout<<"Enter om0:"<<endl;
    cin>>om0;
    cout<<"Enter lam0:"<<endl;
    cin>>lam0;
    cout<<"Enter h:"<<endl;
    cin>>h;
    cout<<"Enter n:"<<endl;
    cin>>n;
    cout<<"Enter astrophysics parameters"<<endl;
    cout<<"Enter fstar:"<<endl;
    cin>>fstar;
    cout<<"Enter fesc:"<<endl;
    cin>>fesc;
    cout<<"Enter nion:"<<endl;
    cin>>nion;
    cout<<"Enter fx:"<<endl;
    cin>>fx;
    cout<<"Enter flya:"<<endl;
    cin>>flya;
    cout<<"Enter popflag:"<<endl;
    cin>>popflag;
    cout<<"Enter xrayflag:"<<endl;
    cin>>xrayflag;
    cout<<"Enter lyaxrayflag:"<<endl;
    cin>>lyaxrayflag;
  }else{
    cin>>s8;
    cin>>omb;
    cin>>om0;
    cin>>lam0;
    cin>>h;
    cin>>n;
    cin>>fstar;
    cin>>fesc;
    cin>>nion;
    cin>>fx;
    cin>>flya;
    cin>>popflag;
    cin>>xrayflag;
    cin>>lyaxrayflag;
  }

  cout<<s8<<"\t"<<omb<<"\t"<<om0<<"\t"<<lam0<<"\t"<<h<<"\t"<<n<<endl;

  Cosmology c(om0,lam0,omb,h,s8,n,omNu);
  
  cout << popflag_in << " " << xin << " " << lyaxray_in << endl;
  Astrophysics a(&c,popflag_in,xin,lyaxray_in,1.0);
  TwentyOneCM tocm(&c,&a);
  a.initAstrophysics(fstar,fesc,nion,fx,flya,popflag,xrayflag,lyaxrayflag, false);

///////////////////////////////////////////////////////////////
// Calculate the global 21 cm signature
///////////////////////////////////////////////////////////////

  double z;
  double *result;
  double tk,lyaflux,xi,xe;
  double ts,tb;
  double tcmb;

  result=dvector(1,3);

  file="xc_history.dat";
  fout.open(file);

  z=40.0;
  while(z>4.0){
    a.getTIGM(z,result);
    tk=result[1];
    xi=result[2];
    xe=result[3];
    lyaflux=a.lyaFlux(z);
    ts=tocm.tSpinGen(z,tk,xe,lyaflux);
    tb=(1.0-xi)*tocm.tBrightGen(z,tk,xe,lyaflux);
    tcmb=a.getTCmb(z);
    z-=0.1;
    cout<<z<<"\t"<<xi<<"\t"<<xe<<"\t"<<lyaflux<<"\t"<<tk<<"\t"<<tcmb<<"\t"<<ts<<"\t"<<tb<<endl;
    fout<<z<<"\t"<<xi<<"\t"<<xe<<"\t"<<lyaflux<<"\t"<<tk<<"\t"<<tcmb<<"\t"<<ts<<"\t"<<tb<<endl;
  }

  fout.close();

  free_dvector(result,1,3);
 

  ///////////////////////////////////////////////////////////////////////
*/

















    map<string,double> params;
    CosmoCalc cosmo(params);
   // params = cosmo.give_fiducial_params();
   // Global21cmInterface g21(params); 

    ofstream fout;
    fout.open("output/interp.dat");
    double z = 6.5;
    while (z <=9.5){
        fout << z << " " << cosmo.delta_Tb_bar_G21(z) << " " << cosmo.delta_Tb_bar(z)/1000.0 << endl;
        z += 0.01;
    }
    fout.close();
  //ofstream outfile;
    //outfile.open("run_history.dat", ios::out | ios::app);
    
    //CosmoWrite writer(params);
    //writer.calculate_bessels_exact(1000);
    //writer.calculate_bessels_cubic(1000);
    //CosmoCalc cosmo(params);
    //cosmo.compare(1000, 0.5, 0.5);
    //cout << cosmo.limber2(500,24) << endl;
    //cout << cosmo.corr_Tb_new(100, 0.5, 0.5, 0.4, 0.6) << endl;
    //cout << cosmo.Cl_simplified(100, 0.5, 0.5)<< endl;
  /*  ofstream file;
    file.open("output/limber.dat");
    double q = 0.99;
    int l = 3000;
    for (int i = 0; i < 2000; i++) {
        q += 0.00001;
        file << q << " " << 2*(l+0.5)*(l+0.5)*cosmo.limber(l, 1,q)/3.14159265359<< endl;
    }
    file.close();*/
    //writer.calculate_Ml(1,2000,0.01,0.2);
    //writer.calculate_Nl(1,200,0.01,0.02);
    //writer.calculate_Cl_full(100, 0.5, 0.08, 1, 0.0001);
    //writer.generate_movie_Cl(1, 100, 0.5, 0.001, 1, 0.0001);
    //writer.generate_movie(30);
    //writer.calculate_integrandsimple(30, 1, 1, 100);
    //writer.calculate_integrandlong(30, 1, 1, 100);
    
    //writer.calculate_qdot();
    //writer.calculate_q();
    //writer.calculate_dTb(5, 20, 100);
    /* writer.calculate_integrandMM(198, 0.03, 0.03, 1000000);
   
    writer.calculate_integrandMN(142, 1, 1, 1000000);
    
    writer.calculate_integrandMN(142, 0.3, 0.3, 1000000);
    writer.calculate_integrandMN(142, 0.3, 1, 1000000);
    writer.calculate_integrandMN(142, 0.03, 1, 1000000);
    writer.calculate_integrandMN(142, 1.32, 2.5, 1000000);
    writer.calculate_integrandMN(142, 0.012, 0.012, 1000000);
    writer.calculate_integrandMN(142, 0.003, 0.003, 1000000);
    writer.calculate_integrandMN(142, 0.015, 0.015, 1000000);

    writer.calculate_integrandMM(199, 0.5, 0.5, 1000000);
    */ 

/*
    clock_t t1, t2;
    string Fl_filepath = "output/Fls_noRSD_full.dat"; 
    Fisher fish(params, Fl_filepath);
    t1 = clock();
    double res = fish.F("ombh2", "ombh2");
    t2 = clock();
    
    //outfile << " ##################### " << endl;
    //outfile << "kstep_Cl = 3" << endl;
    //outfile << "lmax = 15" << endl; 
    cout << "Result is " << res << endl;
    //outfile << "Result is " << res << endl;
    float diff ((float)t2 - (float)t1);
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    //outfile << "runtime was " << diff/CLOCKS_PER_SEC << endl;
  */  
    /*     
    ifstream filesimlpe, filelong;
    filesimlpe.open("output/Fls_k4_simple.dat");
    filelong.open("output/Fls_k4_long.dat");
    vector<double> flsimple, flslong;
    double a, b;
    while (filesimlpe >> a >> b)
    {
        flsimple.push_back(b);
    }
    while (filelong >> a >> b)
    {
        flslong.push_back(b);
    }
    ofstream errorfile("output/error_fls_k4.dat");
    double err;

    for (int i = 0; i < flsimple.size(); ++i)
    {
        err = abs((flsimple[i] - flslong[i]) / flsimple[i]);
        errorfile << i << " " << err << endl;
    }

    */
    (void) argc;
    (void) argv;
    /*
    mat A = randu<mat>(2,2);
    mat B = randu<mat>(2,2);

    A(0,0) = 1;
    A(0,1) = 0.1;
    A(1,0) = 0;
    A(1,1) = 1.1;
    cout << A << endl;
    cout << B << endl;
    cout << A*B  << endl;
    */
    /*
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
       cout << "New integration yields " << integrate_simps(g, 0.0, 1.0, 100) << endl;
       map<string,double> params;

       params["O"] = 0.02;
       params["ombh2"] = 0.02;
       cout << params["O"] << endl;
       CosmoBasis base(params);
    //base.show_params(); 
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

    //CosmoCalc calc(params);
    //calc.show_cosmo_calcs();
    vector<double> v;
    v.clear();
    cout << v.size()<< endl;

    //CosmoWrite writer(params);
    */
    /* 
       clock_t t1, t2;

       Fisher fish(params);
       t1 = clock();
       double res = fish.Cl_loglog_derivative(142, "ombh2", 0.01, 0.01);
       cout << res << endl;

    //fish.write_logder("ombh2", 0.0226, 0.0001, 0.01, 99, 142, 0.01, 0.01,\
    "_l142_0-01_0-01");
    t2 = clock();
    float diff ((float)t2 - (float)t1);
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    */
    //t1 = clock();
    //writer.calculate_Ml(5, 0.1, 0.01, 0.5, 10000); 
    //writer.calculate_distances(10);
    //t2 = clock();
    //float diff ((float)t2 - (float)t1);
    //cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;
    /* 
       t1 = clock();
    //writer.calculate_densities_rho(5000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    //writer.calculate_densities_Omega(10000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    //writer.calculate_H(1000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_P(0.001, 10, 10000, "default", "default");
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_P_CLASS(0.001, 10, 0, 10000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_dTb(5, 20, 100);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_xHI(0, 20, 100);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_Ts(0, 20, 100);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

    t1 = clock();
    writer.calculate_Tk(0, 1000, 1000);
    t2 = clock();
    diff = (float)t2 - (float)t1;
    cout << "runtime was " << diff/CLOCKS_PER_SEC << endl;

*/

    return 0;
}
