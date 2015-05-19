//--------------------------------------------------------------------------
//
// Description:
// 	class Engine : see header file (Engine.hh) for description.
//
//------------------------------------------------------------------------
//-----------------------
// This Class's Header --
//-----------------------
#include "Engine.hpp"
// ------------------------
// Collaborating classes --
//-------------------------
//--------------------
// C++
//--------------------
#include<numeric>
#include<iostream>
#include<stdexcept>
//--------------------
// C 
//----------------

using namespace std;
//---------------
// Constructors --
//----------------
Engine::Engine():_lmax(-1)
{
}
//--------------
// Destructor --
//--------------

//-----------------
// Member functions --
//-----------------

void 
Engine::writeCls(std::ostream &of){

    vector<unsigned> lvec(_lmax-1,1);
    lvec[0]=2;
    partial_sum(lvec.begin(),lvec.end(),lvec.begin());

    vector<double> cltt,clte,clee,clbb,clpp,cltp,clep;
    bool hasLensing=false;
    try{
        getCls(lvec,cltt,clte,clee,clbb);
        hasLensing=getLensing(lvec,clpp,cltp,clep); 
    }
    catch (std::exception &e){
        cout << "GIOSH" << e.what() << endl;
    }

    cout.precision( 16 );
    for (size_t i=0;i<lvec.size();i++) {
        of << lvec[i] << "\t" 
            << lvec[i]*(lvec[i]+1)*cltt[i] << "\t" 
            << clte[i] << "\t" 
            << clee[i] << "\t" 
            << clbb[i];
        if (hasLensing){
            of << "\t" << clpp[i] << "\t" << cltp[i] << "\t" << clep[i];
        }
        of << "\n";
    }

}

void Engine::writePks(std::ostream &o)
{
    double z = 2;
    double k = 0.0001;
    double kstep = k;
    for (int n = 0; n < 10000; ++n) {
        k += kstep;
        o << k << " " << return_Pkz(k,z) << endl; 
    }
}
