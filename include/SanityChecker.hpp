#pragma once

//This class should include all the functions to test the CosmoCalc class.
//ie. ultimately the best function will be determined and this class will 
//be superfluous... ha... like that is ever going to happen.
// This is mainly so I don't lose the overview of CosmoCalc.

#include "CosmologyCalculatorClass.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

class SanityChecker : public CosmoCalc {

    public:
        SanityChecker(map<string, double> params, int *Pk_index, int *Tb_index, int *q_index);
        ~SanityChecker();

        void kappa_integral(int l, double z, double zp, double *out1, double k_low, double k_high, int n);
        double kappa_integrand(int l, double z, double zp, double kappa);
        void Compare_Cl(int l, double k1, double k2, double k_low,\
                double k_high, int n_levin, double *ratio, double *time_r);
        void plot_integrand_z(int l, double k1, double k2, int zsteps, string filename);
        void plot_intjj(int l, double zp, int zsteps, string filename);
        void plot_integrand_zp(int l, double z, double k2, int zsteps, string filename);

        double Cl_gauss(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index);
        double Cl_gauss_fewerZ(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index);

        double corr_Tb_MM(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index);
        double corr_Tb_MN(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index);

        static double M_integrand_gsl(double z, void * p)
        {
            struct my_M_integrand_params * params = (struct my_M_integrand_params *)p;
            int l = (params->l);
            int Pk_index = (params->Pk_index);
            int Tb_index = (params->Tb_index);
            int q_index = (params->q_index);
            double k1 = (params->k1);
            double kappa = (params->kappa);
            SanityChecker* check = (params->checker);
            double r,q;
            r = check->r_interp(z);
            q = check->q_interp(z, q_index);
            //double res = log(1.0*z)/sqrt(z);
            return pow(r,2)* check->Tb_interp(z, Tb_index) * check->sph_bessel_camb(l,k1*r) *\
                check->sph_bessel_camb(l,kappa*q) *\
                sqrt(check->Pk_interp(kappa*check->qs[q_index].h,z, Pk_index)/\
                        pow(check->qs[q_index].h,3)) / (check->Hf_interp(z)*1000.0);
        }

        double M_gsl(int l, double k1, double kappa, int Pk_index, int Tb_index, int q_index)
        {

            struct my_M_integrand_params params = {Pk_index, Tb_index, q_index, l, k1, kappa, this};
            gsl_function F;
            F.function = &M_integrand_gsl;
            F.params = &params;
            gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(100000);
            double result, error;
            size_t steps;
           // gsl_set_error_handler_off();
            gsl_integration_cquad(&F, 7, 9, 0, 1e-6, w, &result, &error, &steps); 
            cout << "error = " << error << ", steps = " << steps << endl;
            gsl_integration_cquad_workspace_free(w);
            return result;
        }
        double Cl_gsl(double k1, double k2);
        struct my_M_integrand_params
        {
            int Pk_index;
            int Tb_index;
            int q_index;
            int l;
            double k1;
            double kappa;
            SanityChecker* checker;
        };

};
