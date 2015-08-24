#pragma once

//This class should include all the functions to test the CosmoCalc class.
//ie. ultimately the best function will be determined and this class will 
//be superfluous... ha... like that is ever going to happen.
// This is mainly so I don't lose the overview of CosmoCalc.

#include "CosmologyCalculatorClass.hpp"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_math.h>


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
            gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);
            double result, error;
            gsl_integration_qag(&F, 7, 9, 0, 1e-5, 1000, 6, w, &result, &error); 
            cout << "error = " << error << endl;
            gsl_integration_workspace_free(w);
            return this->prefactor_Ml * result;
        }
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

        struct my_kappa_integrand_params
        {
            int Pk_index;
            int Tb_index;
            int q_index;
            int l;
            double k1;
            double k2;
            SanityChecker* checker;
        };


        static double Kappa_integrand(double kappa, void *p)
        {
            struct my_kappa_integrand_params * params = (struct my_kappa_integrand_params *)p;
            int l = (params->l);
            int Pk_index = (params->Pk_index);
            int Tb_index = (params->Tb_index);
            int q_index = (params->q_index);
            double k1 = (params->k1);
            double k2 = (params->k2);
            SanityChecker* check = (params->checker);

            if (k1 == k2)
                return pow(kappa,2) * pow(check->M(l,k1,kappa,Pk_index,Tb_index,q_index),2);
            else 
                return pow(kappa,2) * check->M(l,k1,kappa,Pk_index,Tb_index,q_index) *\
                    check->M(l,k2,kappa,Pk_index,Tb_index,q_index);       
        }
        double Cl_gsl(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index)
        {
            double low;
            if (l < 50){
                low = k_low;
            } else if (l < 1000){
                low = (double)l/(1.2*10000);
            } else {
                low = (double)l/(10000);
            }
            double lower_kappa_bound;// = k_low;
            if (low > k_low)
                lower_kappa_bound = low;
            else
                lower_kappa_bound = k_low;

            double higher_kappa_bound = max(k1,k2) + 0.1;
            struct my_kappa_integrand_params params = {Pk_index, Tb_index, q_index, l, k1, k2, this};
            gsl_function F;
            F.function = &Kappa_integrand;
            F.params = &params;
            gsl_integration_workspace *w = gsl_integration_workspace_alloc(100000);
            double result, error;
            gsl_integration_qag(&F, lower_kappa_bound, higher_kappa_bound, 0, 1e-5, 10000, 6, w,\
                    &result, &error); 
            cout << "error = " << error << endl;
            gsl_integration_workspace_free(w);
            return result;

        } 
        struct my_MC_params
        {
            int Pk_index;
            int Tb_index;
            int q_index;
            int l;
            double k1;
            double k2;
            SanityChecker* checker;
        };

        static double g(double *k, size_t dim, void * p)
        {
            struct my_MC_params * params = (struct my_MC_params *)p;
            int l = (params->l);
            int Pk_index = (params->Pk_index);
            int Tb_index = (params->Tb_index);
            int q_index = (params->q_index);
            double k1 = (params->k1);
            double k2 = (params->k2);
            SanityChecker* check = (params->checker);

            double kappa = k[0];
            double z = k[1];
            double zp = k[2];
            double r,q;
            r = check->r_interp(z);
            q = check->q_interp(z, q_index);
            double rp,qp;
            rp = check->r_interp(zp);
            qp = check->q_interp(zp, q_index);

            double hhh = pow(check->qs[q_index].h,3);
            double sP = sqrt(check->Pk_interp(kappa*check->qs[q_index].h,z, Pk_index)/hhh);
            double sPp = sqrt(check->Pk_interp(kappa*check->qs[q_index].h,zp, Pk_index)/hhh);
            double I1 = kappa*kappa * sP * sPp * check->sph_bessel_camb(l,kappa*q) *\
                check->sph_bessel_camb(l,kappa*qp);
            double I2 = rp*rp / (check->Hf_interp(zp)*1000.0) * check->Tb_interp(zp, Tb_index) *\
                check->sph_bessel_camb(l,k2*rp);
            double I3 = r*r / (check->Hf_interp(z)*1000.0) * check->Tb_interp(z, Tb_index) *\
                check->sph_bessel_camb(l,k1*r);

            return pow(check->prefactor_Ml,2) * I1 * I2 * I3;
        }

        double Cl_MC(int l, double k1, double k2, double k_low,\
                double k_high, int Pk_index, int Tb_index, int q_index)
        {
            double low;
            if (l < 50){
                low = k_low;
            } else if (l < 1000){
                low = (double)l/(1.2*10000);
            } else {
                low = (double)l/(10000);
            }
            double lower_kappa_bound;// = k_low;
            if (low > k_low)
                lower_kappa_bound = low;
            else
                lower_kappa_bound = k_low;

            double higher_kappa_bound = max(k1,k2) + 0.1;
            double result, error;
            struct my_MC_params params = {Pk_index, Tb_index, q_index, l, k1, k2, this};
            const gsl_rng_type *T;
            gsl_rng *r;
            double xl[3] = {lower_kappa_bound,7,7};
            double xu[3] = {higher_kappa_bound,9,9};
            gsl_monte_function G  = {&g, 3, &params};

            size_t calls = 100000;
            gsl_rng_env_setup();

            T = gsl_rng_default;
            r = gsl_rng_alloc(T);

            gsl_monte_miser_state *s = gsl_monte_miser_alloc(3);
            gsl_monte_miser_integrate(&G, xl, xu, 3, calls, r, s, &result, &error);
            gsl_monte_miser_free(s);

            return result;

        }
};
