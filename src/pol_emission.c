/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "gsl/gsl_sf_hyperg.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_errno.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

// FUNCTIONS
////////////

// B.15:
double f_m(double X) {
    return 2.011 * exp(-19.78*pow(X, -0.5175)) -
           cos(1./sqrt(X)*39.89) * exp(-pow(X, -0.6) * 70.16) - 0.011 * exp(-1.69*pow(X, -0.5)) +
           (0.011 * exp(-1.69*pow(X, -0.5)) - 0.003135 * pow(X, 4. / 3.)) *
               0.5 * (1. + tanh(10. * log(0.6648*pow(X, -0.5))));
}

// Bessel function approximations:
double bessel_appr(int n, double x) {

    // use taylor expanded version
    if (x < 1. / 5.) {
        if (n == 0)
            return -log(x / 2.) - 0.5772;

        if (n == 1)
            return 1. / x;

        if (n == 2)
            return 2. / x / x;
    }
    // in this case all bessel functions are really small... Theta_e is small,
    // so no emission anyway.?
    else if (x > (1 / 0.00004) && 0) {
        return 1e-100;
    } else // use full versions in between.
        return gsl_sf_bessel_Kn(n, x);

    exit(0);
}

// Planck function
double planck_function(double nu, double THETA_e) {
    double T = THETA_e * ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT /
               BOLTZMANN_CONSTANT;
    return 2. * PLANCK_CONSTANT * nu * nu * nu /
           (SPEED_OF_LIGHT * SPEED_OF_LIGHT) * 1. /
           (exp(PLANCK_CONSTANT * nu / (BOLTZMANN_CONSTANT * T)) - 1.);
}

double DeltaJ_5(double X) {
    return 0.4379 * log(1. + 1.3414 * pow(X, -0.7517));
}

double get_efficiency(double sigma, double beta){
     double eff;
#if (Epsilon == turbulence)
// Ref: Meringolo 2023ApJ Microphysical Plasma Relations from Special-relativistic Turbulence
    if (sigma > 0.1  &&  beta > 1.e-4 && sigma<10.0 && beta < 2.){
        double e0 = 1.;
        double e1 = - 0.23;
        double e2 = 0.5;
        double e3 = - 10.18;
        eff = e0 + e1/sqrt(sigma) + e2*pow(sigma,(1.0/10.0))*tanh(e3*beta*pow(sigma,(1.0/10.0)));
        }
    else eff = 0.;
    if (eff > 1.) eff = 1.;
    if (eff < 0.) eff = 0.;
    return eff;
#elif (Epsilon == reconnection)
// Ref: Ball 2018ApJ Electron and Proton Acceleration in Trans-relativistic Magnetic Reconnection: Dependence on Plasma Beta and Magnetization
     if (sigma > 0.01 && sigma<7.2 && beta < 2.5 && beta>1.e-4){
         double e0 = 1.;
         double e1 = - 1;
         double e2 = 0.64;
         double e3 = - 68;
         eff = e0 + e1/(4.2*pow(sigma, 0.55) + 1.) + e2*pow(sigma, 0.07)*tanh(e3*beta*pow(sigma,(0.13)));
         }
     else eff = 0.;
     if (eff > 1.) eff = 1.;
     if (eff < 0.)  eff = 0.;
     return eff;
#endif
 }

// calculate fraction of nonthermal electrons DF=VAR_POWER
double get_eta(double theta_e,double B,double eff,double beta,double sigma,
        double *power,double *gamma_min, double *gamma_max,double *gamma_br,double *ratio){
    double eta;
    double aa,bb,cc;
#if (Epsilon == turbulence)
//    aa = 1.8 + 0.2/sqrt(sigma);
//    bb = 1.6 * pow(sigma, (-6./10.));
//    cc = 2.25* pow(sigma,(1./3.));
//    *power = aa + bb*tanh(cc*beta);
//    if (*power < 2.1) *power = 2.1;
//    if (*power > 5.8) *power = 5.8;
    *power=Power;
#elif (Epsilon == reconnection)
    aa = 1.8 + 0.7/sqrt(sigma);
    bb = 3.7*pow(sigma,-0.19);
    cc = 23.4*pow(sigma,0.26);
    *power  = aa + bb*tanh(cc*beta);
    if (*power < 2.1) *power = 2.1;
    if (*power > 4.6) *power = 4.6;
    //*power=Power;
#endif
    double f_theta_e=(6+15*theta_e)/(4+5*theta_e);
    *gamma_min=1+f_theta_e*theta_e;
    *gamma_max=Gamma_max;
    *gamma_br=3.9e3*pow(B/100.0, -2.0);
    *ratio=(pow(*gamma_min,1.-*power)-pow(*gamma_br,1.-*power))/(pow(*gamma_min,1.-*power)-pow(*gamma_br,1.-*power)+
               (1.-1./ *power)*(pow(*gamma_br,-*power)-pow(*gamma_max,-*power)));
    if (*ratio>0.99) *ratio=1.;
    if (*ratio<0.01) *ratio=0.;

    double C1,A;
    if (*gamma_br > *gamma_max*0.9 || Power_type == nth){
        C1=1./(pow(*gamma_min,1.-*power)-pow(*gamma_max,1.-*power));
        A=C1*(*power-1.)*(pow(*gamma_min,2.-*power)-pow(*gamma_max,2.-*power))/(*power-2.);
        eta=eff*(*gamma_min-2.)/((1.-eff)*(A-1.)+eff*(*gamma_min-2.));
        }
    else if ( Power_type == cooling ){
        C1=1/(pow(*gamma_min,1.-*power)-pow(*gamma_br,1.-*power)+
             (1.-1./ *power)*(*gamma_br)*(pow(*gamma_br,-*power)-pow(*gamma_max,-*power)));
        A=C1*(*gamma_br*(pow(*gamma_br,1.-*power)-pow(*gamma_max,1.-*power))+
             (*power-1.)*(pow(*gamma_min,2.-*power)-pow(*gamma_br,2.-*power))/(*power-2.));
        eta=eff*(*gamma_min-2.)/((1.-eff)*(A-1.)+eff*(*gamma_min-2.));
        }

    if(eta<0. || *gamma_min>*gamma_br || C1<0. || A<0. || eta>1.){
        printf("nonthermal partition err: eff=%.2e min=%.2e br=%.2e C1=%.2e A=%.2e eta=%e ration=%e\n",
               eff,*gamma_min,*gamma_br,C1,A,eta,*ratio);
        exit(1);
    }
    return eta;
}

double get_w_kappa(double theta_e, double beta, double sigma, double kappa) {

    if (DF == KAPPA)
    {
        return (kappa-3)*theta_e/kappa;
    }
    double eff, w;
    eff = get_efficiency(sigma, beta);
    if (kappa < 3.){
        w = theta_e + (kappa - 3.)*(MPoME)*sigma*eff/(6.*kappa);
    }
    else{
        // Jordy's prescription
        //w = theta_e*(kappa - 3.)/kappa + (kappa - 3.)*(MPoME)*sigma*eff/(6.*kappa);

        // Hongxuan's prescription
        //w=theta_e*(kappa-3.)/(kappa*(1.-eff));

        // Xufan's prescription
        if ((kappa<(2+1.529/(0.529+eff)))) return theta_e*(kappa-3.)/(kappa*(1.-eff));//return to Hongxuan's
        //Calculate A
        double f_theta_e=(6.+15.*theta_e)/(4.+5.*theta_e);
        double gamma_min=1+f_theta_e*theta_e;
        double A=1./(1.-eff)*(3.*theta_e+pow(gamma_min,3.)/(pow(gamma_min,2.)+2*gamma_min*theta_e+2*pow(theta_e,2.)));
        //3rd order equation
        double a=6.*pow(kappa,2.);
        double b=6.*gamma_min*pow(kappa,2.)-2.*kappa*A*(kappa-3.);
        double c=3.*pow(gamma_min,2.)*(kappa-1.)*kappa-2.*A*gamma_min*kappa*(kappa-3.);
        double d=pow(gamma_min,3.)*(kappa-2.)*(kappa-1.)-A*pow(gamma_min,2.)*(kappa-3.)*(kappa-1.);

        double x0, x1, x2;
        int num_roots = gsl_poly_solve_cubic(b/a, c/a, d/a, &x0, &x1, &x2);

        //printf("num_roots=%d x0=%g x1=%g x2=%g\n",num_roots,x0,x1,x2);
        if (num_roots==1) w=x0;
        else if (num_roots==3){
            w=x0;
            if (x1>w) w=x1;
            if (x2>w) w=x2;
        }
		if (w<0){
			//printf("negative w\n");
			w=theta_e*(kappa-3.)/(kappa*(1.-eff));//return to Hongxuan's
			//exit(0);
		}
    }
    return w;
}

double get_kappa(double beta, double sigma){
    double aa, bb, cc, pp;
#if (DF == KAPPA)
	return kappa_const;
#elif (Epsilon == turbulence)
// Ref: Meringolo 2023ApJ Microphysical Plasma Relations from Special-relativistic Turbulence
	if (sigma > 0.1  &&  beta > 1.e-4 && sigma<10.0 && beta < 2.){
		aa = 2.8 + 0.2/sqrt(sigma);
		bb = 1.6 * pow(sigma, (-6./10.));
		cc = 2.25* pow(sigma,(1./3.));
		pp = aa + bb*tanh(cc*beta);
		if (pp < 3.1) pp = 3.1;
		if (pp > 7.5) pp = 7.5;
	}
	else pp=7.5;
    return pp;
#elif (Epsilon == reconnection)
// Ref: Ball 2018ApJ Electron and Proton Acceleration in Trans-relativistic Magnetic Reconnection: Dependence on
	if (sigma > 0.01 && sigma<7.2 && beta < 2.5 && beta>1.e-4){
    aa = 2.8 + 0.7/sqrt(sigma);
    bb = 3.7 * pow(sigma, (-0.19));
    cc = 23.4* pow(sigma,(0.26));
    pp = aa + bb*tanh(cc*beta);
    if (pp < 3.1) pp = 3.1;
    if (pp > 7.5) pp = 7.5;
	}
	else pp=7.5;
    return pp;
#endif
}


///////////////////////////Rho_Q

double rho_Q_kappa(double theta_e, double n_e, double nu, double B,
                   double theta_B, double beta, double sigma) {
    double kappa = get_kappa(beta, sigma);
    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double factor_F = 0, factor_Q = 0;

    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));

    double F_factor_K3_5 =
        (1 - exp(-pow(X_kappa, 0.84) / 30.) -
         sin(X_kappa / 10.) * exp(-3 * pow(X_kappa, 0.471) / 2.));
    double F_factor_K4 = (1 - exp(-pow(X_kappa, 0.84) / 18.) -
                          sin(X_kappa / 6.) * exp(-7 * pow(X_kappa, 0.5) / 4.));
    double F_factor_K4_5 = (1 - exp(-pow(X_kappa, 0.84) / 12.) -
                            sin(X_kappa / 4.) * exp(-2 * pow(X_kappa, 0.525)));
    double F_factor_K5 =
        (1 - exp(-pow(X_kappa, 0.84) / 8.) -
         sin(3 * X_kappa / 8.) * exp(-9 * pow(X_kappa, 0.541) / 4.));

    double Q_factor_K3_5 =
        (17. * w - 3 * pow(w, 1. / 2.) + 7 * pow(w, 1. / 2.) * exp(-5. * w));
    double Q_factor_K4 = ((46. / 3.) * w - (5. / 3.) * pow(w, 1. / 2.) +
                          (17. / 3.) * pow(w, 1. / 2.) * exp(-5. * w));
    double Q_factor_K4_5 = (14. * w - (13. / 8.) * pow(w, 1. / 2.) +
                            (9. / 2.) * pow(w, 1. / 2.) * exp(-5. * w));
    double Q_factor_K5 =
        ((25. / 2.) * w - pow(w, 1. / 2.) + 5 * pow(w, 1. / 2.) * exp(-5. * w));

    if (kappa < 3.5)
    {
        factor_F = F_factor_K3_5;
        factor_Q = Q_factor_K3_5;
    }
    else if (kappa >= 3.5 && kappa < 4.0){
        factor_F = (F_factor_K4 - F_factor_K3_5) / 0.5 * (kappa - 3.5) + F_factor_K3_5;
        factor_Q = (Q_factor_K4 - Q_factor_K3_5) / 0.5 * (kappa - 3.5) + Q_factor_K3_5;
    }
    else if (kappa >= 4.0 && kappa < 4.5){
        factor_F = (F_factor_K4_5 - F_factor_K4) / 0.5 * (kappa - 4.) + F_factor_K4;
        factor_Q = (Q_factor_K4_5 - Q_factor_K4) / 0.5 * (kappa - 4.) + Q_factor_K4;
    }
    else if (kappa >= 4.5 && kappa < 5.0){
        factor_F = (F_factor_K5 - F_factor_K4_5) / 0.5 * (kappa - 4.5) + F_factor_K4_5;
        factor_Q = (Q_factor_K5 - Q_factor_K4_5) / 0.5 * (kappa - 4.5) + Q_factor_K4_5;
    }
    else if (kappa >= 5.0){
        factor_F = F_factor_K5;
        factor_Q = Q_factor_K5;
    }

    return -n_e * pow(ELECTRON_CHARGE, 2) * pow(nuc, 2) * pow(sin(theta_B), 2) /
           (ELECTRON_MASS * SPEED_OF_LIGHT * pow(nu, 3)) * factor_F * factor_Q;
}


double rho_Q_thermal(double theta_e, double n_e, double nu, double B,
                     double theta_B) {
    double nu_c = ELECTRON_CHARGE * B / (2.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double nu_s = (2. / 9.)*nu_c * sin(theta_B) * theta_e * theta_e;
    double X = nu / nu_s;

    double Thetaer = 1. / theta_e;

    return - n_e*ELECTRON_CHARGE*ELECTRON_CHARGE*nu_c*nu_c*sin(theta_B) * sin(theta_B)/\
            (ELECTRON_MASS*SPEED_OF_LIGHT*nu*nu*nu)*f_m(X)*(bessel_appr(1, Thetaer) /\
             bessel_appr(2, Thetaer) + 6. * theta_e);
}

double rho_Q_power(double theta_e, double n_e, double nu, double B,double theta_B,
                   double power,double gamma_min,double gamma_max) {
    double nuc =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double P_perp =
        (n_e * ELECTRON_CHARGE * ELECTRON_CHARGE) /
        (ELECTRON_MASS * SPEED_OF_LIGHT * nuc * sin(theta_B)) * (power - 1) *
        (1. / (pow(gamma_min, 1. - power) - pow(gamma_max, 1. - power)));
    return -P_perp * pow((nuc * sin(theta_B)) / nu, 3.) *
           (pow(gamma_min, 2. - power) / ((power / 2.) - 1.)) *
           (1. -
            pow((2. * nuc * sin(theta_B) * gamma_min * gamma_min) / (3. * nu),
                (power / 2.) - 1.));
}

double rho_Q(double theta_e, double n_e, double nu, double B, double theta_B, double beta, double sigma) {
#if (DF == KAPPA)
    return rho_Q_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
#elif (DF == VAR_KAPPA)
    double eff;
    eff = get_efficiency(sigma, beta);
    double rQ_kappa = rho_Q_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
    if (isnan(rQ_kappa))
    {
        rQ_kappa = 0;
    }
    double nu_break = 2.5e30/(B/100.);
    if (nu > nu_break)
    {
        return ((1. - eff) * rho_Q_thermal(theta_e, n_e, nu, B, theta_B) + eff * rQ_kappa)/sqrt(nu/nu_break);
    }
    else{
        return ((1. - eff) * rho_Q_thermal(theta_e, n_e, nu, B, theta_B) + eff * rQ_kappa);
    }
#elif (DF == POWER)
    return rho_Q_power(theta_e, n_e, nu, B, theta_B,Power,Gamma_min,Gamma_max);
#elif (DF == TH)
    return rho_Q_thermal(theta_e, n_e, nu, B, theta_B);
#elif (DF == VAR_POWER)
    double eff = get_efficiency(sigma, beta);
    if (eff<0.01) return rho_Q_thermal(theta_e, n_e, nu, B, theta_B);
    else{
        double power,gamma_min,gamma_max,gamma_br,eta,ratio;
        eta=get_eta(theta_e,B,eff,beta,sigma,&power,&gamma_min,&gamma_max,&gamma_br,&ratio);
        //if(eta>0.09)printf("eta =%e\n",eta);
        if (Power_type == nth || gamma_br>gamma_max*0.9)
        return (rho_Q_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                rho_Q_power(theta_e, n_e*eta, nu, B, theta_B,power,gamma_min,gamma_max));
        else if (Power_type == cooling)
        return (rho_Q_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                rho_Q_power(theta_e, n_e*eta*ratio, nu, B, theta_B,power,gamma_min,gamma_br)+
                rho_Q_power(theta_e, n_e*eta*(1-ratio), nu, B, theta_B,power+1.,gamma_br,gamma_max));
        }
#endif
}

///////////Rho_V

double rho_V_kappa(double theta_e, double n_e, double nu, double B,
                   double theta_B, double beta, double sigma) {
    double kappa = get_kappa(beta, sigma);

    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double factor_G = 0, factor_V = 0;
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2.) * sin(theta_B));
    double G_factor_K3_5 =
        (1. - 0.17 * log(1 + 0.447 * pow(X_kappa, -1. / 2.)));
    double G_factor_K4 = (1. - 0.17 * log(1 + 0.391 * pow(X_kappa, -1. / 2.)));
    double G_factor_K4_5 =
        (1. - 0.17 * log(1 + 0.348 * pow(X_kappa, -1. / 2.)));
    double G_factor_K5 = (1. - 0.17 * log(1 + 0.313 * pow(X_kappa, -1. / 2.)));

    double V_factor_K3_5 =
        ((pow(w, 2.) + 2. * w + 1.) / ((25. / 8.) * pow(w, 2.) + 4. * w + 1.));
    double V_factor_K4 = ((pow(w, 2.) + 54 * w + 50.) /
                          ((30. / 11.) * pow(w, 2) + 134 * w + 50));
    double V_factor_K4_5 = ((pow(w, 2.) + 43 * w + 38.) /
                            ((7. / 3.) * pow(w, 2) + (185. / 2.) * w + 38));
    double V_factor_K5 = ((w + (13. / 14.)) / (2 * w + (13. / 14.)));

    if (kappa < 3.5)
    {
        factor_G = G_factor_K3_5;
        factor_V = V_factor_K3_5;
    }
    else if (kappa >= 3.5 && kappa < 4.0){
        factor_G = (G_factor_K4 - G_factor_K3_5) / 0.5 * (kappa - 3.5) + G_factor_K3_5;
        factor_V = (V_factor_K4 - V_factor_K3_5) / 0.5 * (kappa - 3.5) + V_factor_K3_5;
    }
    else if (kappa >= 4.0 && kappa < 4.5){
        factor_G = (G_factor_K4_5 - G_factor_K4) / 0.5 * (kappa - 4.) + G_factor_K4;
        factor_V = (V_factor_K4_5 - V_factor_K4) / 0.5 * (kappa - 4.) + V_factor_K4;
    }
    else if (kappa >= 4.5 && kappa < 5.0){
        factor_G = (G_factor_K5 - G_factor_K4_5) / 0.5 * (kappa - 4.5) + G_factor_K4_5;
        factor_V = (V_factor_K5 - V_factor_K4_5) / 0.5 * (kappa - 4.5) + V_factor_K4_5;
    }
    else if (kappa >= 5.0){
        factor_G = G_factor_K5;
        factor_V = V_factor_K5;
    }

    return 2 * n_e * pow(ELECTRON_CHARGE, 2.) * nuc * cos(theta_B) /
           (ELECTRON_MASS * SPEED_OF_LIGHT * pow(nu, 2.)) *
           bessel_appr(0., 1. / w) / (bessel_appr(2., 1. / w)) * factor_G *
           factor_V;
}

double rho_V_thermal(double theta_e, double n_e, double nu, double B,
                     double theta_B) {
    double nu_c = ELECTRON_CHARGE * B / (2.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double nu_s = (2. / 9.)*nu_c * sin(theta_B) * theta_e * theta_e;

    double X = nu / nu_s;
    double Thetaer = 1. / theta_e;

    double k2 = bessel_appr(2, Thetaer);
    double k0 = bessel_appr(0, Thetaer);

    return 2.0 * n_e *ELECTRON_CHARGE*ELECTRON_CHARGE*nu_c / (ELECTRON_MASS*SPEED_OF_LIGHT*nu*nu) * \
            (k0 - DeltaJ_5(X))/k2* cos(theta_B);
}

double rho_V_power(double theta_e, double n_e, double nu, double B,double theta_B,
                   double power,double gamma_min,double gamma_max) {
    double nuc =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double P_perp =
        (n_e * ELECTRON_CHARGE * ELECTRON_CHARGE) /
        (ELECTRON_MASS * SPEED_OF_LIGHT * nuc * sin(theta_B)) * (power - 1.) *
        pow((pow(gamma_min, 1. - power) - pow(gamma_max, 1. - power)), -1.);
    return 2. * P_perp * ((power + 2.) / (power + 1.)) *
           pow(((nuc * sin(theta_B)) / (nu)), 2.) *
           pow(gamma_min, -(power + 1.)) * log(gamma_min) * cos(theta_B) /
           sin(theta_B);
}
double rho_V(double theta_e, double n_e, double nu, double B, double theta_B, double beta, double sigma) {
#if (DF == KAPPA)
    return rho_V_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
#elif (DF == VAR_KAPPA)
    double eff;
    eff = get_efficiency(sigma, beta);
    double rV_kappa = rho_V_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
    if (isnan(rV_kappa))
    {
        rV_kappa = 0;
    }
    double nu_break = 2.5e30/(B/100.);
    if (nu > nu_break)
    {
        return ((1. - eff) * rho_V_thermal(theta_e, n_e, nu, B, theta_B) + eff * rV_kappa)/sqrt(nu/nu_break);
    }
    else{
        return ((1. - eff) * rho_V_thermal(theta_e, n_e, nu, B, theta_B) + eff * rV_kappa);
    }

#elif (DF == POWER)
    return rho_V_power(theta_e, n_e, nu, B, theta_B,Power,Gamma_min,Gamma_max);
#elif (DF == TH)
    return rho_V_thermal(theta_e, n_e, nu, B, theta_B);
#elif (DF == VAR_POWER)
    double eff = get_efficiency(sigma, beta);
    if (eff<0.01) return rho_V_thermal(theta_e, n_e, nu, B, theta_B);
    else{
        double power,gamma_min,gamma_max,gamma_br,eta,ratio;
        eta=get_eta(theta_e,B,eff,beta,sigma,&power,&gamma_min,&gamma_max,&gamma_br,&ratio);
        if (Power_type == nth || gamma_br>gamma_max*0.9)
        return (rho_V_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                rho_V_power(theta_e, n_e*eta, nu, B, theta_B,power,gamma_min,gamma_max));
        else if (Power_type == cooling)
        return (rho_V_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                rho_V_power(theta_e, n_e*eta*ratio, nu, B, theta_B,power,gamma_min,gamma_br)+
                rho_V_power(theta_e, n_e*eta*(1-ratio), nu, B, theta_B,power+1.,gamma_br,gamma_max));
        }
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////J
/// for thermal and kappa

double j_I_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B, double beta, double sigma) {

    double kappa = get_kappa(beta, sigma);
    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double J_low_factor_I = 1;
    double J_high_factor_I = 1;
    double J_S_factor_I = 1;
    double J_x_I = 3 * pow(kappa, -3. / 2.);
    double J_I_low_kappa = pow(X_kappa, 1. / 3.) * sin(theta_B) * 4 * M_PI *
                           tgamma(kappa - 4. / 3.) /
                           (pow(3, 7. / 3.) * tgamma(kappa - 2.)) *
                           J_low_factor_I;
    double J_I_high_kappa = pow(X_kappa, -(kappa - 2) / 2.) * sin(theta_B) *
                            pow(3, (kappa - 1) / 2.) * (kappa - 2.) *
                            (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
                            tgamma(kappa / 4. + 4. / 3.) * J_high_factor_I;
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac *
           pow((pow(J_I_low_kappa, -J_x_I) + pow(J_I_high_kappa, -J_x_I)),
               -1. / J_x_I) *
           J_S_factor_I;
}

double j_I_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B) {

    double nu_c = ELECTRON_CHARGE * B / (2.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double nu_s = (2. / 9.)*nu_c * sin(theta_B) * theta_e * theta_e;

    double x = nu / nu_s;

    double j_I = (n_e*ELECTRON_CHARGE*ELECTRON_CHARGE*nu_c)/SPEED_OF_LIGHT*exp( - pow(x, 1. / 3.)) * sqrt(2.) * M_PI / 27. * sin(theta_B) * \
           pow((sqrt(x) + pow(2., 11./12.) * pow(x, 1./6.)), 2.);

    return j_I;
}

double j_I_power(double theta_e, double n_e, double nu, double B,double theta_B,
                 double power,double gamma_min,double gamma_max) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double j_I_power_factor = 1;
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac * (pow(3, power / 2.) * (power - 1) * sin(theta_B)) /
           (2 * (power + 1) *
            (pow(gamma_min, 1 - power) - pow(gamma_max, 1 - power))) *
           tgamma((3 * power - 1) / 12.) * tgamma((3 * power + 19) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power - 1) / 2.) *
           j_I_power_factor;
}

double j_I(double theta_e, double n_e, double nu, double B, double theta_B, double beta, double sigma) {

#if (DF == KAPPA)
    return j_I_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
#elif (DF == VAR_KAPPA)
    double eff;
    eff = get_efficiency(sigma, beta);
    double jI_kappa = j_I_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
    if (isnan(jI_kappa))
    {
        jI_kappa = 0;
    }
    double nu_break = 2.5e30/(B/100.);
    if (nu > nu_break)
    {
        return ((1. - eff) * j_I_thermal(theta_e, n_e, nu, B, theta_B) + eff * jI_kappa)/sqrt(nu/nu_break);
    }
    else{
        return ((1. - eff) * j_I_thermal(theta_e, n_e, nu, B, theta_B) + eff * jI_kappa);
    }

#elif (DF == POWER)
    return j_I_power(theta_e, n_e, nu, B, theta_B,Power,Gamma_min,Gamma_max);
#elif (DF == TH)
    return j_I_thermal(theta_e, n_e, nu, B, theta_B);
#elif (DF == VAR_POWER)
    double eff = get_efficiency(sigma, beta);
    if (eff<0.01) return j_I_thermal(theta_e, n_e, nu, B, theta_B);
    else{
        double power,gamma_min,gamma_max,gamma_br,eta,ratio;
        eta=get_eta(theta_e,B,eff,beta,sigma,&power,&gamma_min,&gamma_max,&gamma_br,&ratio);
        if (Power_type == nth || gamma_br>gamma_max*0.9)
        return (j_I_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                j_I_power(theta_e, n_e*eta, nu, B, theta_B,power,gamma_min,gamma_max));
        else if (Power_type == cooling)
        return (j_I_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                j_I_power(theta_e, n_e*eta*ratio, nu, B, theta_B,power,gamma_min,gamma_br)+
                j_I_power(theta_e, n_e*eta*(1-ratio), nu, B, theta_B,power+1.,gamma_br,gamma_max));
        }
#endif
}

double j_Q_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B, double beta, double sigma) {
    double kappa = get_kappa(beta, sigma);

    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double J_low_factor_Q = 0.5;
    double J_high_factor_Q = (pow((4. / 5.), 2.) + kappa / 50.);
    double J_S_factor_Q = (-1);
    double J_x_Q = (37. / 10.) * pow(kappa, -8. / 5.);
    double J_Q_low_kappa = pow(X_kappa, 1. / 3.) * sin(theta_B) * 4 * M_PI *
                           tgamma(kappa - 4. / 3.) /
                           (pow(3, 7. / 3.) * tgamma(kappa - 2.)) *
                           J_low_factor_Q;
    double J_Q_high_kappa = pow(X_kappa, -(kappa - 2) / 2.) * sin(theta_B) *
                            pow(3, (kappa - 1) / 2.) * (kappa - 2.) *
                            (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
                            tgamma(kappa / 4. + 4. / 3.) * J_high_factor_Q;
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac *
           pow((pow(J_Q_low_kappa, -J_x_Q) + pow(J_Q_high_kappa, -J_x_Q)),
               -1. / J_x_Q) *
           J_S_factor_Q;
}

double j_Q_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B) {

    double nu_c = ELECTRON_CHARGE * B / (2.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double nu_s = (2. / 9.)*nu_c * sin(theta_B) * theta_e * theta_e;

    double x = nu / nu_s;

    return (n_e*ELECTRON_CHARGE*ELECTRON_CHARGE*nu_c)/SPEED_OF_LIGHT*exp(-pow(x, 1./3.))*\
           (-sqrt(2)*M_PI/27.)*sin(theta_B)*pow((sqrt(x) + (7.*pow(theta_e, 24./25.) + 35.)/\
            (10.*pow(theta_e, 24./25.) + 75.)*pow(2., 11./12.)*pow(x, 1./6.)),2.);
}

double j_Q_power(double theta_e, double n_e, double nu, double B,double theta_B,
                 double power,double gamma_min,double gamma_max) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double j_Q_power_factor = -(power + 1) / (power + 7. / 3.);
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac * (pow(3, power / 2.) * (power - 1) * sin(theta_B)) /
           (2 * (power + 1) *
            (pow(gamma_min, 1 - power) - pow(gamma_max, 1 - power))) *
           tgamma((3 * power - 1) / 12.) * tgamma((3 * power + 19) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power - 1) / 2.) *
           j_Q_power_factor;
}

double j_Q(double theta_e, double n_e, double nu, double B, double theta_B, double beta, double sigma) {
#if (DF == KAPPA)
    return j_Q_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
#elif (DF == VAR_KAPPA)
    double eff;
    eff = get_efficiency(sigma, beta);
    double jQ_kappa = j_Q_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
    if (isnan(jQ_kappa))
    {
        jQ_kappa = 0;
    }
    double nu_break = 2.5e30/(B/100.);
    if (nu > nu_break)
    {
        return ((1. - eff) * j_Q_thermal(theta_e, n_e, nu, B, theta_B) + eff * jQ_kappa)/sqrt(nu/nu_break);
    }
    else{
        return ((1. - eff) * j_Q_thermal(theta_e, n_e, nu, B, theta_B) + eff * jQ_kappa);
    }

#elif (DF == POWER)
    return j_Q_power(theta_e, n_e, nu, B, theta_B,Power,Gamma_min,Gamma_max);
#elif (DF == TH)
    return j_Q_thermal(theta_e, n_e, nu, B, theta_B);
#elif (DF == VAR_POWER)
    double eff = get_efficiency(sigma, beta);
    if (eff<0.01) return j_Q_thermal(theta_e, n_e, nu, B, theta_B);
    else{
        double power,gamma_min,gamma_max,gamma_br,eta,ratio;
        eta=get_eta(theta_e,B,eff,beta,sigma,&power,&gamma_min,&gamma_max,&gamma_br,&ratio);
        if (Power_type == nth || gamma_br>gamma_max*0.9)
        return (j_Q_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                j_Q_power(theta_e, n_e*eta, nu, B, theta_B,power,gamma_min,gamma_max));
        else if (Power_type == cooling)
        return (j_Q_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                j_Q_power(theta_e, n_e*eta*ratio, nu, B, theta_B,power,gamma_min,gamma_br)+
                j_Q_power(theta_e, n_e*eta*(1-ratio), nu, B, theta_B,power+1.,gamma_br,gamma_max));
        }
#endif
}

double j_V_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B, double beta, double sigma) {
    double kappa = get_kappa(beta, sigma);

    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double J_low_factor_V =
        (pow(0.75, 2.) * pow((pow(sin(theta_B), -12. / 5.) - 1), 12. / 25.) *
         pow(kappa, -66. / 125.) * pow(X_kappa, -7. / 20.) / w);
    double J_high_factor_V =
        (pow((7. / 8.), 2.) *
         pow((pow(sin(theta_B), -5. / 2.) - 1), 11. / 25.) *
         pow(kappa, -11. / 25.) * pow(X_kappa, -1. / 2.) / w);
    double J_S_factor_V = (sign(cos(theta_B)));
    double J_x_V = 3 * pow(kappa, -3. / 2.);
    double J_V_low_kappa = pow(X_kappa, 1. / 3.) * sin(theta_B) * 4 * M_PI *
                           tgamma(kappa - 4. / 3.) /
                           (pow(3, 7. / 3.) * tgamma(kappa - 2.)) *
                           J_low_factor_V;
    double J_V_high_kappa = pow(X_kappa, -(kappa - 2) / 2.) * sin(theta_B) *
                            pow(3, (kappa - 1) / 2.) * (kappa - 2.) *
                            (kappa - 1.) / 4. * tgamma(kappa / 4. - 1. / 3.) *
                            tgamma(kappa / 4. + 4. / 3.) * J_high_factor_V;
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac *
           pow((pow(J_V_low_kappa, -J_x_V) + pow(J_V_high_kappa, -J_x_V)),
               -1. / J_x_V) *
           J_S_factor_V;
}

double j_V_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B) {

    double nu_c = ELECTRON_CHARGE * B / (2.0 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double nu_s = (2. / 9.)*nu_c * sin(theta_B) * theta_e * theta_e;

    double x = nu / nu_s;

    return (n_e*ELECTRON_CHARGE*ELECTRON_CHARGE*nu_c)/SPEED_OF_LIGHT*exp(-pow(x, 1./3.))\
           *1./theta_e*cos(theta_B)*(M_PI/3. + M_PI/3.*pow(x, 1./3.) + \
            (2./300.)*sqrt(x) + (2*M_PI)/19.*pow(x, 2./3.));
}

double j_V_power(double theta_e, double n_e, double nu, double B,double theta_B,
                 double power,double gamma_min,double gamma_max) {
    double nuc =
        ELECTRON_CHARGE * B / (2. * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double j_V_power_factor = (171. / 250.) *
                              (pow(power, 49. / 100.) / tan(theta_B)) *
                              pow(nu / (3. * nuc * sin(theta_B)), -1. / 2.);
    double prefac =
        n_e * ELECTRON_CHARGE * ELECTRON_CHARGE * nuc / SPEED_OF_LIGHT;
    return prefac * (pow(3., power / 2.) * (power - 1.) * sin(theta_B)) /
           (2. * (power + 1.) *
            (pow(gamma_min, 1. - power) - pow(gamma_max, 1. - power))) *
           tgamma((3. * power - 1.) / 12.) * tgamma((3. * power + 19.) / 12.) *
           pow(nu / (nuc * sin(theta_B)), -(power - 1.) / 2.) *
           j_V_power_factor;
}

double j_V(double theta_e, double n_e, double nu, double B, double theta_B, double beta, double sigma) {
#if (DF == KAPPA)
    return j_V_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
#elif (DF == VAR_KAPPA)
    double eff;
    eff = get_efficiency(sigma, beta);
    double jV_kappa = j_V_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
    if (isnan(jV_kappa))
    {
        jV_kappa = 0;
    }
    double nu_break = 2.5e30/(B/100.);
    if (nu > nu_break)
    {
        return ((1. - eff) * j_V_thermal(theta_e, n_e, nu, B, theta_B) + eff * jV_kappa)/sqrt(nu/nu_break);
    }
    else{
        return ((1. - eff) * j_V_thermal(theta_e, n_e, nu, B, theta_B) + eff * jV_kappa);
    }

#elif (DF == POWER)
    return j_V_power(theta_e, n_e, nu, B, theta_B,Power,Gamma_min,Gamma_max);
#elif (DF == TH)
    return j_V_thermal(theta_e, n_e, nu, B, theta_B);
#elif (DF == VAR_POWER)
    double eff = get_efficiency(sigma, beta);
    if (eff<0.01) return j_V_thermal(theta_e, n_e, nu, B, theta_B);
    else{
        double power,gamma_min,gamma_max,gamma_br,eta,ratio;
        eta=get_eta(theta_e,B,eff,beta,sigma,&power,&gamma_min,&gamma_max,&gamma_br,&ratio);
        if (Power_type == nth || gamma_br>gamma_max*0.9)
        return (j_V_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                j_V_power(theta_e, n_e*eta, nu, B, theta_B,power,gamma_min,gamma_max));
        else if (Power_type == cooling)
        return (j_V_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B)+
                j_V_power(theta_e, n_e*eta*ratio, nu, B, theta_B,power,gamma_min,gamma_br)+
                j_V_power(theta_e, n_e*eta*(1-ratio), nu, B, theta_B,power+1.,gamma_br,gamma_max));
        }
#endif
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////A
/// for thermal and kappa
double hyp2F1_f(double theta_e, double beta, double sigma, double kappa) {
    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double a = kappa - 1. / 3.;
    double b = kappa + 1.;
    double c = kappa + 2. / 3.;
    double X = kappa * w;
    double z;

    if (X < 1e-4) {
        return 0;
    }
    z = -X;

    // 设置GSL的错误处理机制
    gsl_set_error_handler_off(); // 关闭默认的错误处理

    // 使用try-catch的方式检测GSL函数返回的错误
    gsl_sf_result hyp_result1,hyp_result2;
    int status1 = gsl_sf_hyperg_2F1_e(a, c - b, a - b + 1., 1. / (1. - z), &hyp_result1);
    int status2 = gsl_sf_hyperg_2F1_e(b, c - a, b - a + 1., 1. / (1. - z), &hyp_result2);

    // 检查返回状态
    if (status1 != GSL_SUCCESS || status2 != GSL_SUCCESS) {
        printf("GSL hyp2F1 Error: a=%.4f b=%.4f c=%.4f z=%.4f\n",a,b,c,z);
        exit(1); // 返回默认值或根据需求进行其他处理
    }

    // 如果没有错误，返回计算结果
    return pow(1. - z, -a) * tgamma(c) * tgamma(b - a) /
               (tgamma(b) * tgamma(c - a)) * hyp_result1.val +
           pow(1. - z, -b) * tgamma(c) * tgamma(a - b) /
               (tgamma(a) * tgamma(c - b)) * hyp_result2.val;
}



double a_I_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B, double beta, double sigma) {
    double kappa = get_kappa(beta, sigma);
    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double A_low_factor_I = 1;
    double A_high_factor_I = (pow((3.0 / kappa), 19. / 4.) + (3. / 5.));
    double A_S_factor_I = 1;
    double A_x_I = pow((-7. / 4. + kappa * 8. / 5.), -43. / 50.);
    double A_I_low_kappa = pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) *
                           (10. / 41.) * 2 * M_PI /
                           pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
                           (kappa - 1.) * kappa / (3. * kappa - 1.) *
                           tgamma(5. / 3.) * hyp2F1_f(theta_e, beta, sigma, kappa) * A_low_factor_I;
    double A_I_high_kappa =
        pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
        ((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
        (2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) * A_high_factor_I;
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);
    return prefac *
           pow((pow(A_I_low_kappa, -A_x_I) + pow(A_I_high_kappa, -A_x_I)),
               -1. / A_x_I) *
           A_S_factor_I;
}

double a_I_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_I_thermal) {
    double B_nu = planck_function(nu, theta_e); // Planck function
    return j_I_thermal / B_nu;
}

double a_I_power(double theta_e, double n_e, double nu, double B,double theta_B,
                 double power,double gamma_min,double gamma_max) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double a_I_power_factor = 1;
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);
    return prefac * (pow(3, (power + 1) / 2.) * (power - 1)) /
           (4 * (pow(gamma_min, 1 - power) - pow(gamma_max, 1 - power))) *
           tgamma((3 * power + 2) / 12.) * tgamma((3 * power + 22) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power + 2) / 2.) *
           a_I_power_factor;
}

double a_I(double theta_e, double n_e, double nu, double B, double theta_B,
           double beta, double sigma) {
	 double jI_thermal = j_I_thermal(theta_e, n_e, nu, B, theta_B);
#if (DF == KAPPA)
    return a_I_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
#elif (DF == VAR_KAPPA)
    double eff;
    eff = get_efficiency(sigma, beta);
    double aI_kappa = a_I_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
    if (isnan(aI_kappa)) aI_kappa = 0;
    double nu_break = 2.5e30/(B/100.);
    if (nu > nu_break)
    {
        return ((1. - eff) * a_I_thermal(theta_e, n_e, nu, B, theta_B, jI_thermal) + eff * aI_kappa)/sqrt(nu/nu_break);
    }
    else{
        return ((1. - eff) * a_I_thermal(theta_e, n_e, nu, B, theta_B, jI_thermal) + eff * aI_kappa);
    }

#elif (DF == POWER)
    return a_I_power(theta_e, n_e, nu, B, theta_B,Power,Gamma_min,Gamma_max);
#elif (DF == TH)
    return a_I_thermal(theta_e, n_e, nu, B, theta_B, jI_thermal);
#elif (DF == VAR_POWER)
    double eff = get_efficiency(sigma, beta);
    if (eff<0.01) return a_I_thermal(theta_e, n_e, nu, B, theta_B, jI_thermal);
    else{
        double power,gamma_min,gamma_max,gamma_br,eta,ratio;
        eta=get_eta(theta_e,B,eff,beta,sigma,&power,&gamma_min,&gamma_max,&gamma_br,&ratio);
        if (Power_type == nth || gamma_br>gamma_max*0.9)
        return (a_I_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B, jI_thermal)+
                a_I_power(theta_e, n_e*eta, nu, B, theta_B,power,gamma_min,gamma_max));
        else if (Power_type == cooling)
        return (a_I_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B, jI_thermal)+
                a_I_power(theta_e, n_e*eta*ratio, nu, B, theta_B,power,gamma_min,gamma_br)+
                a_I_power(theta_e, n_e*eta*(1-ratio), nu, B, theta_B,power+1.,gamma_br,gamma_max));
        }
#endif
}


double a_Q_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B, double beta, double sigma) {
    double kappa = get_kappa(beta, sigma);

    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double A_low_factor_Q = (25. / 48.);
    double A_high_factor_Q =
        (pow(21, 2.) * pow(kappa, -pow(12. / 5., 2.)) + (11. / 20.));
    double A_S_factor_Q = (-1);
    double A_x_Q = (7. / 5.) * pow(kappa, -23. / 20.);
    double A_Q_low_kappa = pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) *
                           (10. / 41.) * 2 * M_PI /
                           pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
                           (kappa - 1.) * kappa / (3. * kappa - 1.) *
                           tgamma(5. / 3.) * hyp2F1_f(theta_e, beta, sigma, kappa) * A_low_factor_Q;
    double A_Q_high_kappa =
        pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
        ((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
        (2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) * A_high_factor_Q;
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);

    return prefac *
           pow((pow(A_Q_low_kappa, -A_x_Q) + pow(A_Q_high_kappa, -A_x_Q)),
               -1. / A_x_Q) *
           A_S_factor_Q;
}

double a_Q_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_Q_thermal) {
    double B_nu = planck_function(nu, theta_e); // Planck function
    return j_Q_thermal / B_nu;
}

double a_Q_power(double theta_e, double n_e, double nu, double B,double theta_B,
                 double power,double gamma_min,double gamma_max) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double a_Q_power_factor =
        -pow(((17 * power / 500.) - 43. / 1250.), 43. / 500.);
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);

    return prefac * (pow(3, (power + 1) / 2.) * (power - 1)) /
           (4 * (pow(gamma_min, 1 - power) - pow(gamma_max, 1 - power))) *
           tgamma((3 * power + 2) / 12.) * tgamma((3 * power + 22) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power + 2) / 2.) *
           a_Q_power_factor;
}

double a_Q(double theta_e, double n_e, double nu, double B, double theta_B,
           double beta, double sigma) {
	double jQ_thermal = j_Q_thermal(theta_e, n_e, nu, B, theta_B);
#if (DF == KAPPA)
    return a_Q_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
#elif (DF == VAR_KAPPA)
    double eff;
    eff = get_efficiency(sigma, beta);
    double aQ_kappa = a_Q_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
    if (isnan(aQ_kappa)) aQ_kappa = 0;
    double nu_break = 2.5e30/(B/100.);
    if (nu > nu_break)
    {
        return ((1. - eff) * a_Q_thermal(theta_e, n_e, nu, B, theta_B, jQ_thermal) + eff * aQ_kappa)/sqrt(nu/nu_break);
    }
    else{
        return ((1. - eff) * a_Q_thermal(theta_e, n_e, nu, B, theta_B, jQ_thermal) + eff * aQ_kappa);
    }

#elif (DF == POWER)
    return a_Q_power(theta_e, n_e, nu, B, theta_B,Power,Gamma_min,Gamma_max);
#elif (DF == TH)
    return a_Q_thermal(theta_e, n_e, nu, B, theta_B, jQ_thermal);
#elif (DF == VAR_POWER)
    double eff = get_efficiency(sigma, beta);
    if (eff<0.01) return a_Q_thermal(theta_e, n_e, nu, B, theta_B, jQ_thermal);
    else{
        double power,gamma_min,gamma_max,gamma_br,eta,ratio;
        eta=get_eta(theta_e,B,eff,beta,sigma,&power,&gamma_min,&gamma_max,&gamma_br,&ratio);
        if (Power_type == nth || gamma_br>gamma_max*0.9)
        return (a_Q_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B, jQ_thermal)+
                a_Q_power(theta_e, n_e*eta, nu, B, theta_B,power,gamma_min,gamma_max));
        else if (Power_type == cooling)
        return (a_Q_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B, jQ_thermal)+
                a_Q_power(theta_e, n_e*eta*ratio, nu, B, theta_B,power,gamma_min,gamma_br)+
                a_Q_power(theta_e, n_e*eta*(1-ratio), nu, B, theta_B,power+1.,gamma_br,gamma_max));
        }
#endif
}

double a_V_kappa(double theta_e, double n_e, double nu, double B,
                 double theta_B, double beta, double sigma) {
    double kappa = get_kappa(beta, sigma);

    double w = get_w_kappa(theta_e, beta, sigma, kappa);
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X_kappa = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    double A_low_factor_V =
        ((77. / (100 * w)) *
         pow((pow(sin(theta_B), -114. / 50.) - 1), 223. / 500.) *
         pow(X_kappa, -7. / 20.) * pow(kappa, -7. / 10.));
    double A_high_factor_V =
        ((143. / 10.) * pow(w, -116. / 125.) *
         pow((pow(sin(theta_B), -41. / 20.) - 1), 1. / 2.) *
         (pow(13, 2) * pow(kappa, -8) + (13. / 2500.) * kappa - (263. / 5000.) +
          (47. / (200 * kappa))) *
         pow(X_kappa, -1. / 2.));
    double A_S_factor_V = (sign(cos(theta_B)));
    double A_x_V = ((61. / 50.) * pow(kappa, -142. / 125.) + (7. / 1000.));
    double A_V_low_kappa = pow(X_kappa, -2. / 3.) * pow(3, 1. / 6.) *
                           (10. / 41.) * 2 * M_PI /
                           pow(w * kappa, 10. / 3. - kappa) * (kappa - 2.) *
                           (kappa - 1.) * kappa / (3. * kappa - 1.) *
                           tgamma(5. / 3.) * hyp2F1_f(theta_e, beta, sigma, kappa) * A_low_factor_V;
    double A_V_high_kappa =
        pow(X_kappa, -(1. + kappa) / 2.) * (pow(M_PI, 3. / 2.) / 3.) *
        ((kappa - 2.) * (kappa - 1.) * kappa / pow(w * kappa, 3.)) *
        (2 * tgamma(2. + kappa / 2.) / (2. + kappa) - 1.) * A_high_factor_V;
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);
    return prefac *
           pow((pow(A_V_low_kappa, -A_x_V) + pow(A_V_high_kappa, -A_x_V)),
               -1. / A_x_V) *
           A_S_factor_V;
}


double a_V_thermal(double theta_e, double n_e, double nu, double B,
                   double theta_B, double j_V_thermal) {
    double B_nu = planck_function(nu, theta_e); // Planck function
    return j_V_thermal / B_nu;
}

double a_V_power(double theta_e, double n_e, double nu, double B,double theta_B,
                 double power,double gamma_min,double gamma_max) {
    double nuc =
        ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double a_V_power_factor =
        pow(((71. * power / 100.) + 22. / 625.), 197. / 500.) *
        pow(((31. / 10.) * pow(sin(theta_B), -48. / 25.) - 31. / 10.),
            64. / 125.) *
        pow(nu / (nuc * sin(theta_B)), -1. / 2.) * sign(cos(theta_B));
    double prefac = n_e * ELECTRON_CHARGE * ELECTRON_CHARGE /
                    (ELECTRON_MASS * SPEED_OF_LIGHT * nu);
    return prefac * (pow(3., (power + 1) / 2.) * (power - 1.)) /
           (4. * (pow(gamma_min, 1. - power) - pow(gamma_max, 1. - power))) *
           tgamma((3. * power + 2.) / 12.) * tgamma((3. * power + 22.) / 12.) *
           pow((nu) / (nuc * sin(theta_B)), -(power + 2.) / 2.) *
           a_V_power_factor;
}

double a_V(double theta_e, double n_e, double nu, double B, double theta_B,
           double beta, double sigma) {
	double jV_thermal = j_V_thermal(theta_e, n_e, nu, B, theta_B);
#if (DF == KAPPA)
    return a_V_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
#elif (DF == VAR_KAPPA)
    double eff;
    eff = get_efficiency(sigma, beta);
    double aV_kappa = a_V_kappa(theta_e, n_e, nu, B, theta_B, beta, sigma);
    if (isnan(aV_kappa))
    {
        aV_kappa = 0;
    }
    double nu_break = 2.5e30/(B/100.);
    if (nu > nu_break)
    {
        return ((1. - eff) * a_V_thermal(theta_e, n_e, nu, B, theta_B, jV_thermal) + eff * aV_kappa)/sqrt(nu/nu_break);
    }
    else{
        return ((1. - eff) * a_V_thermal(theta_e, n_e, nu, B, theta_B, jV_thermal) + eff * aV_kappa);
    }

#elif (DF == POWER)
    return a_V_power(theta_e, n_e, nu, B, theta_B,Power,Gamma_min,Gamma_max);
#elif (DF == TH)
    return a_V_thermal(theta_e, n_e, nu, B, theta_B, jV_thermal);
#elif (DF == VAR_POWER)
    double eff = get_efficiency(sigma, beta);
    if (eff<0.01) return a_V_thermal(theta_e, n_e, nu, B, theta_B, jV_thermal);
    else{
        double power,gamma_min,gamma_max,gamma_br,eta,ratio;
        eta=get_eta(theta_e,B,eff,beta,sigma,&power,&gamma_min,&gamma_max,&gamma_br,&ratio);
        if (Power_type == nth || gamma_br>gamma_max*0.9)
        return (a_V_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B, jV_thermal)+
                a_V_power(theta_e, n_e*eta, nu, B, theta_B,power,gamma_min,gamma_max));
        else if (Power_type == cooling)
        return (a_V_thermal(theta_e, n_e*(1.-eta), nu, B, theta_B, jV_thermal)+
                a_V_power(theta_e, n_e*eta*ratio, nu, B, theta_B,power,gamma_min,gamma_br)+
                a_V_power(theta_e, n_e*eta*(1-ratio), nu, B, theta_B,power+1.,gamma_br,gamma_max));
        }
#endif
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
