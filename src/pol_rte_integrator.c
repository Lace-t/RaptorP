/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 *
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>

static void enforce_absorption_bound(double *aI, double *aQ, double *aU,
                                     double *aV);

// FUNCTIONS
////////////

// y contains the 4-position and the 4-velocity for one lightray/particle.
void f_parallel(double y[], double complex f_u[], double fvector[],
                double complex f_u_vector[]) {
    // Create variable (on the stack) for the connection
    double gamma_udd[4][4][4];
    // Einstein summation over indices v and w

    LOOP_ijk gamma_udd[i][j][k] = 0.;

    // Initialize position, four-velocity, and four-acceleration vectors based
    // on values of y
    double X_u[4] = {y[0], y[1], y[2], y[3]}; // X
    double U_u[4] = {y[4], y[5], y[6], y[7]}; // dX/dLambda
    double complex A_u[4] = {0., 0., 0., 0.}; // d^2X/dLambda^2

    // Obtain the Christoffel symbols at the current location
#if (metric == MKSBHAC || metric == MKSHARM || metric == LQMKS)
    connection_udd(X_u, gamma_udd);
#else
    connection_num_udd(X_u, gamma_udd);
#endif

    // Compute 4-acceleration using the geodesic equation
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * U_u[k];
    LOOP_i {
        fvector[i] = U_u[i];
        fvector[i + 4] = A_u[i];
    }

    // Reset A_u
    LOOP_i A_u[i] = 0.;

    // Compute f_u vector acceleration
    LOOP_ijk A_u[i] -= gamma_udd[i][j][k] * U_u[j] * f_u[k];
    LOOP_i { f_u_vector[i] = A_u[i]; }
}

void rk4_step_f(double y[], double complex f_u[], double dt) {
    // Array containing all "update elements" (4 times Nelements because RK4)
    double dx[4 * 2 * 4];
    double complex df[4 * 4];

    // Create a copy of the "y vector" that can be shifted for the
    // separate function calls made by RK4
    double yshift[4 * 2] = {y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7]};
    double complex f_u_shift[4] = {f_u[0], f_u[1], f_u[2], f_u[3]};

    // fvector contains f(yshift), as applied to yshift (the 'current' y
    // during RK steps). It is used to compute the 'k-coefficients' (dx)
    double fvector[4 * 2];
    double complex f_u_vector[4];

    // Compute the RK4 update coefficients ('K_n' in lit., 'dx' here)
    int i, q;
    double complex weights[4] = {0.5, 0.5, 1.,
                                 0.}; // Weights used for updating y
    for (q = 0; q < 4; q++) {
        f_parallel(
            yshift, f_u_shift, fvector,
            f_u_vector); // Apply function f to current y to obtain fvector
        for (i = 0; i < 4 * 2; i++) {
            dx[q * 4 * 2 + i] = dt * fvector[i]; // Use fvector to fill dx
            yshift[i] = y[i] + dx[q * 4 * 2 + i] * weights[q]; // Update y
        }
        for (i = 0; i < 4; i++) {
            df[q * 4 + i] = dt * f_u_vector[i];
            f_u_shift[i] = f_u[i] + df[q * 4 + i] * weights[q];
        }
    }

    // Update the y-vector (light ray)
    for (i = 0; i < 4 * 2; i++) {
        y[i] = y[i] + 1. / 6. *
                          (dx[0 * 4 * 2 + i] + dx[1 * 4 * 2 + i] * 2. +
                           dx[2 * 4 * 2 + i] * 2. + dx[3 * 4 * 2 + i]);
    }

    // Update the f-vector (polarization)
    for (i = 0; i < 4; i++) {
        f_u[i] = f_u[i] + 1. / 6. *
                              (df[0 * 4 + i] + df[1 * 4 + i] * 2. +
                               df[2 * 4 + i] * 2. + df[3 * 4 + i]);
    }
}

void f_tetrad_to_stokes(double Iinv, double Iinv_pol,
                        double complex f_tetrad_u[], double complex S_A[]) {
    S_A[0] = Iinv;
    S_A[1] = Iinv_pol * (cabs(f_tetrad_u[1]) * cabs(f_tetrad_u[1]) -
                         cabs(f_tetrad_u[2]) * cabs(f_tetrad_u[2]));
    S_A[2] = Iinv_pol * (conj(f_tetrad_u[1]) * f_tetrad_u[2] +
                         f_tetrad_u[1] * conj(f_tetrad_u[2]));
    S_A[3] = Iinv_pol * (I * (conj(f_tetrad_u[1]) * f_tetrad_u[2] -
                              f_tetrad_u[1] * conj(f_tetrad_u[2])));
}

void stokes_to_f_tetrad(double complex S_A[], double *Iinv, double *Iinv_pol,
                        double complex f_tetrad_u[]) {

    *Iinv = S_A[0];

    *Iinv_pol = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]);

    double Qnorm = S_A[1] / (*Iinv_pol);
    double Unorm = S_A[2] / (*Iinv_pol);
    double Vnorm = S_A[3] / (*Iinv_pol);

    // source:
    // https://physics.stackexchange.com/questions/238957/converting-stokes-parameters-to-jones-vector
    f_tetrad_u[1] = sqrt((1. + Qnorm) / 2.);

    if (f_tetrad_u[1] == 0)
        f_tetrad_u[2] = 1.;
    else
        f_tetrad_u[2] =
            Unorm / (2. * f_tetrad_u[1]) - I * Vnorm / (2. * f_tetrad_u[1]);
}

// NOTE: works only in Kerr metric
// Ziri's suggestion: construct U vecs
void construct_U_vector(double X_u[], double U_u[]) {
// Obtain relevant metric terms:
#if (metric == CKS)
    double U_KS[4];
    double X_KS[4];

    CKS_to_KS(X_u, X_KS);

    double g_uu[4][4];
    metric_KS_uu(X_KS, g_uu);
    double g_uu00 = g_uu[0][0];
    double g_uu03 = g_uu[0][3];
    double g_uu33 = g_uu[3][3];

#else
    double g_uu[4][4];
    metric_uu(X_u, g_uu);
    double g_uu00 = g_uu[0][0];
    double g_uu03 = g_uu[0][3];
    double g_uu33 = g_uu[3][3];
#endif
    // Observer/plasma wave vector:
    double U_d[4] = {-1., 0., 0., 0.};
    double B__ = -g_uu03 * U_d[0] / g_uu33;
    double C__ = -(1. + g_uu00 * U_d[0] * U_d[0]) / g_uu33;

    // Properly normalize U_u:
    U_d[3] = B__ + sqrt(B__ * B__ + C__);

#if (metric == CKS)
    LOOP_i {
        U_KS[i] = 0.;
        U_u[i] = 0;
    }
    raise_index_KS(X_KS, U_d, U_KS);

    double coordKS[8];
    double coordCKS[8];

    LOOP_i {
        coordKS[i] = X_KS[i];
        coordKS[i + 4] = U_KS[i];
    }

    KS_to_CKS_u(coordKS, coordCKS);

    LOOP_i U_u[i] = coordCKS[i + 4];

#else
    LOOP_i U_u[i] = 0.;
    raise_index(X_u, U_d, U_u);

#endif
}

// NEW FUNCTIONS JUNE 2021
//////////////////////////

// Transform f_tetrad_u to f_u
void f_tetrad_to_f(double complex *f_u, double tetrad_u[][4],
                   double complex *f_tetrad_u) {

    LOOP_i f_u[i] = 0.;
    LOOP_ij f_u[i] += tetrad_u[i][j] * f_tetrad_u[j];
}

// Transform f_u to f_tetrad_u
void f_to_f_tetrad(double complex *f_tetrad_u, double tetrad_d[][4],
                   double complex *f_u) {

    LOOP_i f_tetrad_u[i] = 0.;
    LOOP_ij f_tetrad_u[i] += tetrad_d[j][i] * f_u[j];
}


void evaluate_coeffs_user(double *jI, double *jQ, double *jU, double *jV,
                            double *rQ, double *rU, double *rV, double *aI,
                            double *aQ, double *aU, double *aV, double nu_p,
                            struct GRMHD modvar, double pitch_ang) {

    *jI = (j_I(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *jQ = (j_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *jU = 0.;
    *jV = (j_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));

    *rQ = (rho_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *rU = 0.;
    *rV = rho_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma);

    *aI = (a_I(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *aQ = (a_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *aU = 0;
    *aV = (a_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));

    // *rQ=0.;
    // *rU=0.;
    // *rV=0.;
    // *aI=0.;
    // *aQ=0.;
    // *aU=0.;
    // *aV=0.;

    // Transform to invariant forms
    *jI /= (nu_p * nu_p);
    *jQ /= (nu_p * nu_p);
    *jV /= (nu_p * nu_p);

    *aI *= nu_p;
    *aQ *= nu_p;
    *aV *= nu_p;

    *rQ *= nu_p;
    *rV *= nu_p;

    // somtimes in very specific cells issue with Ipol>S_I, numerical round off
    // issues
    //  and/or scheme handling pure polarizaiton states poorly. renormalizing to
    //  ensure pol_frac<1.0.

    double pol_frac = sqrt((*jQ) * (*jQ) + (*jV) * (*jV)) / (*jI);
    if (pol_frac > 1.) {
        *jQ /= (pol_frac + 0.005);
        *jU /= (pol_frac + 0.005);
        *jV /= (pol_frac + 0.005);
    }

    if (!isfinite(*jI)) *jI = 0.;
    if (!isfinite(*jQ)) *jQ = 0.;
    if (!isfinite(*jU)) *jU = 0.;
    if (!isfinite(*jV)) *jV = 0.;
    if (!isfinite(*rQ)) *rQ = 0.;
    if (!isfinite(*rU)) *rU = 0.;
    if (!isfinite(*rV)) *rV = 0.;
    if (!isfinite(*aI)) *aI = 0.;
    if (!isfinite(*aQ)) *aQ = 0.;
    if (!isfinite(*aU)) *aU = 0.;
    if (!isfinite(*aV)) *aV = 0.;
    enforce_absorption_bound(aI, aQ, aU, aV);
    /*
        pol_frac = sqrt((*aQ) * (*aQ) + (*aV) * (*aV)) / (*aI);
        if (pol_frac > 1.) {
            *aQ *= (pol_frac - 0.005);
            *aU *= (pol_frac - 0.005);
            *aV *= (pol_frac - 0.005);
        }
    */
}
void evaluate_coeffs_single(double *jI, double *jQ, double *jU, double *jV,
                            double *rQ, double *rU, double *rV, double *aI,
                            double *aQ, double *aU, double *aV, double nu_p,
                            struct GRMHD modvar, double pitch_ang) {
    *jI = (j_I(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *jQ = (j_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *jU = 0.;
    *jV = (j_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));

    *rQ = (rho_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *rU = 0.;
    *rV = rho_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma);

    *aI = (a_I(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *aQ = (a_Q(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    *aU = 0;
    *aV = (a_V(modvar.theta_e, modvar.n_e, nu_p, modvar.B, pitch_ang, modvar.beta, modvar.sigma));
    // *rQ=0.;
    // *rU=0.;
    // *rV=0.;
    // *aI=0.;
    // *aQ=0.;
    // *aU=0.;
    // *aV=0.;

    // Transform to invariant forms
    *jI /= (nu_p * nu_p);
    *jQ /= (nu_p * nu_p);
    *jV /= (nu_p * nu_p);

    *aI *= nu_p;
    *aQ *= nu_p;
    *aV *= nu_p;

    *rQ *= nu_p;
    *rV *= nu_p;

    // somtimes in very specific cells issue with Ipol>S_I, numerical round off
    // issues
    //  and/or scheme handling pure polarizaiton states poorly. renormalizing to
    //  ensure pol_frac<1.0.

    double pol_frac = sqrt((*jQ) * (*jQ) + (*jV) * (*jV)) / (*jI);
    if (pol_frac > 1.) {
        *jQ /= (pol_frac + 0.005);
        *jU /= (pol_frac + 0.005);
        *jV /= (pol_frac + 0.005);
    }

    if (!isfinite(*jI)) *jI = 0.;
    if (!isfinite(*jQ)) *jQ = 0.;
    if (!isfinite(*jU)) *jU = 0.;
    if (!isfinite(*jV)) *jV = 0.;
    if (!isfinite(*rQ)) *rQ = 0.;
    if (!isfinite(*rU)) *rU = 0.;
    if (!isfinite(*rV)) *rV = 0.;
    if (!isfinite(*aI)) *aI = 0.;
    if (!isfinite(*aQ)) *aQ = 0.;
    if (!isfinite(*aU)) *aU = 0.;
    if (!isfinite(*aV)) *aV = 0.;
    enforce_absorption_bound(aI, aQ, aU, aV);
    /*
        pol_frac = sqrt((*aQ) * (*aQ) + (*aV) * (*aV)) / (*aI);
        if (pol_frac > 1.) {
            *aQ *= (pol_frac - 0.005);
            *aU *= (pol_frac - 0.005);
            *aV *= (pol_frac - 0.005);
        }
    */
}

double check_stiffness(double jI, double jQ, double jU, double jV,
                       double rQ, double rU, double rV,
                       double aI, double aQ, double aU, double aV,
                       double C) {
    double a2 = rQ*rQ + rV*rV - aQ*aQ - aV*aV;
    double a0 = -2.*aV*aQ*rV*rQ - aQ*aQ*rQ*rQ - aV*aV*rV*rV;

    complex double zplus  = (-a2 + csqrt(a2*a2 - 4.*a0)) / 2.;
    complex double zminus = (-a2 - csqrt(a2*a2 - 4.*a0)) / 2.;

    complex double l[4];
    l[0] = aI + csqrt(zplus);
    l[1] = aI - csqrt(zplus);
    l[2] = aI + csqrt(zminus);
    l[3] = aI - csqrt(zminus);

    double max_rate = 0.0;

    for (int i = 0; i < 4; ++i) {
        double rate = cabs(l[i]);   // |λ_i|
        if (rate > max_rate)
            max_rate = rate;
    }

    if (max_rate <= 0.0)
        return 1e40;   // effectively unlimited

    return EPS_RTE / (C * max_rate);
}

void pol_rte_rk4_step(double jI, double jQ, double jU, double jV, double rQ,
                      double rU, double rV, double aI, double aQ, double aU,
                      double aV, double dl_current, double C,
                      double complex S_A[]) {
    double complex I0 = S_A[0];
    double complex Q0 = S_A[1];
    double complex U0 = S_A[2];
    double complex V0 = S_A[3];

    // RK4 with constant coefficients
    // k1
    double complex Ik1 =
        dl_current * C * jI -
        dl_current * C * (aI * I0 + aQ * Q0 + aU * U0 + aV * V0);
    double complex Qk1 =
        dl_current * C * jQ -
        dl_current * C * (aQ * I0 + aI * Q0 + rV * U0 - rU * V0);
    double complex Uk1 =
        dl_current * C * jU -
        dl_current * C * (aU * I0 - rV * Q0 + aI * U0 + rQ * V0);
    double complex Vk1 =
        dl_current * C * jV -
        dl_current * C * (aV * I0 + rU * Q0 - rQ * U0 + aI * V0);

    // k2
    double complex Ik2 = dl_current * C * jI -
                         dl_current * C *
                             (aI * (I0 + 0.5 * Ik1) + aQ * (Q0 + 0.5 * Qk1) +
                              aU * (U0 + 0.5 * Uk1) + aV * (V0 + 0.5 * Vk1));
    double complex Qk2 = dl_current * C * jQ -
                         dl_current * C *
                             (aQ * (I0 + 0.5 * Ik1) + aI * (Q0 + 0.5 * Qk1) +
                              rV * (U0 + 0.5 * Uk1) - rU * (V0 + 0.5 * Vk1));
    double complex Uk2 = dl_current * C * jU -
                         dl_current * C *
                             (aU * (I0 + 0.5 * Ik1) - rV * (Q0 + 0.5 * Qk1) +
                              aI * (U0 + 0.5 * Uk1) + rQ * (V0 + 0.5 * Vk1));
    double complex Vk2 = dl_current * C * jV -
                         dl_current * C *
                             (aV * (I0 + 0.5 * Ik1) + rU * (Q0 + 0.5 * Qk1) -
                              rQ * (U0 + 0.5 * Uk1) + aI * (V0 + 0.5 * Vk1));

    // k3
    double complex Ik3 = dl_current * C * jI -
                         dl_current * C *
                             (aI * (I0 + 0.5 * Ik2) + aQ * (Q0 + 0.5 * Qk2) +
                              aU * (U0 + 0.5 * Uk2) + aV * (V0 + 0.5 * Vk2));
    double complex Qk3 = dl_current * C * jQ -
                         dl_current * C *
                             (aQ * (I0 + 0.5 * Ik2) + aI * (Q0 + 0.5 * Qk2) +
                              rV * (U0 + 0.5 * Uk2) - rU * (V0 + 0.5 * Vk2));
    double complex Uk3 = dl_current * C * jU -
                         dl_current * C *
                             (aU * (I0 + 0.5 * Ik2) - rV * (Q0 + 0.5 * Qk2) +
                              aI * (U0 + 0.5 * Uk2) + rQ * (V0 + 0.5 * Vk2));
    double complex Vk3 = dl_current * C * jV -
                         dl_current * C *
                             (aV * (I0 + 0.5 * Ik2) + rU * (Q0 + 0.5 * Qk2) -
                              rQ * (U0 + 0.5 * Uk2) + aI * (V0 + 0.5 * Vk2));

    // k4
    double complex Ik4 =
        dl_current * C * jI - dl_current * C *
                                  (aI * (I0 + Ik3) + aQ * (Q0 + Qk3) +
                                   aU * (U0 + Uk3) + aV * (V0 + Vk3));
    double complex Qk4 =
        dl_current * C * jQ - dl_current * C *
                                  (aQ * (I0 + Ik3) + aI * (Q0 + Qk3) +
                                   rV * (U0 + Uk3) - rU * (V0 + Vk3));
    double complex Uk4 =
        dl_current * C * jU - dl_current * C *
                                  (aU * (I0 + Ik3) - rV * (Q0 + Qk3) +
                                   aI * (U0 + Uk3) + rQ * (V0 + Vk3));
    double complex Vk4 =
        dl_current * C * jV - dl_current * C *
                                  (aV * (I0 + Ik3) + rU * (Q0 + Qk3) -
                                   rQ * (U0 + Uk3) + aI * (V0 + Vk3));

    S_A[0] = I0 + 1. / 6. * (Ik1 + 2. * Ik2 + 2. * Ik3 + Ik4);
    S_A[1] = Q0 + 1. / 6. * (Qk1 + 2. * Qk2 + 2. * Qk3 + Qk4);
    S_A[2] = U0 + 1. / 6. * (Uk1 + 2. * Uk2 + 2. * Uk3 + Uk4);
    S_A[3] = V0 + 1. / 6. * (Vk1 + 2. * Vk2 + 2. * Vk3 + Vk4);
}

void pol_rte_trapezoid_step(double jI, double jQ, double jU, double jV,
                            double rQ, double rU, double rV, double aI,
                            double aQ, double aU, double aV, double dl_current,
                            double C, double complex S_A[]) {
    double complex I0 = S_A[0];
    double complex Q0 = S_A[1];
    double complex U0 = S_A[2];
    double complex V0 = S_A[3];

    double u11 = 1. + 0.5 * dl_current * C * aI;
    
    // Check for numerical issues in diagonal elements
    // if (fabs(u11) < 1e-10) {
    //     FILE *fp = fopen("diagnosis.txt", "a");
    //     fprintf(fp, "WARNING: u11 too small: %g, dl_current=%g, C=%g, aI=%g\n", u11, dl_current, C, aI);
    //     fclose(fp);
    // }
    
    double u12 = 0.5 * dl_current * C * aQ;
    double u14 = 0.5 * dl_current * C * aV;
    double l21 = 0.5 * dl_current * C * aQ / u11;
    double u22 = 1. + 0.5 * dl_current * C * aI - l21 * u12;
    
    // if (fabs(u22) < 1e-10) {
    //     FILE *fp = fopen("diagnosis.txt", "a");
    //     fprintf(fp, "WARNING: u22 too small: %g\n", u22);
    //     fclose(fp);
    // }
    
    double u23 = 0.5 * dl_current * C * rV;
    double u24 = -l21 * u14;
    double l32 = -0.5 * dl_current * C * rV / u22;
    double u33 = 1. + 0.5 * dl_current * C * aI - l32 * u23;
    
    // if (fabs(u33) < 1e-10) {
    //     FILE *fp = fopen("diagnosis.txt", "a");
    //     fprintf(fp, "WARNING: u33 too small: %g\n", u33);
    //     fclose(fp);
    // }
    
    double u34 = 0.5 * dl_current * C * rQ - l32 * u24;
    double l41 = 0.5 * dl_current * C * aV / u11;
    double l42 = -l41 * u12 / u22;
    double l43 = (-0.5 * dl_current * C * rQ - l42 * u23) / u33;
    double u44 =
        1. + 0.5 * dl_current * C * aI - l41 * u14 - l42 * u24 - l43 * u34;
    
    // if (fabs(u44) < 1e-10) {
    //     FILE *fp = fopen("diagnosis.txt", "a");
    //     fprintf(fp, "WARNING: u44 too small: %g\n", u44);
    //     fclose(fp);
    // }

    // Construct b-vector.
    double b1 =
        I0 + dl_current * C / 2. * (2. * jI - (aI * I0 + aQ * Q0 + aV * V0));
    double b2 =
        Q0 + dl_current * C / 2. * (2. * jQ - (aQ * I0 + aI * Q0 + rV * U0));
    double b3 =
        U0 + dl_current * C / 2. * (2. * jU - (-rV * Q0 + aI * U0 + rQ * V0));
    double b4 =
        V0 + dl_current * C / 2. * (2. * jV - (aV * I0 - rQ * U0 + aI * V0));

    // Construct y.
    double y1 = b1;
    double y2 = b2 - l21 * y1;
    double y3 = b3 - l32 * y2;
    double y4 = b4 - l41 * y1 - l42 * y2 - l43 * y3;

    // Construct x.
    double x4 = y4 / u44;
    double x3 = (y3 - u34 * x4) / u33;
    double x2 = (y2 - u23 * x3 - u24 * x4) / u22;
    double x1 = (y1 - u12 * x2 - u14 * x4) / u11;

    S_A[0] = x1;
    S_A[1] = x2;
    S_A[2] = x3;
    S_A[3] = x4;
    
    // Check if result is NaN
    // if (isnan(creal(S_A[0]))) {
    //     FILE *fp = fopen("diagnosis.txt", "a");
    //     fprintf(fp, "ERROR in trapezoid_step: S_A[0] is NaN\n");
    //     fprintf(fp, "  Input: jI=%g jQ=%g jV=%g aI=%g aQ=%g aV=%g\n", jI, jQ, jV, aI, aQ, aV);
    //     fprintf(fp, "  Input: rQ=%g rV=%g dl_current=%g C=%g\n", rQ, rV, dl_current, C);
    //     fprintf(fp, "  Initial: I0=%g Q0=%g U0=%g V0=%g\n", creal(I0), creal(Q0), creal(U0), creal(V0));
    //     fprintf(fp, "  u11=%g u22=%g u33=%g u44=%g\n", u11, u22, u33, u44);
    //     fprintf(fp, "  x1=%g x2=%g x3=%g x4=%g\n", x1, x2, x3, x4);
    //     fclose(fp);
    // }
}

void pol_rte_euler_step(double jI, double jQ, double jU, double jV,
                        double rQ, double rU, double rV,
                        double aI, double aQ, double aU, double aV,
                        double dl_current, double C,
                        double complex S_A[]){
    double complex I0 = S_A[0];
    double complex Q0 = S_A[1];
    double complex U0 = S_A[2];
    double complex V0 = S_A[3];

    double complex dI =
        C * ( jI - (aI * I0 + aQ * Q0 + aU * U0 + aV * V0) );

    double complex dQ =
        C * ( jQ - (aQ * I0 + aI * Q0 + rV * U0 - rU * V0) );

    double complex dU =
        C * ( jU - (aU * I0 - rV * Q0 + aI * U0 + rQ * V0) );

    double complex dV =
        C * ( jV - (aV * I0 + rU * Q0 - rQ * U0 + aI * V0) );

    S_A[0] = I0 + dl_current * dI;
    S_A[1] = Q0 + dl_current * dQ;
    S_A[2] = U0 + dl_current * dU;
    S_A[3] = V0 + dl_current * dV;
}

void f_to_stokes(double complex f_u[], double complex f_tetrad_u[],
                 double tetrad_d[][4], double complex S_A[], double Iinv,
                 double Iinv_pol) {
    f_to_f_tetrad(f_tetrad_u, tetrad_d, f_u);
    f_tetrad_to_stokes(Iinv, Iinv_pol, f_tetrad_u, S_A);
}

void stokes_to_f(double complex f_u[], double complex f_tetrad_u[],
                 double tetrad_u[][4], double complex S_A[], double *Iinv,
                 double *Iinv_pol) {
    stokes_to_f_tetrad(S_A, Iinv, Iinv_pol, f_tetrad_u);

    f_tetrad_to_f(f_u, tetrad_u, f_tetrad_u);
}

static inline void advance_polarized_step_position(double X_u[],
                                                   const double k_u_old[],
                                                   double dl_accumulated,
                                                   double affine_scale) {
    if (dl_accumulated <= 0.)
        return;

    LOOP_i X_u[i] += k_u_old[i] * (dl_accumulated / affine_scale);
}

static const double MAX_ADAPTIVE_SUBSTEPS = 1.e4;

static inline double capped_adaptive_dl(double dl_remain, double dl_opt) {
    if (!(dl_opt > 0.))
        return dl_remain;

    if (dl_remain / dl_opt > MAX_ADAPTIVE_SUBSTEPS)
        return dl_remain / MAX_ADAPTIVE_SUBSTEPS;

    return dl_opt;
}

static inline double distance_to_next_coeff_refresh(double dl_remain,
                                                    double dl_total,
                                                    double dl_mark) {
    if (!(dl_mark > 0.))
        return dl_remain;

    double next_refresh_remain = dl_mark * dl_total;

    if (dl_remain > next_refresh_remain)
        return dl_remain - next_refresh_remain;

    return dl_remain;
}

static int stokes_has_nonfinite(const double complex S_A[]) {
    for (int i = 0; i < 4; ++i) {
        if (!isfinite(creal(S_A[i])) || !isfinite(cimag(S_A[i])))
            return 1;
    }

    return 0;
}

static void dump_bad_stokes_state(const char *branch, const char *solver,
                                  int polarization_active,
                                  double dl_current, double C,
                                  double complex S_before[],
                                  double complex S_after[],
                                  double Iinv_before, double Iinv_pol_before,
                                  double jI, double jQ, double jV,
                                  double rQ, double rV,
                                  double aI, double aQ, double aV,
                                  double X_u[], struct GRMHD modvar_local,
                                  double nu_p, double pitch_ang) {
    printf("bad_stokes branch=%s solver=%s pol=%d dl=%.2g dlC=%.2g\n",
           branch, solver, polarization_active, dl_current, dl_current * C);
    printf("r=%.2g X=%.2g %.2g %.2g n_e=%.2g theta_e=%.2g b=%.2g nu_p=%.2g pitch_ang=%.2g\n",
           sqrt(X_u[1] * X_u[1] + X_u[2] * X_u[2]), X_u[1], X_u[2], X_u[3],
           modvar_local.n_e, modvar_local.theta_e, modvar_local.B, nu_p,
           pitch_ang);
    printf("Iinv_before=%.2g Iinv_pol_before=%.2g\n",
           Iinv_before, Iinv_pol_before);
    printf("S_before=%.2g %.2g %.2g %.2g\n",
           creal(S_before[0]), creal(S_before[1]), creal(S_before[2]),
           creal(S_before[3]));
    printf("jI=%.2e jQ=%.2e jV=%.2e rQ=%.2e rV=%.2e aI=%.2e aQ=%.2e aV=%.2e\n",
           jI, jQ, jV, rQ, rV, aI, aQ, aV);
    printf("S_after=%.2g %.2g %.2g %.2g\n",
           creal(S_after[0]), creal(S_after[1]), creal(S_after[2]),
           creal(S_after[3]));
}

static void enforce_absorption_bound(double *aI, double *aQ, double *aU,
                                     double *aV) {
    if (!isfinite(*aI) || *aI <= 0.) {
        *aI = 0.;
        *aQ = 0.;
        *aU = 0.;
        *aV = 0.;
        return;
    }

    double a_pol = sqrt((*aQ) * (*aQ) + (*aU) * (*aU) + (*aV) * (*aV));
    if (!isfinite(a_pol) || a_pol <= 0.)
        return;

    if (a_pol > *aI) {
        double scale = 0.999 * (*aI) / a_pol;
        *aQ *= scale;
        *aU *= scale;
        *aV *= scale;
    }
}

static int pol_rte_analytic_step(double jI, double jQ, double jU, double jV,
                                 double rQ, double rU, double rV, double aI,
                                 double aQ, double aU, double aV,
                                 double dl_current, double C,
                                 double complex S_A[]) {
    double H_data[25] = {0.};
    double expH_data[25] = {0.};
    gsl_matrix_view H = gsl_matrix_view_array(H_data, 5, 5);
    gsl_matrix_view expH = gsl_matrix_view_array(expH_data, 5, 5);
    gsl_error_handler_t *old_handler = gsl_set_error_handler_off();

    const double step_scale = dl_current * C;

    // Propagate one constant-coefficient segment exactly via exp(dl * H).
    gsl_matrix_set(&H.matrix, 0, 0, -step_scale * aI);
    gsl_matrix_set(&H.matrix, 0, 1, -step_scale * aQ);
    gsl_matrix_set(&H.matrix, 0, 2, -step_scale * aU);
    gsl_matrix_set(&H.matrix, 0, 3, -step_scale * aV);
    gsl_matrix_set(&H.matrix, 0, 4,  step_scale * jI);

    gsl_matrix_set(&H.matrix, 1, 0, -step_scale * aQ);
    gsl_matrix_set(&H.matrix, 1, 1, -step_scale * aI);
    gsl_matrix_set(&H.matrix, 1, 2, -step_scale * rV);
    gsl_matrix_set(&H.matrix, 1, 3,  step_scale * rU);
    gsl_matrix_set(&H.matrix, 1, 4,  step_scale * jQ);

    gsl_matrix_set(&H.matrix, 2, 0, -step_scale * aU);
    gsl_matrix_set(&H.matrix, 2, 1,  step_scale * rV);
    gsl_matrix_set(&H.matrix, 2, 2, -step_scale * aI);
    gsl_matrix_set(&H.matrix, 2, 3, -step_scale * rQ);
    gsl_matrix_set(&H.matrix, 2, 4,  step_scale * jU);

    gsl_matrix_set(&H.matrix, 3, 0, -step_scale * aV);
    gsl_matrix_set(&H.matrix, 3, 1, -step_scale * rU);
    gsl_matrix_set(&H.matrix, 3, 2,  step_scale * rQ);
    gsl_matrix_set(&H.matrix, 3, 3, -step_scale * aI);
    gsl_matrix_set(&H.matrix, 3, 4,  step_scale * jV);

    int status = gsl_linalg_exponential_ss(&H.matrix, &expH.matrix,
                                           GSL_PREC_DOUBLE);
    gsl_set_error_handler(old_handler);

    if (status != GSL_SUCCESS)
        return 0;

    double y_re0[5] = {creal(S_A[0]), creal(S_A[1]), creal(S_A[2]),
                       creal(S_A[3]), 1.0};
    double y_im0[5] = {cimag(S_A[0]), cimag(S_A[1]), cimag(S_A[2]),
                       cimag(S_A[3]), 0.0};

    for (int i = 0; i < 4; ++i) {
        double re1 = 0.;
        double im1 = 0.;
        for (int j = 0; j < 5; ++j) {
            double Tij = gsl_matrix_get(&expH.matrix, i, j);
            re1 += Tij * y_re0[j];
            im1 += Tij * y_im0[j];
        }
        S_A[i] = re1 + I * im1;
    }

    return 1;
}

void pol_integration_step(struct GRMHD modvar, double frequency,
                          double *dl_current0, double C, double X_u[],
                          double k_u[], double k_d[], int *POLARIZATION_ACTIVE,
                          double complex f_u[], double complex f_tetrad_u[],
                          double tetrad_d[][4], double tetrad_u[][4],
                          double complex S_A[], double *Iinv, double *Iinv_pol,
                          double *tau, double *tauF) {

    double jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV;
    double pitch_ang, nu_p;
    double k_u_old[4],X_u_old[4];
    struct GRMHD modvar_local = modvar;

	double affine_scale=(ELECTRON_MASS*SPEED_OF_LIGHT*SPEED_OF_LIGHT)/(PLANCK_CONSTANT*frequency);
	double dl_opt=1e40,dl_current;
    // Unpolarized: 1) Create light path by integration. 2) For each
    // step in lightpath, perform one radiative transfer step.
    // Polarized:   1) Create light path by integration. 2) For each
    // step in lightpath, perform one radiative transfer step, AND,
    // OUTSIDE in_volume loop, do a spacetime propagation step.

    // TRANSFER STEP
    ////////////////

    // Obtain pitch angle: still no units (geometric)
    pitch_ang = pitch_angle(X_u, k_u, modvar_local.B_u, modvar_local.U_u);

    // perfect field alignment, no emission
    if (fmod(pitch_ang, M_PI) == 0)
        return;

    // CGS UNITS USED FROM HERE ON OUT
    //////////////////////////////////
    LOOP_i {
        k_u_old[i] = k_u[i];
        X_u_old[i] = X_u[i];
    }

    // Scale the wave vector to correct energy
    LOOP_i k_u[i] *= PLANCK_CONSTANT * frequency /
                     (ELECTRON_MASS * SPEED_OF_LIGHT * SPEED_OF_LIGHT);

    // Convert distance dlambda accordingly
    *dl_current0 *= affine_scale ;

    // lower the index of the wavevector
    lower_index(X_u, k_u, k_d);

    // Compute the photon frequency in the plasma frame:
    nu_p = freq_in_plasma_frame(modvar_local.U_u, k_d);


	//adaptive optical depth stepsize
	#if (EMISUSER)
    evaluate_coeffs_user(&jI, &jQ, &jU, &jV, &rQ, &rU, &rV, &aI,
                         &aQ, &aU, &aV, nu_p, modvar_local, pitch_ang);
    #else
    evaluate_coeffs_single(&jI, &jQ, &jU, &jV, &rQ, &rU, &rV, &aI,
                           &aQ, &aU, &aV, nu_p, modvar_local, pitch_ang);
    #endif
    dl_opt=check_stiffness(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV, C);
    double step_ratio = (dl_opt > 0.) ? (*dl_current0 / dl_opt) : INFINITY;

	if (dl_opt>*dl_current0){
	//if (1) {
	    dl_current=*dl_current0;
        double complex S_before_step[4] = {S_A[0], S_A[1], S_A[2], S_A[3]};
        const char *solver_used = "trapezoid";
        double Iinv_before = *Iinv;
        double Iinv_pol_before = *Iinv_pol;

        // Create tetrad, needed whether POLARIZATION_ACTIVE is true or
        // false.
        create_observer_tetrad(X_u, k_u, modvar_local.U_u, modvar_local.B_u, tetrad_u);
        LOOP_ij if (isnan(tetrad_u[i][j])) return;
        
        create_tetrad_d(X_u, tetrad_u, tetrad_d);

        // FILE *fp=fopen("diagnosis.txt", "a");
        // fprintf(fp,"tetrad_u=%g+%gi %g+%gi %g+%gi tetrad_d=%g+%gi %g+%gi %g+%gi %g+%gi\n",
        //             creal(tetrad_u[0][0]),cimag(tetrad_u[0][0]),creal(tetrad_u[0][1]),cimag(tetrad_u[0][1]),
        //             creal(tetrad_u[0][2]),cimag(tetrad_u[0][2]),creal(tetrad_u[0][3]),cimag(tetrad_u[0][3]),
        //             creal(tetrad_d[0][0]),cimag(tetrad_d[0][0]),creal(tetrad_d[0][1]),cimag(tetrad_d[0][1]),
        //             creal(tetrad_d[0][2]),cimag(tetrad_d[0][2]),creal(tetrad_d[0][3]),cimag(tetrad_d[0][3]));
        // fclose(fp);

        // FROM F VECTOR TO STOKES (when applicable)
        ////////////////////////////////////////////

        // If (POLARIZATION_ACTIVE), get Stokes params from f_u and p.
        // (Otherwise, never been in volume before; we simply use
        // S_I_current)


        if (*POLARIZATION_ACTIVE) {
            f_to_stokes(f_u, f_tetrad_u, tetrad_d, S_A, *Iinv, *Iinv_pol);
        }
        
        // Given Stokes params and plasma coeffs, compute NEW Stokes params
        // after plasma step.

        // If both rotation coeffs (times dlambda) are smaller than
        // threshold, take an RK4 step; otherwise, implicit Euler.
        // if (fabs(rQ) < THRESH && fabs(rV) < THRESH) {
        
        // Use RK4 for more robust numerical integration, especially when
        // coefficients are very small or when trapezoid might be unstable
        double trap_threshold = 0.5 * dl_opt;
        double param_magnitude = fabs(aI) + fabs(aQ) + fabs(aV) + fabs(rQ) + fabs(rV);
        
        // Use RK4 for non-stiff segments. Once dl*C*coeff becomes O(1),
        // prefer the implicit trapezoid step to avoid explicit blow-ups.
        if (*dl_current0 < trap_threshold ||
            (fabs(aI * dl_current * C) <= 0.5 &&
             param_magnitude * dl_current * C <= 1.0) ||
            fabs(aI) < 1e-20) {
            solver_used = "rk4";
            pol_rte_rk4_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV,
                             dl_current, C, S_A);
        } else {
            pol_rte_trapezoid_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV,
                                   dl_current, C, S_A);
        }

        if (stokes_has_nonfinite(S_A)) {
            dump_bad_stokes_state("single", solver_used, *POLARIZATION_ACTIVE,
                                  dl_current, C, S_before_step, S_A,
                                  Iinv_before, Iinv_pol_before, jI, jQ, jV,
                                  rQ, rV, aI, aQ, aV, X_u, modvar_local, nu_p,
                                  pitch_ang);
        }

        // FROM STOKES TO F VECTOR
        ///////////////////////////
        // somtimes in very specific cells issue with Ipol>S_I, numerical round off
        // issues? renormalizing.
        
        // Check for NaN before proceeding
        // if (isnan(creal(S_A[0])) || isnan(creal(S_A[1])) || 
        //     isnan(creal(S_A[2])) || isnan(creal(S_A[3]))) {
        //     FILE *fp = fopen("diagnosis.txt", "a");
        //     fprintf(fp, "ERROR: NaN detected in S_A after integration step\n");
        //     fprintf(fp, "  RK4 was used: %d\n", (fabs(aI * dl_current * C) > 0.5 || param_magnitude * dl_current * C > 1.0));
        //     fprintf(fp, "  aI=%g aQ=%g aV=%g\n", aI, aQ, aV);
        //     fprintf(fp, "  rQ=%g rV=%g\n", rQ, rV);
        //     fprintf(fp, "  S_A[0]=%g+%gi S_A[1]=%g+%gi\n", creal(S_A[0]), cimag(S_A[0]), creal(S_A[1]), cimag(S_A[1]));
        //     fclose(fp);
        //     // Reset to previous state
        //     S_A[0] = S_A_initial[0];
        //     S_A[1] = S_A_initial[1];
        //     S_A[2] = S_A_initial[2];
        //     S_A[3] = S_A_initial[3];
        // }
        
        double pol_frac =
            sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]) /
            sqrt(S_A[0] * S_A[0]);

        if (isnan(pol_frac) || isinf(pol_frac)) {
            // Handle NaN/Inf in polarization fraction
            S_A[1] = 0.;
            S_A[2] = 0.;
            S_A[3] = 0.;
        } else if (pol_frac > 1.) {
            S_A[1] /= (pol_frac + 0.005);
            S_A[2] /= (pol_frac + 0.005);
            S_A[3] /= (pol_frac + 0.005);
        }

        *Iinv = S_A[0];
        *Iinv_pol = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]);

        // We have now updated the Stokes vector using plasma at current
        // position. Only do stuff below this line IF S_A[0] > 1.e-40. If
        // not, POLARIZATION_ACTIVE is set to FALSE and we reset S_A[i] = 0
        if (*Iinv_pol > 1.e-100) {
            stokes_to_f(f_u, f_tetrad_u, tetrad_u, S_A, Iinv, Iinv_pol);
            *tau += aI * (dl_current) * C;
            *tauF += fabs(rV) * (dl_current) * C;

            // Set POLARIZATION_ACTIVE to true; we are, after all,
            // in_volume.
            *POLARIZATION_ACTIVE = 1;

        } else {
            *POLARIZATION_ACTIVE = 0;
            S_A[1] = 0.;
            S_A[2] = 0.;
            S_A[3] = 0.;
        }
	}
    // Use sub-steps when the optical-depth step is only moderately stiff.
	else if (step_ratio <= MAX_ADAPTIVE_SUBSTEPS) {
        // if (*dl_current0/dl_opt>1.e5){
        //     printf ("N=%e x1=%.2e x2=%.2e x3=%.2e jI=%.2e jQ=%.2e jV=%.2e rQ=%.2e rV=%.2e aI=%.2e aQ=%.2e aV=%.2e C=%.2e\n",
        //             (*dl_current0/dl_opt),X_u[1],X_u[2],X_u[3],jI,jQ,jV,rQ,rV,aI,aQ,aV,C);
        // }
        // printf("[Adaptive Optical Depth stepsize] x1=%.2e x2=%.2e x3=%.2e dl_current=%.2e dl_opt=%.2e\n",
        //       X_u[1],X_u[2],X_u[3],*dl_current0,dl_opt);
        double dl_remain = *dl_current0;
        double dl_mark = 1.; // mark for updating radiation coefficients
        double dl_since_refresh = 0.;
        double param_magnitude_adaptive =
            fabs(aI) + fabs(aQ) + fabs(aV) + fabs(rQ) + fabs(rV);
        int update_point = 1;
        int tetrad_initialized = 0;
        int stokes_dirty = 0;

        while (dl_remain > 0.) {
            double complex S_before_step[4] = {S_A[0], S_A[1], S_A[2], S_A[3]};
            const char *solver_used = "trapezoid";
            double Iinv_before = *Iinv;
            double Iinv_pol_before = *Iinv_pol;
            // check whether update radiation coefficients and fluid parameters
            if (dl_mark > 0. && dl_remain <= dl_mark * (*dl_current0)) {
                dl_mark -= 0.1;
                update_point = 1;
            } else {
                update_point = 0;
            }

            if (update_point || !tetrad_initialized) {
                if (stokes_dirty && *POLARIZATION_ACTIVE && tetrad_initialized) {
                    stokes_to_f(f_u, f_tetrad_u, tetrad_u, S_A, Iinv, Iinv_pol);
                    stokes_dirty = 0;
                }
                if (tetrad_initialized) {
                    advance_polarized_step_position(X_u, k_u_old,
                                                    dl_since_refresh,
                                                    affine_scale);
                    dl_since_refresh = 0.;
                }

                create_observer_tetrad(X_u, k_u, modvar_local.U_u, modvar_local.B_u, tetrad_u);
                LOOP_ij if (isnan(tetrad_u[i][j])) return;
                create_tetrad_d(X_u, tetrad_u, tetrad_d);

                if (*POLARIZATION_ACTIVE) {
                    f_to_stokes(f_u, f_tetrad_u, tetrad_d, S_A, *Iinv, *Iinv_pol);
                }
                tetrad_initialized = 1;
            }

            double dl_to_refresh =
                distance_to_next_coeff_refresh(dl_remain, *dl_current0, dl_mark);
            double dl_eff = capped_adaptive_dl(dl_remain, dl_opt);
            dl_current = MIN(dl_to_refresh, dl_eff);

            // Given Stokes params and plasma coeffs, compute NEW Stokes params
            // after plasma step.
            if ((fabs(aI * dl_current * C) <= 0.3 &&
                 param_magnitude_adaptive * dl_current * C <= 0.5) ||
                fabs(aI) < 1e-20) {
                solver_used = "rk4";
                pol_rte_rk4_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV, dl_current, C, S_A);
            } else {
                pol_rte_trapezoid_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV, dl_current, C, S_A);
            }

            if (stokes_has_nonfinite(S_A)) {
                dump_bad_stokes_state("adaptive", solver_used,
                                      *POLARIZATION_ACTIVE, dl_current, C,
                                      S_before_step, S_A, Iinv_before,
                                      Iinv_pol_before, jI, jQ, jV, rQ, rV,
                                      aI, aQ, aV, X_u, modvar_local, nu_p,
                                      pitch_ang);
            }

            double pol_frac =
                sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]) /
                sqrt(S_A[0] * S_A[0]);

            if (isnan(pol_frac) || isinf(pol_frac)) {
                S_A[1] = 0.;
                S_A[2] = 0.;
                S_A[3] = 0.;
            } else if (pol_frac > 1.) {
                S_A[1] /= (pol_frac + 0.005);
                S_A[2] /= (pol_frac + 0.005);
                S_A[3] /= (pol_frac + 0.005);
            }

            *Iinv = S_A[0];
            *Iinv_pol = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]);

            if (*Iinv_pol > 1.e-100) {
                *POLARIZATION_ACTIVE = 1;
                *tau += aI * (dl_current) * C;
                *tauF += fabs(rV) * (dl_current) * C;
                stokes_dirty = 1;
            } else {
                *POLARIZATION_ACTIVE = 0;
                S_A[1] = 0.;
                S_A[2] = 0.;
                S_A[3] = 0.;
                stokes_dirty = 0;
            }

            dl_since_refresh += dl_current;

            if (update_point) {
                advance_polarized_step_position(X_u, k_u_old, dl_since_refresh,
                                                affine_scale);
                dl_since_refresh = 0.;
                get_fluid_params(X_u, &modvar_local);
                pitch_ang = pitch_angle(X_u, k_u_old, modvar_local.B_u, modvar_local.U_u);
                lower_index(X_u, k_u, k_d);
                nu_p = freq_in_plasma_frame(modvar_local.U_u, k_d);

                #if (EMISUSER)
                evaluate_coeffs_user(&jI, &jQ, &jU, &jV, &rQ, &rU, &rV, &aI,
                                     &aQ, &aU, &aV, nu_p, modvar_local, pitch_ang);
                #else
                evaluate_coeffs_single(&jI, &jQ, &jU, &jV, &rQ, &rU, &rV, &aI,
                                       &aQ, &aU, &aV, nu_p, modvar_local, pitch_ang);
                #endif

                dl_opt = check_stiffness(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV, C);
                param_magnitude_adaptive =
                    fabs(aI) + fabs(aQ) + fabs(aV) + fabs(rQ) + fabs(rV);
            }
            dl_remain -= dl_current;
        }
        advance_polarized_step_position(X_u, k_u_old, dl_since_refresh,
                                        affine_scale);
        if (stokes_dirty && *POLARIZATION_ACTIVE) {
            stokes_to_f(f_u, f_tetrad_u, tetrad_u, S_A, Iinv, Iinv_pol);
        }
	}
    // For extremely stiff segments, propagate the whole cell with one
    // constant-coefficient analytic step.
    else {
        dl_current = *dl_current0;
        double complex S_before_step[4] = {S_A[0], S_A[1], S_A[2], S_A[3]};
        const char *solver_used = "analytic";
        double Iinv_before = *Iinv;
        double Iinv_pol_before = *Iinv_pol;

        create_observer_tetrad(X_u, k_u, modvar_local.U_u, modvar_local.B_u, tetrad_u);
        LOOP_ij if (isnan(tetrad_u[i][j])) return;

        create_tetrad_d(X_u, tetrad_u, tetrad_d);

        if (*POLARIZATION_ACTIVE) {
            f_to_stokes(f_u, f_tetrad_u, tetrad_d, S_A, *Iinv, *Iinv_pol);
        }

        if (!pol_rte_analytic_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU,
                                   aV, dl_current, C, S_A)) {
            double param_magnitude =
                fabs(aI) + fabs(aQ) + fabs(aV) + fabs(rQ) + fabs(rV);
            if ((fabs(aI * dl_current * C) <= 0.5 &&
                 param_magnitude * dl_current * C <= 1.0) ||
                fabs(aI) < 1e-20) {
                solver_used = "rk4";
                pol_rte_rk4_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU, aV,
                                 dl_current, C, S_A);
            } else {
                solver_used = "trapezoid";
                pol_rte_trapezoid_step(jI, jQ, jU, jV, rQ, rU, rV, aI, aQ, aU,
                                       aV, dl_current, C, S_A);
            }
        }

        if (stokes_has_nonfinite(S_A)) {
            dump_bad_stokes_state("stiff", solver_used, *POLARIZATION_ACTIVE,
                                  dl_current, C, S_before_step, S_A,
                                  Iinv_before, Iinv_pol_before, jI, jQ, jV,
                                  rQ, rV, aI, aQ, aV, X_u, modvar_local, nu_p,
                                  pitch_ang);
        }

        double pol_frac =
            sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]) /
            sqrt(S_A[0] * S_A[0]);

        if (isnan(pol_frac) || isinf(pol_frac)) {
            S_A[1] = 0.;
            S_A[2] = 0.;
            S_A[3] = 0.;
        } else if (pol_frac > 1.) {
            S_A[1] /= (pol_frac + 0.005);
            S_A[2] /= (pol_frac + 0.005);
            S_A[3] /= (pol_frac + 0.005);
        }

        *Iinv = S_A[0];
        *Iinv_pol = sqrt(S_A[1] * S_A[1] + S_A[2] * S_A[2] + S_A[3] * S_A[3]);

        if (*Iinv_pol > 1.e-100) {
            stokes_to_f(f_u, f_tetrad_u, tetrad_u, S_A, Iinv, Iinv_pol);
            *tau += aI * dl_current * C;
            *tauF += fabs(rV) * dl_current * C;
            *POLARIZATION_ACTIVE = 1;
        } else {
            *POLARIZATION_ACTIVE = 0;
            S_A[1] = 0.;
            S_A[2] = 0.;
            S_A[3] = 0.;
        }
	}
    LOOP_i X_u[i] = X_u_old[i];
    //     printf("n_e=%g theta_e=%g B=%g\n",modvar_local.n_e,modvar_local.theta_e,modvar_local.B);
    //     printf("jI=%g jQ=%g jU=%g jV=%g\n",jI,jQ,jU,jV);
    //     printf("aI=%g aQ=%g aU=%g aV=%g\n",aI,aQ,aU,aV);
    //     printf("rQ=%g rU=%g rV=%g\n",rQ,rU,rV);
    //     printf("S_A=%g %g %g %g\n",S_A[0],S_A[1],S_A[2],S_A[3]);
    //     exit(1);
    // }
    
    // FILE *fp = fopen("diagnosis.txt", "a");
    // fprintf(fp,"r=%g X=%g %g %g dl=%g n_e=%g theta_e=%g b=%g S_A=%g %g %g %g act=%d\n",
    //         sqrt(X_u[1]*X_u[1]+X_u[2]*X_u[2]),X_u[1],X_u[2],X_u[3],*dl_current0/affine_scale,modvar_local.n_e,modvar_local.theta_e,modvar_local.B,
    //         S_A[0],S_A[1],S_A[2],S_A[3],*POLARIZATION_ACTIVE);
    // fclose(fp);
    if (!isfinite(creal(S_A[0])) || creal(S_A[0])<0. || creal(S_A[0])>1.e42){
        printf("r=%.2g X=%.2g %.2g %.2g n_e=%.2g theta_e=%.2g b=%.2g nu_p=%.2g pitch_ang=%.2g dl=%.2g dlC=%.2g\n",
                sqrt(X_u[1]*X_u[1]+X_u[2]*X_u[2]),X_u[1],X_u[2],X_u[3],
                modvar_local.n_e,modvar_local.theta_e,modvar_local.B,nu_p,pitch_ang,
                dl_current,dl_current*C);
        printf("jI=%.2e jQ=%.2e jV=%.2e rQ=%.2e rV=%.2e aI=%.2e aQ=%.2e aV=%.2e S_A=%.2g %.2g %.2g %.2g \n",
                jI,jQ,jV,rQ,rV,aI,aQ,aV,creal(S_A[0]),creal(S_A[1]),creal(S_A[2]),creal(S_A[3]));
        exit(1);
    }
}

void construct_f_obs_tetrad_u(double *X_u, double *k_u, double complex *f_u,
                              double complex *f_obs_tetrad_u) {

    double cam_up_u[4] = {0., 0., 0., -1.};
    double U_obs_u[4] = {0., 0., 0., 0.};
    double obs_tetrad_u[4][4], obs_tetrad_d[4][4];
    LOOP_ij obs_tetrad_u[i][j] = 0.;
    LOOP_ij obs_tetrad_d[i][j] = 0.;

    construct_U_vector(X_u, U_obs_u);
    create_observer_tetrad(X_u, k_u, U_obs_u, cam_up_u, obs_tetrad_u);
    create_tetrad_d(X_u, obs_tetrad_u, obs_tetrad_d);

    // Convert f_u to f_obs_tetrad_u
    LOOP_i f_obs_tetrad_u[i] = 0.;
    LOOP_ij f_obs_tetrad_u[i] += obs_tetrad_d[j][i] * f_u[j];
}

void radiative_transfer_polarized(double *lightpath, int steps,
                                  double frequency, double *f_x, double *f_y,
                                  double *p, int PRINT_POLAR, double *IQUV,
                                  double *tau, double *tauF) {
    int path_counter;
    double dl_current;

    double X_u[4], k_u[4], k_d[4];

    double Iinv, Iinv_pol;
    int POLARIZATION_ACTIVE = 0;

    double tetrad_u[4][4], tetrad_d[4][4];
    LOOP_ij tetrad_u[i][j] = 0.;
    LOOP_ij tetrad_d[i][j] = 0.;

    double photon_u_current[8] = {0., 0., 0., 0., 0., 0., 0., 0.};
    double complex f_tetrad_u[4] = {0., 0., 0., 0.};
    double complex f_u[4] = {0., 0., 0., 0.};
    double complex S_A[4] = {0., 0., 0., 0.};

    struct GRMHD modvar;
    modvar.B = 0;
    modvar.n_e = 0.;
    modvar.theta_e = 0;

    LOOP_i {
        modvar.B_u[i] = 0;
        modvar.U_u[i] = 0;
        modvar.B_d[i] = 0;
        modvar.U_d[i] = 0;
    }
    modvar.igrid_c = -1;

    // Move backward along constructed lightpath
    for (path_counter = steps - 1; path_counter > 0; path_counter--) {
        // Current position, wave vector, and dlambda
        LOOP_i {
            X_u[i] = lightpath[path_counter * 9 + i];
            k_u[i] = lightpath[path_counter * 9 + 4 + i];
        }
        dl_current = fabs(lightpath[(path_counter - 1) * 9 + 8]);

        // check normalization of k vectors.
        if (fabs(four_velocity_norm(X_u, k_u)) > 1e-6 && sqrt(X_u[1] * X_u[1] + X_u[2] * X_u[2] ) > 2.)
            normalize_null(X_u, k_u);

        // PLASMA INTEGRATION STEP
        //////////////////////////

        double r_current = sqrt(X_u[1] * X_u[1] + X_u[2] * X_u[2] );

        // Check whether the ray is currently in the GRMHD simulation volume
        if (get_fluid_params(X_u, &modvar) && r_current < RT_OUTER_CUTOFF) {
            pol_integration_step(modvar, frequency, &dl_current, C_CONST, X_u,
                                 k_u, k_d, &POLARIZATION_ACTIVE, f_u,
                                 f_tetrad_u, tetrad_d, tetrad_u, S_A, &Iinv,
                                 &Iinv_pol, tau, tauF);
        } // End of if(IN_VOLUME)
        // else {
        //     FILE *fp = fopen("diagnosis.txt", "a");
        //     fprintf(fp,"r=%g X=%g %g %g dl=%g S_A=%g act=%d\n",
        //             sqrt(X_u[1]*X_u[1]+X_u[2]*X_u[2]),X_u[1],X_u[2],X_u[3],dl_current,S_A[0],POLARIZATION_ACTIVE);
        //     fclose(fp);
        // }

        // SPACETIME-INTEGRATION STEP
        /////////////////////////////

        // If we HAVE been in-volume before, transport f_u (which is now
        // defined) one step. The final time this is done will be when
        // path_counter = 1; dl_current will then be at index 0 (path_counter -
        // 1).
        if (POLARIZATION_ACTIVE && path_counter > 0) {
            // Obtain the right k-vector, pointing back to observer, and
            // associated position. Pop into photon_u_current.
            LOOP_i {
                photon_u_current[i] = X_u[i];
                photon_u_current[i + 4] = k_u[i];
            }

            // One step: parallel transport of polarization vector.
            rk4_step_f(photon_u_current, f_u, dl_current);
        }
    } // End of for(path_counter...

    // CONSTRUCT FINAL (NON-INVARIANT) STOKES PARAMS SEEN BY OBSERVER
    /////////////////////////////////////////////////////////////////

    // Construct the observer tetrad.
    // X_u_current and k_u_current are simply the initial position and wave
    // vector. Note that k_u_current points INTO the camera sensor plane.
    LOOP_i {
        X_u[i] = lightpath[i];
        k_u[i] = lightpath[4 + i];
    }

    double complex f_obs_tetrad_u[4] = {0., 0., 0., 0.};
    construct_f_obs_tetrad_u(X_u, k_u, f_u, f_obs_tetrad_u);

    LOOP_i IQUV[i] = 0.;
    if (POLARIZATION_ACTIVE) {
        f_tetrad_to_stokes(Iinv, Iinv_pol, f_obs_tetrad_u, S_A);

        // Construct final (NON-INVARIANT) Stokes params.
        LOOP_i {
            // printf("S_A[%d]=%g ",i,S_A[i]);
            IQUV[i] = S_A[i] * pow(frequency, 3.);
        }
    }
}
