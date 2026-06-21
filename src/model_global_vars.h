/*
 * model_global_vars.h
 *
 * Please note that most of the code for the harm3d model was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 * Modifications were made in August 2022 by J. Davelaar
 */

#ifndef MODEL_GLOBAL_VARS_H
#define MODEL_GLOBAL_VARS_H

extern double ****p;

extern int N1, N2, N3;

extern double R_HIGH, R_LOW, gam;

extern double R0, Rin, Rout, a, hslope,P;
extern double startx[NDIM], stopx[NDIM];
extern double *x1l,*x1r,*x2l,*x2r,*x3l,*x3r;
extern double dx1min,dx1max,dx2min,dx2max,dx3min,dx3max;

extern double L_unit, T_unit;
extern double RHO_unit, U_unit, B_unit;
extern double Ne_unit, Thetae_unit;

#endif // MODEL_GLOBAL_VARS_H
