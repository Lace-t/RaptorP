/*
 * Radboud Polarized Integrator
 * Copyright 2014-2021 Black Hole Cam (ERC Synergy Grant)
 * Authors: Thomas Bronzwaer, Jordy Davelaar, Monika Moscibrodzka, Ziri Younsi
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

// GLOBAL VARS
//////////////

char RMHD_FILE[256];

double MBH, M_UNIT, TIME_INIT, INCLINATION,azimuth;
double R_HIGH, R_LOW, Tur_or_Rec, P, Amin, Risco, Router,SHIFT;
double FREQS_PER_DEC, FREQ_MIN, FREQ_MAX;

double SOURCE_DIST; // Distance to M87 (cm); for Sgr A* use (2.47e22)

int IMG_WIDTH, IMG_HEIGHT;
double CAM_SIZE_X, CAM_SIZE_Y;
double STEPSIZE,STEPSIZE_MAX,STEPSIZE_MIN;

// FUNCTIONS
////////////

// Read model parameters from model.in
void read_model(char *argv[]) {
    char temp[100], temp2[100];
    FILE *input;
    char inputfile[100];
	int re;

    sscanf(argv[1], "%s", inputfile);
    fprintf(stdout, "\nUsing model parameter file %s\n", inputfile);

    input = fopen(inputfile, "r");
    if (input == NULL) {
        fprintf(stderr, "Can't read file %s! Aborting", inputfile);
        exit(1);
    }

    // Model parameters
    re=fscanf(input, "%s %s %lf", temp, temp2, &MBH);
    re=fscanf(input, "%s %s %lf", temp, temp2, &SOURCE_DIST);
    re=fscanf(input, "%s %s %lf", temp, temp2, &M_UNIT);
    re=fscanf(input, "%s %s %lf", temp, temp2, &R_LOW);
    re=fscanf(input, "%s %s %lf", temp, temp2, &R_HIGH);
    //re=fscanf(input, "%s %s %lf", temp, temp2, &Tur_or_Rec);
    //re=fscanf(input, "%s %s %lf", temp, temp2, &P);
    //re=fscanf(input, "%s %s %lf", temp, temp2, &Amin);
    //re=fscanf(input, "%s %s %lf", temp, temp2, &Risco);
    //re=fscanf(input, "%s %s %lf", temp, temp2, &Router);//XF: useless
	Tur_or_Rec=1.;P=0.;Amin=0.;Risco=3.9882;Router=2.22174;
    re=fscanf(input, "%s %s %lf", temp, temp2, &INCLINATION);
    re=fscanf(input, "%s %s %lf", temp, temp2, &azimuth);

    // Observer parameters
    re=fscanf(input, "%s %s %d", temp, temp2, &IMG_WIDTH);
    re=fscanf(input, "%s %s %d", temp, temp2, &IMG_HEIGHT);
	re=fscanf(input, "%s %s %lf", temp, temp2, &SHIFT);//XF: Translation image center along the z-axis
    re=fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_X);
    re=fscanf(input, "%s %s %lf", temp, temp2, &CAM_SIZE_Y);

    re=fscanf(input, "%s %s %lf", temp, temp2, &FREQS_PER_DEC);
    re=fscanf(input, "%s %s %lf", temp, temp2, &FREQ_MIN);
    re=fscanf(input, "%s %s %lf", temp, temp2, &STEPSIZE);
    re=fscanf(input, "%s %s %lf", temp, temp2, &STEPSIZE_MAX);
    re=fscanf(input, "%s %s %lf", temp, temp2, &STEPSIZE_MIN);
    //re=fscanf(input, "%s %s %d", temp, temp2, &max_level);
    max_level=1;

    // Second argument: RMHD file
    sscanf(argv[2], "%s", RMHD_FILE);
    sscanf(argv[3], "%lf", &TIME_INIT);

    fprintf(stderr, "\nModel parameters:\n\n");
    fprintf(stderr, "MBH \t\t= %g Msun\n", MBH);
    fprintf(stderr, "DISTANCE \t= %g kpc\n", SOURCE_DIST);
    fprintf(stderr, "M_UNIT \t\t= %g grams\n", M_UNIT);
    fprintf(stderr, "R_LOW \t\t= %g \n", R_LOW);
    fprintf(stderr, "R_HIGH \t\t= %g \n", R_HIGH);
    //fprintf(stderr, "P \t\t= %g \n", P);
    fprintf(stderr, "INCLINATION \t= %g deg\n", INCLINATION);

    fprintf(stderr, "METRIC \t\t= ");
#if (metric == MKSBHAC)
    fprintf(stderr, "MKS BHAC\n");
#elif (metric == CKS)
    //fprintf(stderr, "CKS BHAC\n");
    fprintf(stderr, "CAR PLUTO\n");
#elif (metric == MKSHARM)
    fprintf(stderr, "MKS HARM3D\n");
#endif

    fprintf(stderr, "\nObserver parameters:\n\n");
    fprintf(stderr, "IMG_WIDTH \t= %d \n", IMG_WIDTH);
    fprintf(stderr, "IMG_HEIGHT \t= %d \n", IMG_HEIGHT);
    fprintf(stderr, "CAM_SIZE_X \t= %g GM/c2\n", CAM_SIZE_X);
    fprintf(stderr, "CAM_SIZE_Y \t= %g GM/c2\n", CAM_SIZE_Y);
    fprintf(stderr, "FREQS \t= %lf \n", FREQS_PER_DEC);
    fprintf(stderr, "FREQ \t= %g Hz\n", FREQ_MIN);
    fprintf(stderr, "STEPSIZE \t= %g \n", STEPSIZE);
    fprintf(stderr, "STEPSIZE_MAX \t= %g \n", STEPSIZE_MAX);
    fprintf(stderr, "STEPSIZE_MIN \t= %g \n", STEPSIZE_MIN);

    // to cgs units
    MBH *= MSUN;
    SOURCE_DIST *= KPCTOCM;

    fclose(input);
}

// For a single block this function will iterate over the pixels and call
// geodesic integrations as well as radiation transfer
void calculate_image_block(struct Camera *intensityfield,
                           double frequencies[num_frequencies]) {
#if (LIGHTPATH)
	FILE *fp = fopen("lightpaths.txt", "w");
    fprintf(fp,"# step alpha beta x1 x2 x3 R dlambda\n");
#endif

#pragma omp parallel for shared(frequencies, intensityfield, p)                \
    schedule(static, 1)
    for (int pixel = 0; pixel < tot_pixels; pixel++) {
        int steps = 0;

        double *lightpath2 = malloc(9 * max_steps * sizeof(double));

#if (POL)
        double f_x = 0.;
        double f_y = 0.;
        double p = 0.;
#endif
        // INTEGRATE THIS PIXEL'S GEODESIC
        integrate_geodesic((*intensityfield).alpha[pixel],
                           (*intensityfield).beta[pixel], lightpath2, &steps,
                           CUTOFF_INNER);

    #if (LIGHTPATH)// write x1,x2,x3,dlambda for analysis
		int thread_id = omp_get_thread_num();
        #pragma omp critical//to keep writing safe
		{
		for (int i = 0; i < steps; i++) {
            double r=sqrt(lightpath2[i*9+1]*lightpath2[i*9+1]+lightpath2[i*9+2]*lightpath2[i*9+2]);
		    fprintf(fp, "%d %.4g %.4g %.4g %.4g %.4g %.4g %.4g\n",i,(*intensityfield).alpha[pixel],(*intensityfield).beta[pixel],
                    lightpath2[i*9+1],lightpath2[i*9+2],lightpath2[i*9+3],r,lightpath2[i*9+8]);
			}
		}
    #endif
        // PERFORM RADIATIVE TRANSFER AT DESIRED FREQUENCIES, STORE RESULTS
#if (POL)
        for (int f = 0; f < num_frequencies; f++) {

            radiative_transfer_polarized(lightpath2, steps, frequencies[f],
                                         &f_x, &f_y, &p, 0,
                                         (*intensityfield).IQUV[pixel][f],
                                         &(*intensityfield).tau[pixel][f],
                                         &(*intensityfield).tauF[pixel][f]);
        }

#else
        radiative_transfer_unpolarized(lightpath2, steps, frequencies,
                                       (*intensityfield).IQUV[pixel],
                                       &(*intensityfield).tau[pixel]);
        for (int f = 0; f < num_frequencies; f++) {
            (*intensityfield).IQUV[pixel][f][0] *= pow(frequencies[f], 3.);
        }
#endif
        free(lightpath2);
    }
#pragma omp barrier

#if (LIGHTPATH)
    fclose(fp);
#endif
}

// Functions that computes a spectrum at every frequency
// by integrating over the image struct
void compute_spec(struct Camera *intensityfield,
                  double energy_spectrum[num_frequencies][nspec]) {
    double dA, S_I, S_Q, S_U, S_V;

    for (int block = 0; block < tot_blocks; block++) {
        dA = (intensityfield)[block].dx[0] * (intensityfield)[block].dx[1];
        for (int pixel = 0; pixel < tot_pixels; pixel++) {
            for (int freq = 0; freq < num_frequencies; freq++) {
#if (POL)

                S_I = (intensityfield)[block].IQUV[pixel][freq][0];
                S_Q = (intensityfield)[block].IQUV[pixel][freq][1];
                S_U = (intensityfield)[block].IQUV[pixel][freq][2];
                S_V = (intensityfield)[block].IQUV[pixel][freq][3];

                // Stokes I
                energy_spectrum[freq][0] += S_I * dA;

                // Stokes Q
                energy_spectrum[freq][1] += S_Q * dA;

                // Stokes U
                energy_spectrum[freq][2] += S_U * dA;

                // stokes V
                energy_spectrum[freq][3] += S_V * dA;

#else
                energy_spectrum[freq][0] +=
                    (intensityfield)[block].IQUV[pixel][freq][0] * dA;
#endif
            }
        }
    }
}
