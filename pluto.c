/*
 * model file for performing changes for non-AMR PLUTO data
 *
 * Written by Xufan Hu 2025
 */
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

#include "pluto.h"

//Symmetric with respect to the plane at z=0
void init_axis_data(char *fname){
    int i,j,k,m,n,success;
    double devi;
    char line[256];
	char* signal;
    int token,re;
    fpos_t file_pos;
    FILE *fp;

    //reading grid.out
    fp = fopen("grid.out","r");
    if (fp == NULL){
        perror("grid.out not found !\n");
        exit(1);
    }

    //Move file pointer to the first line of grid.out that does not begin with a "#".
    success = 0;
    while(!success){
        fgetpos(fp, &file_pos);
        signal=fgets(line,256,fp);
        if (line[0] != '#') success = 1;
        }
    fsetpos(fp, &file_pos);

    //read the left and right sides of x1, x2, x3
    re=fscanf(fp,"%d \n",&(N1));
    x1l=(double *)malloc(N1*sizeof(double));
    x1r=(double *)malloc(N1*sizeof(double));
    for(i=0;i<N1;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x1l[i], &x1r[i]);
        if(dx1min>(x1r[i]-x1l[i])) dx1min=x1r[i]-x1l[i];
        if(dx1max<(x1r[i]-x1l[i])) dx1max=x1r[i]-x1l[i];
        }

    re=fscanf(fp,"%d \n",&(N2));
    x2l=(double *)malloc(N2*sizeof(double));
    x2r=(double *)malloc(N2*sizeof(double));
    for(i=0;i<N2;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x2l[i], &x2r[i]);
        if(dx2min>(x2r[i]-x2l[i])) dx2min=x2r[i]-x2l[i];
        if(dx2max<(x2r[i]-x2l[i])) dx2max=x2r[i]-x2l[i];
        }

    re=fscanf(fp,"%d \n",&(N3));
    x3l=(double *)malloc(2*N3*sizeof(double));
    x3r=(double *)malloc(2*N3*sizeof(double));
	for(i=0;i<N3;i++) {
		re=fscanf(fp,"%d  %lf %lf\n", &token, &x3l[N3+i], &x3r[N3+i]);
		if(dx3min>(x3r[N3+i]-x3l[N3+i])) dx3min=x3r[N3+i]-x3l[N3+i];
		if(dx3max<(x3r[N3+i]-x3l[N3+i])) dx3max=x3r[N3+i]-x3l[N3+i];
		}

    fclose(fp);

    //generate coordinates on the other side
    for(i=0;i<N3;i++){
        x3l[N3-i-1]=-x3r[N3+i];
        x3r[N3-i-1]=-x3l[N3+i];
    }
    //calculate header data
    gam=4./3.;
    a=0.;
    startx[1]=x1l[0];
    startx[2]=x2l[0];
    startx[3]=x3l[0];
    hslope=0.;
    Rin=0.;
    Rout=0.;
    stopx[0] = 1.;
    stopx[1] = x1r[N1-1];
    stopx[2] = x2r[N2-1];
    stopx[3] = x3r[2*N3-1];

    //read .dbl file
    fp=fopen(fname, "rb");
    double val;

    if (fp == NULL) {
        fprintf(stderr, "\nCan't open sim data file... Abort!\n");
        exit(1);
    } else {
        fprintf(stderr, "\nSuccessfully opened %s. \n\nReading", fname);
    }

    p = (double ****)malloc(8 * sizeof(double ***));
    for (i = 0; i < 8; i++) {
        p[i] = (double ***)malloc(N1 * sizeof(double **));
        for (j = 0; j < N1; j++) {
            p[i][j] = (double **)malloc(N2 * sizeof(double *));
            for (k = 0; k < N2; k++) {
                p[i][j][k] = (double *)malloc(2*N3 * sizeof(double));
            }
        }
    }

    for(m=0;m<8;m++){
        if(m==0)n=0;//perform order exchange
        else if(m==7)n=1;
        else n=m+1;
        for(k=0;k<N3;k++){
            for(j=0;j<N2;j++){
                for(i=0;i<N1;i++){
                    if (fread (&val,sizeof(double), 1, fp) != 1){
                        printf ("! InputDataReadSlice(): error reading data\n");
                        break;
                      }
                      if(n==4 || n==7)p[n][i][j][N3-k-1]=-val;
                      else p[n][i][j][N3-k-1]=val;
                      p[n][i][j][k+N3]=val;
                }
            }
        }
        printf(".");
    }
    fclose(fp);
    //calculate Internal energy
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<2*N3;k++)p[UU][i][j][k]=p[UU][i][j][k]/(gam-1.);
        }
    }
    N3*=2;//update current N3

//    printf("%d %d %d %lf %lf %lf %lf %lf %lf %lf\n",N1,N2,N3,gam,
//           startx[1],startx[2],startx[3],stopx[1],stopx[2],stopx[3]);

//    i=250;j=200;k=19;
//    printf("\n %e %e %e %e %e %e %e %e\n", p[KRHO][i][j][k],
//                       p[UU][i][j][k], p[U1][i][j][k], p[U2][i][j][k],
//                       p[U3][i][j][k], p[B1][i][j][k], p[B2][i][j][k],
//                       p[B3][i][j][k]);
//
//    i=250;j=200;k=1180;
//    printf("\n %e %e %e %e %e %e %e %e\n", p[KRHO][i][j][k],
//                       p[UU][i][j][k], p[U1][i][j][k], p[U2][i][j][k],
//                       p[U3][i][j][k], p[B1][i][j][k], p[B2][i][j][k],
//                       p[B3][i][j][k]);
}

//use a tracer exclude the ambient environment
void init_trace_data(char *fname){
    int i,j,k,m,n,success;
    double devi;
    char line[256];
	char *signal;
    int token,re;
    fpos_t file_pos;
    FILE *fp;

    //reading grid.out
    fp = fopen("grid.out","r");
    if (fp == NULL){
        perror("grid.out not found !\n");
        exit(1);
    }

    //Move file pointer to the first line of grid.out that does not begin with a "#".
    success = 0;
    while(!success){
        fgetpos(fp, &file_pos);
        signal=fgets(line,256,fp);
        if (line[0] != '#') success = 1;
        }
    fsetpos(fp, &file_pos);

    //read the left and right sides of x1, x2, x3
    re=fscanf(fp,"%d \n",&(N1));
    x1l=(double *)malloc(N1*sizeof(double));
    x1r=(double *)malloc(N1*sizeof(double));
    for(i=0;i<N1;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x1l[i], &x1r[i]);
        if(dx1min>(x1r[i]-x1l[i])) dx1min=x1r[i]-x1l[i];
        if(dx1max<(x1r[i]-x1l[i])) dx1max=x1r[i]-x1l[i];
        }

    re=fscanf(fp,"%d \n",&(N2));
    x2l=(double *)malloc(N2*sizeof(double));
    x2r=(double *)malloc(N2*sizeof(double));
    for(i=0;i<N2;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x2l[i], &x2r[i]);
        if(dx2min>(x2r[i]-x2l[i])) dx2min=x2r[i]-x2l[i];
        if(dx2max<(x2r[i]-x2l[i])) dx2max=x2r[i]-x2l[i];
        }

    re=fscanf(fp,"%d \n",&(N3));
    x3l=(double *)malloc(N3*sizeof(double));
    x3r=(double *)malloc(N3*sizeof(double));
    for(i=0;i<N3;i++) {
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x3l[i], &x3r[i]);
        if(dx3min>(x3r[i]-x3l[i])) dx3min=x3r[i]-x3l[i];
        if(dx3max<(x3r[i]-x3l[i])) dx3max=x3r[i]-x3l[i];
        }

    fclose(fp);

    //perform translation to save imagesize
    devi=(x3r[N3-1]+x3l[0])/2;
    for(i=0;i<N3;i++){
        x3l[i]=x3l[i]-devi-SHIFT;
        x3r[i]=x3r[i]-devi-SHIFT;
    }
    //calculate header data
    gam=4./3.;
    a=0.;
    startx[1]=x1l[0];
    startx[2]=x2l[0];
    startx[3]=x3l[0];
    hslope=0.;
    Rin=0.;
    Rout=0.;
    stopx[0] = 1.;
    stopx[1] = x1r[N1-1];
    stopx[2] = x2r[N2-1];
    stopx[3] = x3r[N3-1];

    //read .dbl file
    fp=fopen(fname, "rb");
    double val;

    if (fp == NULL) {
        fprintf(stderr, "\nCan't open sim data file... Abort!\n");
        exit(1);
    } else {
        fprintf(stderr, "\nSuccessfully opened %s. \n\nReading", fname);
    }
    //initialize p
    p = (double ****)malloc(8 * sizeof(double ***));
    for (i = 0; i < 8; i++) {
        p[i] = (double ***)malloc(N1 * sizeof(double **));
        for (j = 0; j < N1; j++) {
            p[i][j] = (double **)malloc(N2 * sizeof(double *));
            for (k = 0; k < N2; k++) {
                p[i][j][k] = (double *)malloc(N3 * sizeof(double));
            }
        }
    }

    for(m=0;m<8;m++){
        if(m==0)n=0;//perform order exchange
        else if(m==7)n=1;
        else n=m+1;
        for(k=0;k<N3;k++){
            for(j=0;j<N2;j++){
                for(i=0;i<N1;i++){
                    if (fread (&val,sizeof(double), 1, fp) != 1){
                        printf ("! InputDataReadSlice(): error reading data\n");
                        break;
                      }
                      p[n][i][j][k]=val;
                }
            }
        }
        printf(".");
    }
    //load tracer
    double ***trace;
    trace = (double ***)malloc(N1 * sizeof(double **));
    for (j = 0; j < N1; j++) {
        trace[j] = (double **)malloc(N2 * sizeof(double *));
        for (k = 0; k < N2; k++) {
            trace[j][k] = (double *)malloc(N3 * sizeof(double));
        }
    }

    for(k=0;k<N3;k++){
        for(j=0;j<N2;j++){
            for(i=0;i<N1;i++){
                if (fread (&val,sizeof(double), 1, fp) != 1){
                    printf ("! InputDataReadSlice(): cannot find tracer!\n");
                    break;
                  }
                  trace[i][j][k]=val;
            }
        }
    }

    fclose(fp);
//    printf("%lf",p[KRHO][15][215][200]);
//    printf(" %lf\n",trace[15][215][200]);
    //calculate Internal energy
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
                p[UU][i][j][k]/=(gam-1.);
                if(k>10){//exclude the ambient environment, keep some grids to avoid interpolation error
                    p[KRHO][i][j][k]*=trace[i][j][k];
                    p[UU][i][j][k]*=trace[i][j][k];
                }

            }
        }
    }
//    printf("%lf",p[KRHO][15][215][200]);
}

void init_axis_trace_data(char *fname){
    int i,j,k,m,n,success;
    double devi;
    char line[256];
	char *signal;
    int token,re;
    fpos_t file_pos;
    FILE *fp;

    //reading grid.out
    fp = fopen("grid.out","r");
    if (fp == NULL){
        perror("grid.out not found !\n");
        exit(1);
    }

    //Move file pointer to the first line of grid.out that does not begin with a "#".
    success = 0;
    while(!success){
        fgetpos(fp, &file_pos);
        signal=fgets(line,256,fp);
        if (line[0] != '#') success = 1;
        }
    fsetpos(fp, &file_pos);

    //read the left and right sides of x1, x2, x3
    re=fscanf(fp,"%d \n",&(N1));
    x1l=(double *)malloc(N1*sizeof(double));
    x1r=(double *)malloc(N1*sizeof(double));
    for(i=0;i<N1;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x1l[i], &x1r[i]);
        if(dx1min>(x1r[i]-x1l[i])) dx1min=x1r[i]-x1l[i];
        if(dx1max<(x1r[i]-x1l[i])) dx1max=x1r[i]-x1l[i];
        }

    re=fscanf(fp,"%d \n",&(N2));
    x2l=(double *)malloc(N2*sizeof(double));
    x2r=(double *)malloc(N2*sizeof(double));
    for(i=0;i<N2;i++){
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x2l[i], &x2r[i]);
        if(dx2min>(x2r[i]-x2l[i])) dx2min=x2r[i]-x2l[i];
        if(dx2max<(x2r[i]-x2l[i])) dx2max=x2r[i]-x2l[i];
        }

    re=fscanf(fp,"%d \n",&(N3));
    x3l=(double *)malloc(2*N3*sizeof(double));
    x3r=(double *)malloc(2*N3*sizeof(double));
	for(i=0;i<N3;i++) {
		re=fscanf(fp,"%d  %lf %lf\n", &token, &x3l[N3+i], &x3r[N3+i]);
		if(dx3min>(x3r[N3+i]-x3l[N3+i])) dx3min=x3r[N3+i]-x3l[N3+i];
		if(dx3max<(x3r[N3+i]-x3l[N3+i])) dx3max=x3r[N3+i]-x3l[N3+i];
		}

    fclose(fp);

    //generate coordinates on the other side
    for(i=0;i<N3;i++){
        x3l[N3-i-1]=-x3r[N3+i];
        x3r[N3-i-1]=-x3l[N3+i];
    }
    //calculate header data
    gam=4./3.;
    a=0.;
    startx[1]=x1l[0];
    startx[2]=x2l[0];
    startx[3]=x3l[0];
    hslope=0.;
    Rin=0.;
    Rout=0.;
    stopx[0] = 1.;
    stopx[1] = x1r[N1-1];
    stopx[2] = x2r[N2-1];
    stopx[3] = x3r[2*N3-1];
    //read .dbl file
    fp=fopen(fname, "rb");
    double val;

    if (fp == NULL) {
        fprintf(stderr, "\nCan't open sim data file... Abort!\n");
        exit(1);
    } else {
        fprintf(stderr, "\nSuccessfully opened %s. \n\nReading", fname);
    }

    p = (double ****)malloc(8 * sizeof(double ***));
    for (i = 0; i < 8; i++) {
        p[i] = (double ***)malloc(N1 * sizeof(double **));
        for (j = 0; j < N1; j++) {
            p[i][j] = (double **)malloc(N2 * sizeof(double *));
            for (k = 0; k < N2; k++) {
                p[i][j][k] = (double *)malloc(2*N3 * sizeof(double));
            }
        }
    }

    for(m=0;m<8;m++){
        if(m==0)n=0;//perform order exchange
        else if(m==7)n=1;
        else n=m+1;
        for(k=0;k<N3;k++){
            for(j=0;j<N2;j++){
                for(i=0;i<N1;i++){
                    if (fread (&val,sizeof(double), 1, fp) != 1){
                        printf ("! InputDataReadSlice(): error reading data\n");
                        break;
                      }
                      if(n==4 || n==7)p[n][i][j][N3-k-1]=-val;
                      else p[n][i][j][N3-k-1]=val;
                      p[n][i][j][k+N3]=val;
                }
            }
        }
        printf(".");
    }
    //load tracer
    double ***trace;
    trace = (double ***)malloc(N1 * sizeof(double **));
    for (j = 0; j < N1; j++) {
        trace[j] = (double **)malloc(N2 * sizeof(double *));
        for (k = 0; k < N2; k++) {
            trace[j][k] = (double *)malloc(N3 * sizeof(double));
        }
    }

    for(k=0;k<N3;k++){
        for(j=0;j<N2;j++){
            for(i=0;i<N1;i++){
                if (fread (&val,sizeof(double), 1, fp) != 1){
                    printf ("! InputDataReadSlice(): cannot find tracer!\n");
                    break;
                  }
                  trace[i][j][k]=val;
            }
        }
    }
    fclose(fp);
    //calculate Internal energy
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
                p[UU][i][j][k+N3]=p[UU][i][j][k+N3]/(gam-1.)*trace[i][j][k];
                p[UU][i][j][N3-k-1]=p[UU][i][j][N3-k-1]/(gam-1.)*trace[i][j][k];
                p[KRHO][i][j][k+N3]=p[KRHO][i][j][k+N3]*trace[i][j][k];
                p[KRHO][i][j][N3-k-1]=p[KRHO][i][j][N3-k-1]*trace[i][j][k];
            }
        }
    }
    N3*=2;//update current N3
}
