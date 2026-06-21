/*
 * model file for reading non-AMR PLUTO data directly
 *
 * Written by Xufan Hu 2025
 */
#include <complex.h>
#include </home/titan/miniconda3/include/hdf5.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>

#define KRHO 0
#define UU 1
#define U1 2
#define U2 3
#define U3 4
#define B1 5
#define B2 6
#define B3 7

double ****p;

int N1, N2, N3;

double gam;

double R0, Rin, Rout, a, hslope, P;
double startx[4], stopx[4];
double *x1l,*x1r,*x2l,*x2r,*x3l,*x3r;

void init_grmhd_data(char *fname);
int binary_search(double arr[],double arr2[], int n, double target);

int main(){
    //printf("\nStarting read in of PLUTO RMHD data...\n");
    init_grmhd_data("../df/data.0096.dbl");
    //printf("%d",(int)(int)floor((4.947963+5)/0.02-0.5));
}

void init_grmhd_data(char *fname){
    int i,j,k,m,n,success;
    double devi;
    char line[256];
    char *token;
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
    fgets(line,512,fp);
    if (line[0] != '#') success = 1;
    }
    fsetpos(fp, &file_pos);

    //read the left and right sides of x1, x2, x3
    fscanf(fp,"%d \n",&(N1));
    x1l=(double *)malloc(N1*sizeof(double));
    x1r=(double *)malloc(N1*sizeof(double));
    for(i=0;i<N1;i++) fscanf(fp,"%d  %lf %lf\n", &token, &x1l[i], &x1r[i]);

    fscanf(fp,"%d \n",&(N2));
    x2l=(double *)malloc(N2*sizeof(double));
    x2r=(double *)malloc(N2*sizeof(double));
    for(i=0;i<N2;i++) fscanf(fp,"%d  %lf %lf\n", &token, &x2l[i], &x2r[i]);

    fscanf(fp,"%d \n",&(N3));
    x3l=(double *)malloc(N3*sizeof(double));
    x3r=(double *)malloc(N3*sizeof(double));
    for(i=0;i<N3;i++) fscanf(fp,"%d  %lf %lf\n", &token, &x3l[i], &x3r[i]);

    fclose(fp);

    //perform translation to save imagesize
    devi=(x3r[N3-1]+x3l[0])/2;
    for(i=0;i<N3;i++){
        x3l[i]-=devi;
        x3r[i]-=devi;
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

    printf("%d %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n",N1,N2,N3,gam,a,startx[1],startx[2],startx[3],
           stopx[1],stopx[2],stopx[3]);
    double x=5.295567e2;
    int index;
//    index=binary_search(x3l,x3r,N3,x);
//    printf("i= %d x= %lf xl= %lf xr= %lf",index,x,x3l[index],x3r[index]);

//    //read .dbl file
    fp=fopen(fname, "rb");
    double val;

    if (fp == NULL) {
        fprintf(stderr, "\nCan't open sim data file... Abort!\n");
        exit(1);
    } else {
        fprintf(stderr, "\nSuccessfully opened %s. \n\nReading", fname);
    }

    p = (double ****)malloc(8* sizeof(double ***));
    double *buffer = malloc(N1*N2*N3*sizeof(double));

    for (int n = 0; n < 8; n++) {
        double *data_block = (double *)malloc(N1 * N2 * N3 * sizeof(double));
        p[n] = (double ***)malloc(N1 * sizeof(double **));

        for (int i = 0; i < N1; i++) {
            p[n][i] = (double **)malloc(N2 * sizeof(double *));

            for (int j = 0; j < N2; j++) {
                p[n][i][j] = &data_block[i*N2*N3 + j*N3];
                }
            }
        }

    for(m=0;m<8;m++){
        if(m==0)n=0;//perform order exchange
        else if(m==7)n=1;
        else n=m+1;

        size_t elements_read = fread(buffer, sizeof(double), N3*N2*N1, fp);

        #pragma omp parallel for collapse(2)
        for (int i = 0; i < N1; i++) {
            for (int j = 0; j < N2; j++) {
                double ***p_n = p[n];
                for (int k = 0; k < N3; k++) {
                    size_t src_index = k*N2*N1 + j*N1 + i;
                    p_n[i][j][k] = buffer[src_index];
                    }
                }
            }
        fprintf(stderr,".");
        }
    free(buffer);
    fclose(fp);


//    //calculate Internal energy
//    for(i=0;i<N1;i++){
//        for(j=0;j<N2;j++){
//            for(k=0;k<N3;k++)p[UU][i][j][k]=p[UU][i][j][k]/(gam-1.);
//        }
//    }
//
    i=250;j=200;k=250;
    printf("\n %e %e %e %e %e %e %e %e\n", p[KRHO][i][j][k],
                       p[UU][i][j][k], p[U1][i][j][k], p[U2][i][j][k],
                       p[U3][i][j][k], p[B1][i][j][k], p[B2][i][j][k],
                       p[B3][i][j][k]);
}

int binary_search(double arr[],double arr2[], int n, double target) {
    int low = 0, high = n - 1;
	if (target<arr[0])return -1;
	if (target>arr[n-1])return n;
    while (low <= high) {
        int mid = low+(high - low)/2;
        if (arr[mid] < target) low = mid + 1;
        else if (arr[mid] > target) high = mid - 1;
        else return (mid-1);
    }
    if(target>(arr[high]+arr2[high])/2) return high;
    else return (high-1);//smaller than middle, move to the left grid
}
