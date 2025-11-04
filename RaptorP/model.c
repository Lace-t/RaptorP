/*
 * model file for HARM3D data
 *
 * Please note that most of the code for the harm3d model was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 * Modifications were made in August 2022 by J. Davelaar
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"
#include "pluto.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

// GLOBAL VARS
//////////////
double ****p;

int N1, N2, N3;

double gam;

double R0, Rin, Rout, a, hslope;
double startx[4], stopx[4];
double *x1l,*x1r,*x2l,*x2r,*x3l,*x3r;
double dx1min=1.e40,dx1max=0.,dx2min=1.e40,dx2max=0.,dx3min=1.e40,dx3max=0.;

double L_unit, T_unit;
double RHO_unit, U_unit, B_unit;
double Ne_unit, Thetae_unit;

// FUNCTIONS
////////////

void init_model() {

    set_units(M_UNIT);

    fprintf(stderr, "\nStarting read in of PLUTO RMHD data...\n");

    init_rmhd_data(RMHD_FILE);
    //init_axis_data(RMHD_FILE);//z-axisymmetric
    //init_trace_data(RMHD_FILE);//include tracer
    //init_axis_trace_data(RMHD_FILE);//z-axisymmetry and include tracer
}

void init_rmhd_data(char *fname) {
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
    x3l=(double *)malloc(N3*sizeof(double));
    x3r=(double *)malloc(N3*sizeof(double));
    for(i=0;i<N3;i++) {
        re=fscanf(fp,"%d  %lf %lf\n", &token, &x3l[i], &x3r[i]);
        if(dx3min>(x3r[i]-x3l[i])) dx3min=x3r[i]-x3l[i];
        if(dx3max<(x3r[i]-x3l[i])) dx3max=x3r[i]-x3l[i];
        }
    fclose(fp);

    //perform translation to save imagesize & shift image center
    devi=(x3r[N3-1]+x3l[0])/2.;
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
    init_storage();
    fp=fopen(fname, "rb");
    double *buffer = malloc(N1*N2*N3*sizeof(double));

    if (fp == NULL) {
        fprintf(stderr, "\nCan't open sim data file... Abort!\n");
        exit(1);
    } else {
        fprintf(stderr, "\nSuccessfully opened %s. \n\nReading", fname);
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
    //calculate Internal energy
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++)p[UU][i][j][k]=p[UU][i][j][k]/(gam-1.);
        }
    }
}

int get_fluid_params(double X[NDIM], struct GRMHD *modvar) {

    int i, j, k;
    double del[NDIM];
    double rho, uu;
    double Bp[NDIM], V_u[NDIM], Vfac, VdotV, UdotBp;
    double g_uu[NDIM][NDIM], g_dd[NDIM][NDIM], coeff[4];
    double bsq, beta, beta_trans, b2, trat, Th_unit, two_temp_gam;

    if (X[1] < startx[1] || X[1] > stopx[1] || X[2] < startx[2] ||
        X[2] > stopx[2]) {
        (*modvar).n_e = 0.;
        return 0;
    }//beyond boundary return 0

    Xtoijk(X, &i, &j, &k, del);

    metric_uu(X, g_uu);
    metric_dd(X, g_dd);

    coeff[1] = del[1];
    coeff[2] = del[2];
    coeff[3] = del[3];

    rho = interp_scalar(p[KRHO], i, j, k, coeff);
    if(fabs(rho+1)<1e-4){
        printf("X= %e %e %e\n",X[1],X[2],X[3]);
        exit(1);
    }
    uu = interp_scalar(p[UU], i, j, k, coeff);
    (*modvar).n_e = rho * Ne_unit + 1e-40;

    Bp[1] = interp_scalar(p[B1], i, j, k, coeff);
    Bp[2] = interp_scalar(p[B2], i, j, k, coeff);
    Bp[3] = interp_scalar(p[B3], i, j, k, coeff);

    V_u[1] = interp_scalar(p[U1], i, j, k, coeff);
    V_u[2] = interp_scalar(p[U2], i, j, k, coeff);
    V_u[3] = interp_scalar(p[U3], i, j, k, coeff);

    VdotV = 0.;
    for (i = 1; i < NDIM; i++)
        for (j = 1; j < NDIM; j++)
            VdotV += g_dd[i][j] * V_u[i] * V_u[j];
    Vfac = sqrt(-1. / g_uu[0][0] * (1. + fabs(VdotV)));
    (*modvar).U_u[0] = -Vfac * g_uu[0][0];
    for (i = 1; i < NDIM; i++)
        (*modvar).U_u[i] = V_u[i] - Vfac * g_uu[0][i];

    lower_index(X, (*modvar).U_u, (*modvar).U_d);

    double Utot = 0;
    for (int i = 0; i < NDIM; i++)
        Utot += (*modvar).U_u[i] * (*modvar).U_d[i];

    UdotBp = 0.;
    for (i = 1; i < NDIM; i++)
        UdotBp += (*modvar).U_d[i] * Bp[i];
    (*modvar).B_u[0] = UdotBp;
    for (i = 1; i < NDIM; i++)
        (*modvar).B_u[i] =
            (Bp[i] + (*modvar).U_u[i] * UdotBp) / (*modvar).U_u[0];

    lower_index(X, (*modvar).B_u, (*modvar).B_d);

    bsq = (*modvar).B_u[0] * (*modvar).B_d[0] +
          (*modvar).B_u[1] * (*modvar).B_d[1] +
          (*modvar).B_u[2] * (*modvar).B_d[2] +
          (*modvar).B_u[3] * (*modvar).B_d[3];

    (*modvar).B = sqrt(bsq) * B_unit + 1e-40;

    (*modvar).beta = uu * (gam - 1.) / (0.5 * (bsq + 1.e-40));
    (*modvar).sigma = bsq / (rho + 1.e-40);

    beta = (*modvar).beta;

    beta_trans = 1.;
    b2 = pow(beta / beta_trans, 2);
    double Rhigh = R_HIGH;
    double Rlow = R_LOW;

    trat = Rhigh * b2 / (1. + b2) + Rlow / (1. + b2);

    Thetae_unit = 1. / 3. * (MPoME) / (trat + 1);

    (*modvar).theta_e = (uu / rho) * Thetae_unit;

//    beta = uu * (gam - 1.) / 0.5 / bsq;
//
//    beta_trans = 1.;
//    b2 = pow(beta / beta_trans, 2);
//
//    trat = 3.;
//    two_temp_gam = 0.5 * ((1. + 2. / 3. * (trat + 1.) / (trat + 2.)) + gam);
//    Th_unit = (1.4444444444 - 1.) * (PROTON_MASS / ELECTRON_MASS) / (1. + trat);
//
//    (*modvar).theta_e =
//        (2. / 15.) * (uu / rho) * (PROTON_MASS / ELECTRON_MASS) + 1e-40;

#if (DEBUG)
    if (uu < 0)
        fprintf(stderr, "U %e %e\n", uu, p[UU][i][j][k]);
    ;

    if ((*modvar).theta_e < 0)
        fprintf(stderr, "Te %e\n", (*modvar).theta_e);
    if ((*modvar).B < 0)
        fprintf(stderr, "B %e\n", (*modvar).B);
    if ((*modvar).n_e < 0)
        fprintf(stderr, "Ne %e %e\n", (*modvar).n_e, p[KRHO][i][j][k]);
#endif

//    if (bsq / rho > 1. || exp(X[1]) > 50.) {
//        (*modvar).n_e = 0;
//        return 0;
//    }

    return 1;
}

void set_units(double M_unit_) {

    L_unit = GGRAV * MBH / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    T_unit = L_unit / SPEED_OF_LIGHT;

    RHO_unit = M_unit_ / pow(L_unit, 3);
    U_unit = RHO_unit * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    B_unit = SPEED_OF_LIGHT * sqrt(4. * M_PI * RHO_unit);

    Ne_unit = RHO_unit / (PROTON_MASS + ELECTRON_MASS);
}

double interp_scalar(double ***var, int i, int j, int k, double coeff[4]) {

    double interp;
    int ip1, jp1, kp1;
    double b1, b2, b3, del[NDIM];

    del[1] = coeff[1];
    del[2] = coeff[2];
    del[3] = coeff[3];
    if (del[1] > 1+1e-4 || del[2] > 1+1e-4 || del[3] > 1+1e-4 || del[1] < -1e-4 || del[2] < -1e-4 || del[3] < -1e-4){
        printf("Over boundary: i=%d j=%d k=%d del[1]=%.4e del[2]=%.4e del[3]=%.4e\n",i,j,k,del[1],del[2],del[3]);
        return -1;
        }

    ip1 = i + 1;
    jp1 = j + 1;
    kp1 = k + 1;

    b1 = 1. - del[1];
    b2 = 1. - del[2];
    b3 = 1. - del[3];
    //turn off interpolation outside the boundary
    if (k==-1)return var[0][0][0];
    else if (k==N3)return var[0][0][N3-1];

    interp = var[i][j][k] * b1 * b2 + var[i][jp1][k] * b1 * del[2] +
             var[ip1][j][k] * del[1] * b2 + var[ip1][jp1][k] * del[1] * del[2];

    /* Now interpolate above in x3 */
    interp = b3 * interp + del[3] * (var[i][j][kp1] * b1 * b2 +
                                     var[i][jp1][kp1] * b1 * del[2] +
                                     var[ip1][j][kp1] * del[1] * b2 +
                                     var[ip1][jp1][kp1] * del[1] * del[2]);

    return interp;
}

int binary_search(double arr[],double arr2[], int low,int high, double target) {
	if (target<arr[0])return -1;
	if (target>arr2[high])return high;
    while (low <= high) {
        int mid = low+(high - low)/2;
        if ((arr[mid]+arr2[mid])/2 < target) low = mid + 1;
        else if ((arr[mid]+arr2[mid])/2 > target) high = mid - 1;
        else return (mid-1);
    }
    if(target>(arr[high]+arr2[high])/2) return high;
    else return (high-1);//smaller than middle, move to the left grid
}

void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]) {
    int low,high;
  	if (fabs(dx1max-dx1min)<1.e-4) *i=(int)floor((X[1]-x1l[0])/dx1min-0.5);
    else {// predict coordinate
      	low=(int)floor((X[1]-x1l[0])/dx1max-0.5);
      	high=(int)floor((X[1]-x1l[0])/dx1min+0.5);
      	*i = binary_search(x1l,x1r,low,MIN(high,N1-1),X[1]);
        }
  	if (fabs(dx2max-dx2min)<1.e-4) *j=(int)floor((X[2]-x2l[0])/dx2min-0.5);
    else {
      	low=(int)floor((X[2]-x2l[0])/dx2max-0.5);
      	high=(int)floor((X[2]-x2l[0])/dx2min+0.5);
      	*j = binary_search(x2l,x2r,low,MIN(high,N2-1),X[2]);
        }
  	if (fabs(dx3max-dx3min)<1.e-4) *k=(int)floor((X[3]-x3l[0])/dx3min-0.5);
    else {
      	low=(int)floor((X[3]-x3l[0])/dx3max-0.5);
      	high=(int)floor((X[3]-x3l[0])/dx3min+0.5);
      	*k = binary_search(x3l,x3r,low,MIN(high,N3-1),X[3]);
      	//if ((X[3]-(x3l[*k]+x3r[*k])/2)<0) {printf("%d %d %d\n",low,high,*k);}
        }

    if (*i < 0) {
        *i = 0;
        del[1] = 0.;
    } else if (*i > N1 - 2) {
        *i = N1 - 2;
        del[1] = 1.;
    } else {
        del[1] = (X[1]-(x1l[*i]+x1r[*i])/2)/((x1r[*i]-x1l[*i]+x1r[*i+1]-x1l[*i+1])/2);
    }

    if (*j < 0) {
        *j = 0;
        del[2] = 0.;
    } else if (*j > N2 - 2) {
        *j = N2 - 2;
        del[2] = 1.;
    } else {
        del[2] = (X[2]-(x2l[*j]+x2r[*j])/2) /((x2r[*j]-x2l[*j]+x2r[*j+1]-x2l[*j+1])/2);
    }

    if (*k < 0) {
        *k = -1;//xf: as a flag
        del[3] = 0.;
    } else if (*k > N3 - 2) {
        *k = N3 ;//xf: as a flag
        del[3] = 1.;
    } else {
        del[3] = (X[3]-(x3l[*k]+x3r[*k])/2) /((x3r[*k]-x3l[*k]+x3r[*k+1]-x3l[*k+1])/2);
    }
}

static void *malloc_rank1(int n1, int alloc_size) {
    void *A;

    if ((A = (void *)malloc(n1 * alloc_size)) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank1\n");
        exit(123);
    }

    return A;
}

double ***malloc_rank3(int n1, int n2, int n3) {
    double ***A;
    double *space;
    int i, j;

    space =(double*) malloc_rank1(n1 * n2 * n3, sizeof(double));

    A =(double***)  malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++) {
        A[i] =(double**) malloc_rank1(n2, sizeof(double *));
        for (j = 0; j < n2; j++) {
            A[i][j] = &(space[n3 * (j + n2 * i)]);
        }
    }

    return A;
}

void init_storage(void) {
    p = (double ****)malloc(NPRIM * sizeof(double ***));

    for (int n = 0; n < NPRIM; n++) {
        double *data_block = (double *)malloc(N1 * N2 * N3 * sizeof(double));
        p[n] = (double ***)malloc(N1 * sizeof(double **));

        for (int i = 0; i < N1; i++) {
            p[n][i] = (double **)malloc(N2 * sizeof(double *));

            for (int j = 0; j < N2; j++) {
                p[n][i][j] = &data_block[i*N2*N3 + j*N3];
                }
            }
        }
}

void compute_spec_user(struct Camera *intensityfield,
                       double energy_spectrum[num_frequencies][nspec]) {

    return;
}
