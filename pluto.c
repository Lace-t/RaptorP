#define _FILE_OFFSET_BITS 64
#define _POSIX_C_SOURCE 200809L

/*
 * model file for performing changes for non-AMR PLUTO data
 *
 * Written by Xufan Hu 2025
 */
#include <complex.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <omp.h>
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
    fp = fopen("../grid.out","r");
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
        fprintf(stderr,".");
    }
    fclose(fp);
    //calculate Internal energy
    double inv_gam_minus1 = 1.0/(gam-1.0);
    #pragma omp parallel for collapse(3) schedule(static)
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<2*N3;k++){
                p[UU][i][j][k] = p[UU][i][j][k] * inv_gam_minus1;
            }
        }
    }
    N3*=2;//update current N3
}

//use a tracer exclude the ambient environment
// The location of the tracer is not fixed in .dbl. Default: 'trc' follows 'prs', so jump=0. 
// You may need to check dbl.out to set jump correctly.
void init_trace_data(char *fname,int jump){
    int i,j,k,m,n,success;
    double devi;
    char line[256];
	char *signal;
    int token,re;
    fpos_t file_pos;
    FILE *fp;

    //reading grid.out
    fp = fopen("../grid.out","r");
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
        fprintf(stderr,".");
    }
    //load tracer
    double *buffer = malloc(N1*N2*N3*sizeof(double));
    double ***trace;
    trace = (double ***)malloc(N1 * sizeof(double **));
    for (j = 0; j < N1; j++) {
        trace[j] = (double **)malloc(N2 * sizeof(double *));
        for (k = 0; k < N2; k++) {
            trace[j][k] = (double *)malloc(N3 * sizeof(double));
        }
    }
    //skip useless data
    for (i=0;i<jump;i++){
        size_t elements_read = fread(buffer, sizeof(double), N3*N2*N1, fp);
    }

    for(k=0;k<N3;k++){
        for(j=0;j<N2;j++){
            for(i=0;i<N1;i++){
                if (fread (&val,sizeof(double), 1, fp) != 1){
                    printf ("! InputDataReadSlice(): cannot find tracer!\n");
                    break;
                  }
                  //if (val<0.5)val=0.;
                  trace[i][j][k]=val;
            }
        }
    }

    fclose(fp);

    //calculate Internal energy
    double inv_gam_minus1_t = 1.0/(gam-1.0);
    #pragma omp parallel for collapse(3) schedule(static)
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
                p[UU][i][j][k] = p[UU][i][j][k] * inv_gam_minus1_t;
                if(k>2){ //exclude the ambient environment, keep some grids to avoid interpolation error
                    p[KRHO][i][j][k] *= trace[i][j][k];
                    p[UU][i][j][k] *= trace[i][j][k];
                }

            }
        }
    }
    free(trace);
}

void init_axis_trace_data(char *fname,int jump){
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
        fprintf(stderr,".");
    }
    //load tracer
    double *buffer = malloc(N1*N2*N3*sizeof(double));
    double ***trace;
    trace = (double ***)malloc(N1 * sizeof(double **));
    for (j = 0; j < N1; j++) {
        trace[j] = (double **)malloc(N2 * sizeof(double *));
        for (k = 0; k < N2; k++) {
            trace[j][k] = (double *)malloc(N3 * sizeof(double));
        }
    }
    //skip useless data
     for (i=0;i<jump;i++){
        size_t elements_read = fread(buffer, sizeof(double), N3*N2*N1, fp);
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
    double inv_gam_minus1_a = 1.0/(gam-1.0);
    #pragma omp parallel for collapse(3) schedule(static)
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
                p[UU][i][j][k+N3] = p[UU][i][j][k+N3] * inv_gam_minus1_a * trace[i][j][k];
                p[UU][i][j][N3-k-1] = p[UU][i][j][N3-k-1] * inv_gam_minus1_a * trace[i][j][k];
                p[KRHO][i][j][k+N3] *= trace[i][j][k];
                p[KRHO][i][j][N3-k-1] *= trace[i][j][k];
            }
        }
    }
    N3*=2;//update current N3
}

/* Locate the snapshot number in the basename: data.<digits>... */
static int get_data_number(const char *fname, int *number,
                           const char **digits, const char **digits_end) {
    const char *basename = strrchr(fname, '/');
    const char *marker;
    const char *end;
    char *parse_end;
    long value;

    basename = basename == NULL ? fname : basename + 1;
    marker = strstr(basename, "data.");
    if (marker == NULL) return 0;

    marker += strlen("data.");
    if (!isdigit((unsigned char)*marker)) return 0;

    errno = 0;
    value = strtol(marker, &parse_end, 10);
    end = parse_end;
    if (errno == ERANGE || value > INT_MAX || value < INT_MIN || end == marker)
        return 0;

    if (number != NULL) *number = (int)value;
    if (digits != NULL) *digits = marker;
    if (digits_end != NULL) *digits_end = end;
    return 1;
}

// Replace the number following data. while preserving its zero-padded width.
void replace_number(const char *fname, int new_num, char *newname) {
    const char *digits;
    const char *digits_end;
    int width;
    int written;

    if (!get_data_number(fname, NULL, &digits, &digits_end)) {
        fprintf(stderr, "Invalid PLUTO data filename '%s': expected data.<number>\n", fname);
        exit(1);
    }

    width = (int)(digits_end - digits);
    written = snprintf(newname, 256, "%.*s%0*d%s", (int)(digits - fname),
                       fname, width, new_num, digits_end);
    if (written < 0 || written >= 256) {
        fprintf(stderr, "PLUTO data filename is too long: '%s'\n", fname);
        exit(1);
    }
}


/* One maximal run of adjacent y-columns that comes from the same snapshot.
   j_end is exclusive, so the file contains
   (j_end - j_begin) * N1 consecutive doubles for this range. */
typedef struct {
    int k;
    int j_begin;
    int j_end;  /* exclusive */
} RetardRange;

/* All ranges assigned to one snapshot.  Keeping a separate plan per snapshot
   lets us open that snapshot once and consume every required block. */
typedef struct {
    RetardRange *ranges;
    size_t count;
    size_t capacity;
} RetardSnapshotPlan;

/* Grow a snapshot's range list as the (k,j) map is compressed. */
static void append_retard_range(RetardSnapshotPlan *plan, int k,
                                int j_begin, int j_end) {
    if (plan->count == plan->capacity) {
        size_t new_capacity = plan->capacity == 0 ? 64 : 2 * plan->capacity;
        RetardRange *new_ranges = realloc(plan->ranges,
                                          new_capacity * sizeof(*new_ranges));
        if (new_ranges == NULL) {
            fprintf(stderr, "Cannot allocate slow-light read plan\n");
            exit(1);
        }
        plan->ranges = new_ranges;
        plan->capacity = new_capacity;
    }

    plan->ranges[plan->count].k = k;
    plan->ranges[plan->count].j_begin = j_begin;
    plan->ranges[plan->count].j_end = j_end;
    plan->count++;
}

static void read_retard_block(FILE *fp, const char *fname, int variable,
                              const RetardRange *range, double *buffer) {
    /* PLUTO's on-disk order is [variable][k][j][i].  All N1 x-cells
       belonging to adjacent j-columns are therefore one contiguous block. */
    off_t cell_offset =
        (off_t)variable * (off_t)N1 * (off_t)N2 * (off_t)N3 +
        (off_t)range->k * (off_t)N1 * (off_t)N2 +
        (off_t)range->j_begin * (off_t)N1;
    /* off_t/fseeko keep offsets valid for snapshots larger than 2 GiB. */
    off_t byte_offset = cell_offset * (off_t)sizeof(double);
    size_t value_count =
        (size_t)(range->j_end - range->j_begin) * (size_t)N1;

    if (fseeko(fp, byte_offset, SEEK_SET) != 0) {
        fprintf(stderr,
                "Seek error in %s (variable=%d, k=%d, j=[%d,%d)): %s\n",
                fname, variable, range->k, range->j_begin, range->j_end,
                strerror(errno));
        exit(1);
    }
    if (fread(buffer, sizeof(double), value_count, fp) != value_count) {
        fprintf(stderr,
                "Read error in %s (variable=%d, k=%d, j=[%d,%d))\n",
                fname, variable, range->k, range->j_begin, range->j_end);
        exit(1);
    }
}

static void read_retard_snapshot(const char *fname,
                                 const RetardSnapshotPlan *plan,
                                 int jump, double ***trace) {
    FILE *fp;
    double *buffer;
    size_t max_value_count = 0;
    size_t r;
    int m, n, i, j;

    /* A single reusable buffer is sized for this snapshot's largest range;
       no complete variable or snapshot is held a second time in memory. */
    for (r = 0; r < plan->count; r++) {
        size_t value_count =
            (size_t)(plan->ranges[r].j_end - plan->ranges[r].j_begin) *
            (size_t)N1;
        if (value_count > max_value_count) max_value_count = value_count;
    }
    if (max_value_count == 0) return;

    buffer = malloc(max_value_count * sizeof(*buffer));
    if (buffer == NULL) {
        fprintf(stderr, "Cannot allocate slow-light read buffer for %s\n", fname);
        exit(1);
    }

    fp = fopen(fname, "rb");
    if (fp == NULL) {
        fprintf(stderr, "%s not found: %s\n", fname, strerror(errno));
        exit(1);
    }

    /* PLUTO stores complete variables consecutively.  Keeping the variable
       loop outside the range loop makes file offsets increase monotonically
       and avoids jumping between multi-gigabyte variable planes per column. */
    for (m = 0; m < 8; m++) {
        if (m == 0) n = 0;
        else if (m == 7) n = 1;
        else n = m + 1;

        for (r = 0; r < plan->count; r++) {
            const RetardRange *range = &plan->ranges[r];
            read_retard_block(fp, fname, m, range, buffer);
            /* The file block is [j][i], whereas Raptor stores p[n][i][j][k].
               Scatter the contiguous read buffer into Raptor's layout. */
            for (j = range->j_begin; j < range->j_end; j++) {
                size_t row_offset =
                    (size_t)(j - range->j_begin) * (size_t)N1;
                for (i = 0; i < N1; i++) {
                    p[n][i][j][range->k] = buffer[row_offset + (size_t)i];
                }
            }
        }
    }

    /* A non-negative jump selects the optional tracer plane after the eight
       primitive variables.  It uses exactly the same spatial read plan. */
    if (jump >= 0) {
        int tracer_variable = 8 + jump;
        for (r = 0; r < plan->count; r++) {
            const RetardRange *range = &plan->ranges[r];
            read_retard_block(fp, fname, tracer_variable, range, buffer);
            for (j = range->j_begin; j < range->j_end; j++) {
                size_t row_offset =
                    (size_t)(j - range->j_begin) * (size_t)N1;
                for (i = 0; i < N1; i++) {
                    trace[i][j][range->k] = buffer[row_offset + (size_t)i];
                }
            }
        }
    }

    fclose(fp);
    free(buffer);
}

void init_retard_data(char *fname,int jump) {
    int i,j,k,success;
    double devi;
    char line[256];
	char* signal;
    int token,re;
    fpos_t file_pos;
    FILE *fp;

    // initialize slow light region
    if (azimuth != 270){
        fprintf(stderr,"To use slow-light, azimuth must be 270 degree !\n");
        exit(1);
    }
    int dn=2;//Difference in Data File Numbers
    double dt=2.; //time interval in code unit of data update 
    double region[4]={-3.,3.,40.,60.};//determine the dynamical region {ymin,ymax,zmin,zmax}

    double rectangle = atan((region[1]-region[0]) / (region[3]-region[2])); // radians
    double diagonal = sqrt(pow(region[1]-region[0], 2.) + pow(region[3]-region[2], 2.));
    /* Ensure INCLINATION (degrees) is converted to radians before comparing with rectangle */
    int num;
    if (INCLINATION<=90.){
        num = (int)( cos((INCLINATION * M_PI/180.0) - rectangle) * diagonal / dt ); // numbers of data in need
    }
    else {
        num = (int)( - cos((INCLINATION * M_PI/180.0) +rectangle) * diagonal / dt );
    }

    int start_num;
    if (!get_data_number(fname, &start_num, NULL, NULL)) {
        fprintf(stderr, "Invalid PLUTO data filename '%s': expected data.<number>\n", fname);
        exit(1);
    }
    fprintf(stderr,"\n\nRetard: up to %d files in need\n",num + 1);

    //reading grid.out
    fp = fopen("../grid.out","r");
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
    region[2]=region[2]-devi-SHIFT;
    region[3]=region[3]-devi-SHIFT;

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
    fprintf(stderr,"Begin with %s to load slow light data...\n",fname);
    char newname[256]={0};

    //initialize tracer
    double ***trace = NULL;
    if (jump >= 0) {
        trace = (double ***)malloc(N1 * sizeof(double **));
        for (j = 0; j < N1; j++) {
            trace[j] = (double **)malloc(N2 * sizeof(double *));
            for (k = 0; k < N2; k++) {
                trace[j][k] = (double *)malloc(N3 * sizeof(double));
            }
        }
    }

    /* First assign every (k,j) column to a snapshot.  Then compress each
       fixed-k row into maximal, contiguous j intervals for that snapshot.
       The map costs N2*N3 integers (about 2 MiB for a 500x1000 yz grid). */
    if (num < 0) {
        fprintf(stderr, "Invalid negative slow-light snapshot index: %d\n", num);
        exit(1);
    }
    size_t map_size = (size_t)N2 * (size_t)N3;
    int *snapshot_map = malloc(map_size * sizeof(*snapshot_map));
    RetardSnapshotPlan *plans =
        calloc((size_t)num + 1, sizeof(*plans));
    if (snapshot_map == NULL || plans == NULL) {
        fprintf(stderr, "Cannot allocate slow-light snapshot map\n");
        exit(1);
    }

    double inclination = INCLINATION * M_PI / 180.0;
    double cos_inclination = cos(inclination);
    double sin_inclination = sin(inclination);
    for (k = 0; k < N3; k++) {
        double z_center = (x3l[k] + x3r[k]) / 2.0;
        for (j = 0; j < N2; j++) {
            double y_center = (x2l[j] + x2r[j]) / 2.0;
            double timer;
            int snapshot;

            if (INCLINATION <= 90.0) {
                timer = cos_inclination * (z_center - region[2]) -
                        sin_inclination * (y_center - region[1]);
            } else {
                timer = cos_inclination * (z_center - region[3]) -
                        sin_inclination * (y_center - region[1]);
            }

            /* Preserve the original boundary rules exactly: non-positive
               delays use the starting file, and delays beyond the available
               window use the final snapshot (index num). */
            if (timer <= 0.0) {
                snapshot = 0;
            } else {
                snapshot = (int)floor(timer / dt + 1e-12);
                if (snapshot >= num) snapshot = num;
            }
            snapshot_map[(size_t)k * (size_t)N2 + (size_t)j] = snapshot;
        }
    }

    /* Run-length encode each k-row.  Because the retardation time is linear
       in y for fixed k, a snapshot normally contributes one j interval. */
    for (k = 0; k < N3; k++) {
        int j_begin = 0;
        int snapshot = snapshot_map[(size_t)k * (size_t)N2];
        for (j = 1; j <= N2; j++) {
            int next_snapshot =
                j < N2 ? snapshot_map[(size_t)k * (size_t)N2 + (size_t)j]
                       : -1;
            if (j == N2 || next_snapshot != snapshot) {
                append_retard_range(&plans[snapshot], k, j_begin, j);
                j_begin = j;
                snapshot = next_snapshot;
            }
        }
    }
    free(snapshot_map);

    size_t total_ranges = 0;
    int used_snapshots = 0;
    /* Snapshot indices are inclusive: 0..num.  Process one file at a time so
       only one descriptor is open and every required range is read before it
       is closed.  Empty snapshot plans are skipped. */
    for (int snapshot = 0; snapshot <= num; snapshot++) {
        if (plans[snapshot].count == 0) continue;
        used_snapshots++;
        total_ranges += plans[snapshot].count;

        if (snapshot == 0) {
            int written = snprintf(newname, sizeof(newname), "%s", fname);
            if (written < 0 || written >= (int)sizeof(newname)) {
                fprintf(stderr, "PLUTO data filename is too long: '%s'\n", fname);
                exit(1);
            }
        } else {
            replace_number(fname, start_num + snapshot * dn, newname);
        }
        fprintf(stderr, "Retard: reading %s in %zu contiguous ranges\n",
                newname, plans[snapshot].count);
        read_retard_snapshot(newname, &plans[snapshot], jump, trace);
    }
    fprintf(stderr,
            "Retard: read %d snapshots using %zu contiguous ranges\n",
            used_snapshots, total_ranges);

    for (int snapshot = 0; snapshot <= num; snapshot++) {
        free(plans[snapshot].ranges);
    }
    free(plans);

    //calculate Internal energy & deal with trace
    #pragma omp parallel for collapse(3) schedule(static)
    for(i=0;i<N1;i++){
        for(j=0;j<N2;j++){
            for(k=0;k<N3;k++){
                p[UU][i][j][k]=p[UU][i][j][k]/(gam-1.);
                if (jump>=0){
                    p[KRHO][i][j][k]*=trace[i][j][k];
                    p[UU][i][j][k]*=trace[i][j][k];
                }
            }
        }
    }
    if (jump >= 0) {
        for (j = 0; j < N1; j++) {
            for (k = 0; k < N2; k++) free(trace[j][k]);
            free(trace[j]);
        }
        free(trace);
    }
}
