#include "string.h"
#include "math.h"
#include "stdlib.h"
#include "stdio.h"
#include "util.h"

extern const int npars_common = 19;
extern const int npars_unique = 3;


const double PI = 3.14159265358979323846;
const double SQRT_2PI = 2.5066282746;
const double G = 6.6743e-8; //cgs
const double C = 2.99792e10;
const double MSUN = 1.9885e33;
const double RSUN = 6.955e10;
const double SEC_DAY = 86400.0;
const double BIG_NUM = 1.e30;

void Free_2d(double** arr, int size)
{
    int i;
    for (i = 0; i < size; i++) { free(arr[i]); }
    free(arr);
}

void Free_3d(double*** arr, int size1, int size2)
{
    int i, j;
    for (i = 0; i < size1; i++)
    {
        for (j = 0; j < size2; j++) { free(arr[i][j]); }
        free(arr[i]);
    }
    free(arr);
}

// Standard gaussian
double Gaussian(double x, double mean, double sigma)
{
    return (1 / sigma / SQRT_2PI) * exp(-pow((x - mean) / sigma, 2.) / 2.);
}


// Function to check if (parameter) file exists
int exists(const char* fname)
{
    FILE* file = fopen(fname, "r");
    if (file != NULL)
    {
        fclose(file);
        return 1; // File exists
    }
    return 0; // File does not exist
}


void Get_Datafile_Name(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                       char path[])
{
    char fname[256] = "";
    if (secular_drift_flag == 0)
    {

        sprintf(fname, "%d_sector_%d_run_%d.txt", tic, sector, run_id);
        strcat(path, fname);
    }
    else
    {
        sprintf(fname, "%d_sector_%d_run_%d_drift.txt", tic, sector, run_id);
        strcat(path, fname);
    }
    return;
}


int Load_MCMC_Constants(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                        int* py_niter, int* py_nchains, int* py_npars,
                        int* py_nsectors, int* py_npast, double* py_dtemp, long int* buffer_size)
{
    char path[1024] = "data/py_initialize/";;
    Get_Datafile_Name(tic, sector, run_id, secular_drift_flag, path);
    printf("Reading constants \n");

    if (exists(path) != 1)
    {
        printf("ERROR: Data file does not exist: %s\n", path);
        return 0;
    }

    FILE* data_file = fopen(path, "r");
    *buffer_size = 0;
    int temp_int;
    double temp_dbl;

    // The header of the file contains NITER, NCHAINS, NPARS, NSECTORS and dtemp
    fscanf(data_file, "%d\t", &temp_int);
    (*py_niter) = temp_int;
    fscanf(data_file, "%d\t", &temp_int);
    (*py_nchains) = temp_int;
    fscanf(data_file, "%d\t", &temp_int);
    (*py_npars) = temp_int;
    fscanf(data_file, "%d\t", &temp_int);
    (*py_nsectors) = temp_int;
    fscanf(data_file, "%d\t", &temp_int);
    (*py_npast) = temp_int;
    fscanf(data_file, "%lf\n", &temp_dbl);
    (*py_dtemp) = temp_dbl;

    printf("Read the following input parameters: \n");
    printf("NITER: %d NCHAINS: %d NPARS: %d NSECTORS: %d NPAST: %d dtemp: %f \n", *py_niter,
           *py_nchains, *py_npars, *py_nsectors, *py_npast, *py_dtemp);

    // Get the number of bytes read
    *buffer_size = ftell(data_file);
    fclose(data_file);
    return 1;
}

int Load_MCMC_Parameter_Info(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                             const int NPARS, long int* buffer_size,
                             bounds* limits, bounds* limited, gauss_bounds* gauss_pars,
                             double* X_init, double* sigma)
{
    char path[1024] = "data/py_initialize/";;
    Get_Datafile_Name(tic, sector, run_id, secular_drift_flag, path);
    printf("Reading parameter information \n");

    if (exists(path) != 1)
    {
        printf("ERROR: Data file does not exist: %s\n", path);
        return 0;
    }

    FILE* data_file = fopen(path, "r");
    fseek(data_file, *buffer_size, SEEK_SET);

    double par_min, par_max, par_mean, par_jump, prior_gauss_mean, prior_gauss_std,
           bc_buffer;
    int boundary_condition;

    for (int ipar = 0; ipar < NPARS; ipar++)
    {
        fscanf(data_file, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &par_mean, &par_min, &par_max,
               &bc_buffer, &par_jump, &prior_gauss_mean, &prior_gauss_std);
        boundary_condition = (int)bc_buffer;
        limits[ipar].lo = par_min;
        limits[ipar].hi = par_max;
        limited[ipar].lo = boundary_condition;
        limited[ipar].hi = limited[ipar].lo;
        gauss_pars[ipar].mean = prior_gauss_mean;
        gauss_pars[ipar].sigma = prior_gauss_std;
        X_init[ipar] = par_mean;
        sigma[ipar] = par_jump;
    }

    printf("Read parameter arrays \n");
    *buffer_size = ftell(data_file);
    fclose(data_file);
    return 1;
}

int Load_MCMC_Sector_Points(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                            const int NSECTORS, long int* buffer_size, long int* points_per_sector,
                            long int* py_npoints)
{
    char path[1024] = "data/py_initialize/";;
    Get_Datafile_Name(tic, sector, run_id, secular_drift_flag, path);
    printf("Reading sector information \n");

    if (exists(path) != 1)
    {
        printf("ERROR: Data file does not exist: %s\n", path);
        return 0;
    }

    FILE* data_file = fopen(path, "r");
    fseek(data_file, *buffer_size, SEEK_SET);

    // Read points per sector
    int temp_int;
    *py_npoints = 0;

    for (int i = 0; i < NSECTORS - 1; i++)
    {
        fscanf(data_file, "%d\t", &temp_int);
        points_per_sector[i] = temp_int;
        *py_npoints += temp_int;
    }
    fscanf(data_file, "%d\n", &temp_int);
    points_per_sector[NSECTORS - 1] = temp_int;
    *py_npoints += temp_int;

    printf("Data has %ld points \n", *py_npoints);
    *buffer_size = ftell(data_file);
    fclose(data_file);
    return 1;
}

int Load_MCMC_Data_Arrays(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                          const int NPOINTS, long int* buffer_size, double* times, double* fluxes,
                          double* errors, double* magdata, double* magerr)
{
    char path[1024] = "data/py_initialize/";;
    Get_Datafile_Name(tic, sector, run_id, secular_drift_flag, path);
    printf("Reading lightcurve and color data \n");

    if (exists(path) != 1)
    {
        printf("ERROR: Data file does not exist: %s\n", path);
        return 0;
    }

    FILE* data_file = fopen(path, "r");
    fseek(data_file, *buffer_size, SEEK_SET);

    double time, flux, err;
    for (int i = 0; i < NPOINTS; i++)
    {
        fscanf(data_file, "%lf\t%lf\t%lf\n", &time, &flux, &err);
        times[i] = time;
        fluxes[i] = flux;
        errors[i] = err;
    }

    double dist, gmag, vb, bg, gt;
    fscanf(data_file, "%lf\t%lf\t%lf\t%lf\t%lf\n", &dist, &gmag, &vb, &bg, &gt);
    double gmag_err, vb_err, bg_err, gt_err;
    fscanf(data_file, "%lf\t%lf\t%lf\t%lf", &gmag_err, &vb_err, &bg_err, &gt_err);

    magdata[0] = dist;
    magdata[1] = gmag;
    magdata[2] = vb;
    magdata[3] = bg;
    magdata[4] = gt;

    magerr[0] = gmag_err;
    magerr[1] = vb_err;
    magerr[2] = bg_err;
    magerr[3] = gt_err;

    printf("Distance to the source is %lf pc\n", dist);

    *buffer_size = ftell(data_file);
    fclose(data_file);
    return 1;
}
