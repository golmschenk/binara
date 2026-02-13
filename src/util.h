#ifndef UTIL_H_
#define UTIL_H_
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

extern const double PI;
extern const double SQRT_2PI;
extern const double G; // cgs
extern const double C;
extern const double MSUN;
extern const double RSUN;
extern const double SEC_DAY;
extern const double BIG_NUM;


#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))


struct bounds
{
    double lo;
    double hi;
};

struct gauss_bounds
{
    double mean;
    double sigma;
};

typedef struct bounds bounds;
typedef struct gauss_bounds gauss_bounds;

void Free_2d(double** arr, int size);
void Free_3d(double*** arr, int size1, int size2);
double Gaussian(double x, double mean, double sigma);
int exists(const char* fname);

// Loading data
void Load_MCMC_Constants(
    int* py_niter, int* py_nchains, int* py_npars,
    int* py_nsectors, int* py_npast, double* py_dtemp, long int* buffer_size);
void Load_MCMC_Parameter_Info(
    const int NPARS, long int* buffer_size,
    bounds* limits, bounds* limited, gauss_bounds* gauss_pars,
    double* X_init, double* sigma);
void Load_MCMC_Sector_Points(
    const int NSECTORS, long int* buffer_size, long int* points_per_sector,
    long int* py_npoints);
void Load_MCMC_Data_Arrays(
    const int NPOINTS, long int* buffer_size, double* times, double* fluxes,
    double* errors, double* magdata, double* magerr);

#endif
