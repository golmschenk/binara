#ifndef _UTILHD_342
  #define _UTILHD_342 (1)
  #include "stdio.h"
  #include "stdlib.h"
  #include "math.h"
  #include "string.h"
  #include <unistd.h>

  #define PI 3.14159265358979323846
  #define SQRT_2PI 2.5066282746
  #define G 6.6743e-8 //cgs
  #define C 2.99792e10
  #define AU 1.495979e13
  #define MSUN 1.9885e33
  #define RSUN 6.955e10
  #define SEC_DAY 86400.0
  #define BIG_NUM 1.e30
  #define MAXSECTORS 100
  #define BUF_SIZE 65536

  #define NTAB 32
  #define IM1 2147483563
  #define IM2 2147483399
  #define AM (1.0/IM1)
  #define IMM1 (IM1-1)
  #define IA1 40014
  #define IA2 40692
  #define IQ1 53668
  #define IQ2 52774
  #define IR1 12211
  #define IR2 3791
  #define NDIV (1+IMM1/NTAB)
  #define eps 1.2e-7
  #define RNMX (1.0-eps)

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


  // structure to hold random number values
  struct RNG_Vars
  {
    long idum2; // some seed
    long iy;
    long iv[NTAB];
    int iset;
    double gset;
    // Number of times the rng is called
    long cts;
  };

  typedef struct bounds bounds;
  typedef struct gauss_bounds gauss_bounds;
  typedef struct RNG_Vars RNG_Vars;

  // Key variables that will be used in the MCMC_WRAPPER
  int points_per_sector[MAXSECTORS];

  void Set_Limits_Intialize_Proposals(bounds limited[], bounds limits[], 
    gauss_bounds gauss_pars[], double *sigma, int NSECTORS);
  double Ran2_Parallel(long *idum, RNG_Vars* state);
  double Gasdev2_Parallel(long *idum, RNG_Vars *state);
  void Free_1d(double *arr);
  void Free_2d(double **arr, int size);
  void Free_3d(double ***arr, int size1, int size2);
  double Gaussian(double x, double mean, double sigma);
  int exists(const char *fname);
  
  // Loading data
  void Get_Datafile_Name(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                      char path[]);
  int Load_MCMC_Constants(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                      int *py_niter, int *py_nchains, int *py_npars, 
                      int *py_nsectors, int *py_npast, double *py_dtemp, long int *buffer_size);
  int Load_MCMC_Parameter_Info(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                      const int NPARS,  long int *buffer_size,
                      bounds *limits, bounds *limited, gauss_bounds *gauss_pars, 
                      double *X_init, double *sigma);
  int Load_MCMC_Sector_Points(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                      const int NSECTORS, long int *buffer_size, long int *points_per_sector,
                      long int *py_npoints);
  int Load_MCMC_Data_Arrays(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                      const int NPOINTS, long int *buffer_size, double *times, double *fluxes,
                      double *errors, double *magdata, double *magerr);
  int count_lines(char* fname, long int *line_size);
#endif