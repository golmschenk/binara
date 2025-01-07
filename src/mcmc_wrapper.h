#ifndef _MCMC_WRAPPER_H_
  #define _MCMC_WRAPPER_H_
  #include "likelihood.h"
  #include "omp.h"

  #define ENABLE_OPENMP (1)
  
  extern const int npars_common, npars_unique;

  void Uniform_Proposal(double *x, long *seed, bounds limits[], double *y, RNG_Vars *state, const int NPARS);

  void Gaussian_Proposal(double *x, long *seed, double *sigma, double scale, double temp, 
    double *y, RNG_Vars *state, const int NPARS);

  void Differential_Evolution_Proposal(double *x, long *seed, double **history, double *y, RNG_Vars *state, 
    const int NPARS, const int NPAST, const double GAMMA);

  void Ptmcmc(int *index, double temp[], double logL[], double logP[], const int NCHAINS);

  #ifdef __cplusplus
  extern "C" {
  #endif

    void Run_MCMC(int, int, int, int, int, int);

  #ifdef __cplusplus
  }
  #endif

  void Make_Files(const int tic, const int sector, const int run_id, const int gmag_flag, const int color_flag,
                  const int secular_drift_flag, char *chainname, char *outname, char *parfile);

  void Log_Data(char *chainname, char *outname, char *parname, int iter, double **x, double *logLx, int *index,
              long int *points_per_sector, double all_sector_phases[], double all_sector_fluxes[], 
              double all_sector_uncertainties[], const int NPARS, const int NSECTORS, const int NCHAINS, 
              const int secular_drift_flag);

  void Read_Parameters(char *chainname, double **X, const int NPARS, const int NCHAINS);
#endif