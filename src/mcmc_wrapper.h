#ifndef _MCMC_WRAPPER_H_
#define _MCMC_WRAPPER_H_
#include "likelihood.h"
#include <omp.h>
#include "random_generator.h"

#define ENABLE_OPENMP (1)

extern const int npars_common, npars_unique;

void Gaussian_Proposal(double* x, double* sigma, double scale, double temp, double* y, const int NPARS,
                       RandomGenerator** random_generators_for_chains, int chain_number);

void Differential_Evolution_Proposal(double* x, double** history, double* y, const int NPARS, const int NPAST,
                                     const double GAMMA, RandomGenerator** random_generators_for_chains,
                                     int chain_number);

void Ptmcmc(int* index, double temp[], double logL[], double logP[], const int NCHAINS,
            RandomGenerator* random_generator);

#ifdef __cplusplus
extern "C" {
#endif

void Run_MCMC(int, int);

#ifdef __cplusplus
}
#endif

void Log_Data(int iter, double** x, double* logLx, int* index,
              long int* points_per_sector, double all_sector_phases[], double all_sector_fluxes[],
              double all_sector_uncertainties[], const int NPARS, const int NSECTORS, const int NCHAINS);

void Read_Parameters(double** X, const int NPARS, const int NCHAINS);
#endif
