#ifndef _MCMC_WRAPPER_H_
#define _MCMC_WRAPPER_H_
#include "likelihood.h"
#include <omp.h>
#include "random_generator.h"

#define ENABLE_OPENMP (1)

extern const int npars_common, npars_unique;

void Gaussian_Proposal(double* x, double* sigma, double scale, double temp, double* y, int NPARS,
                       RandomGenerator** random_generators_for_chains, int chain_number);

void Differential_Evolution_Proposal(double* x, double** history, double* y, int NPARS, int NPAST,
                                     double GAMMA, RandomGenerator** random_generators_for_chains,
                                     int chain_number);

void Ptmcmc(int* index, double temp[], double logL[], double logP[], int NCHAINS,
            RandomGenerator* random_generator);

void Run_MCMC(int, int);

void Log_Data(int iter, double** x, double* logLx, int* index,
              long int* points_per_sector, double all_sector_phases[], double all_sector_fluxes[],
              double all_sector_uncertainties[], int NPARS, int NSECTORS, int NCHAINS, double&
              best_recorded_log_likelihood);

void Read_Parameters(double** X, int NPARS, int NCHAINS);
#endif
