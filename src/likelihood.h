#ifndef _LIKELIHOOD_H_
#define _LIKELIHOOD_H_
#include "eclipse_method_enum.h"
#include "util.h"

extern const int npars_common, npars_unique;
void Trajectory(double* times, double* traj_pars, double* d_arr,
                double* Z1_arr, double* Z2_arr, double* rr_arr,
                double* ff_arr, int Nt);
void Calculate_Lightcurve(double* times, size_t Nt, double* pars, double* template_, EclipseMethod eclipse_method);
double calculate_log_likelihood(const double all_sector_phases[], const double all_sector_fluxes[],
                                const double all_sector_uncertainties[], const long int points_per_sector[],
                                int NSECTORS, const double all_parameters[],
                                const double mag_data[], const double mag_err[]);
void Swap(double* a, double* b);
double Log_Prior(int NPARS, double* parameter_values, gauss_bounds* gauss_pars);
#endif
