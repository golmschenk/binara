#ifndef _LIKELIHOOD_H_
#define _LIKELIHOOD_H_
#include "util.h"

extern const int npars_common, npars_unique;
void Trajectory(double* times, double* traj_pars, double* d_arr,
                double* Z1_arr, double* Z2_arr, double* rr_arr,
                double* ff_arr, int Nt);
#ifdef __cplusplus
extern "C" {
#endif
void Calculate_Lightcurve(double* times, long Nt, double* pars, double* template_);
double Log_Likelihood(double all_sector_phases[], double all_sector_fluxes[],
                      double all_sector_uncertainties[], long int points_per_sector[],
                      const int NSECTORS, double all_parameters[],
                      double mag_data[], double mag_err[], const int gmag_flag,
                      const int color_flag, const int secular_drift_flag);
#ifdef __cplusplus
}
#endif
void Swap(double* a, double* b);

double Log_Prior(const int NPARS, double* parameter_values, gauss_bounds* gauss_pars);
#endif
