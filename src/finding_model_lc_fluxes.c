#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"

void Calculate_Lightcurve(double *times, long Nt, double *pars, double *template_);

int main(void) {
    FILE *fp = fopen("data/pars/pars.110602878_sector_34_gmag_OMP_1.out", "r");
    double pars[21];
    for (int i = 0; i < 21; i++) {
        fscanf(fp, "%lf", &pars[i]);
    }
    fclose(fp);

    FILE *ft = fopen("data/lightcurves/folded_lightcurves/110602878_sector_34.txt", "r");
    //FILE *ft = fopen("data/lightcurves/original_folded_lightcurves/110602878_sector_34_binned150.txt", "r");
    long Nt;
    double P_days_file;
    fscanf(ft, "%ld %lf", &Nt, &P_days_file);

    double *t_days= malloc(sizeof(double) * Nt);
    double *flux = malloc(sizeof(double) * Nt);
    double *err  =malloc(sizeof(double) * Nt);

    for (long i = 0; i < Nt; i++) {
        fscanf(ft, "%lf %lf %lf", &t_days[i], &flux[i], &err[i]);
    }
    fclose(ft);

    double *model = malloc(sizeof(double) * Nt);
    Calculate_Lightcurve(t_days, Nt, pars, model);

    FILE *fo = fopen("model_folded.csv", "w");
    fprintf(fo, "time_days,flux_model\n");
    for (long i = 0; i < Nt; i++) {
        fprintf(fo, "%.10f,%.10f,%.10f\n", t_days[i], model[i], flux[i]);
    }
    fclose(fo);

    free(t_days);
    free(flux);
    free(err);
    free(model);

    return 0;
}