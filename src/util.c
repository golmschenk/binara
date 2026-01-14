#include "util.h"

const int npars_common = 15;//19;
const int npars_unique = 7;//3;


double Ran2_Parallel(long *idum, RNG_Vars* state)
{
  int j;
  long k;
  double temp;

  // Update counter
  state-> cts += 1;
  
  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    state->idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) state->iv[j] = *idum;
    }
    state->iy=state->iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=state->idum2/IQ2;
  state->idum2=IA2*(state->idum2-k*IQ2)-k*IR2;
  if (state->idum2 < 0) state->idum2 += IM2;
  j=state->iy/NDIV;
  state->iy=state->iv[j]-state->idum2;
  state->iv[j] = *idum;
  if (state->iy < 1) state->iy += IMM1;

  if ((temp=AM*state->iy) > RNMX) 
  {
    return RNMX;
  }
  else 
  {
    return temp;
  }
}


// gaussian random number
double Gasdev2_Parallel(long *idum, RNG_Vars *state)
{
  double fac, rsq, v1, v2;
  
  if(*idum < 0) state->iset = 0;
  if(state->iset == 0){
    do{
      v1 = 2.0 * Ran2_Parallel(idum, state)-1.0;
      v2 = 2.0 * Ran2_Parallel(idum, state)-1.0;
      rsq = v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    state->gset = v1 * fac;
    state->iset = 1;
    return(v2*fac);
  } else {
    state->iset = 0;
    return (state->gset);
    }
}

/* Free memory from arrays */
void Free_1d(double *arr)
{
  free(arr);
}

void Free_2d(double **arr, int size)
{
  int i;
  for (i=0;i<size;i++){free(arr[i]);}
  free(arr);
}

void Free_3d(double ***arr, int size1, int size2)
{
  int i,j;
  for (i=0;i<size1;i++){
    for (j=0;j<size2;j++){free(arr[i][j]);}
    free(arr[i]);
  }
  free(arr);
}

// Standard gaussian
double Gaussian(double x, double mean, double sigma)
{
  return (1 / sigma / SQRT_2PI) * exp(- pow((x - mean) / sigma, 2.) / 2.);
}


// Function to check if (parameter) file exists
int exists(const char *fname)
{
    FILE *file;
    if (access(fname, R_OK) == 0){
      return 1;
    }
    else {return 0;}
}

void Get_Datafile_Name(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                      char path[])
{
  char fname[256] = "";
  if (secular_drift_flag == 0)
  {
    if (sector == -1)
    {
      sprintf(fname, "%d_sector_all_run_%d.txt", tic, run_id);
    }
    else
    {
      sprintf(fname, "%d_sector_%d_run_%d.txt", tic, sector, run_id);
    }
    strcat(path, fname);
  }
  else
  {
    sprintf(fname, "%d_sector_all_run_%d_drift.txt", tic, run_id);
    strcat(path, fname);
  }
  return;
}


int Load_MCMC_Constants(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                      int *py_niter, int *py_nchains, int *py_npars, 
                      int *py_nsectors, int *py_npast, double *py_dtemp, long int *buffer_size)
{
  char path[1024] = "/home/siddhant/scratch/HB_MCMC/data/py_initialize/";;
  Get_Datafile_Name(tic, sector, run_id, secular_drift_flag, path);
  printf("Reading constants \n");

  if (exists(path) != 1)
  {
    printf("ERROR: Data file does not exist: %s\n", path);
    return 0;
  }
  
  FILE *data_file = fopen(path, "r");
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
                      const int NPARS,  long int *buffer_size,
                      bounds *limits, bounds *limited, gauss_bounds *gauss_pars, 
                      double *X_init, double *sigma)
{
  char path[1024] = "/home/siddhant/scratch/HB_MCMC/data/py_initialize/";;
  Get_Datafile_Name(tic, sector, run_id, secular_drift_flag, path);
  printf("Reading parameter information \n");

  if (exists(path) != 1)
  {
    printf("ERROR: Data file does not exist: %s\n", path);
    return 0;
  }
  
  FILE *data_file = fopen(path, "r");
  fseek(data_file, *buffer_size, SEEK_SET);

  double par_min, par_max, par_mean, par_jump, prior_gauss_mean, prior_gauss_std,
          bc_buffer;
  int boundary_condition;

  for (int ipar=0; ipar<NPARS; ipar++)
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
                      const int NSECTORS, long int *buffer_size, long int *points_per_sector,
                      long int *py_npoints)
{

 char path[1024] = "/home/siddhant/scratch/HB_MCMC/data/py_initialize/";;
  Get_Datafile_Name(tic, sector, run_id, secular_drift_flag, path);
  printf("Reading sector information \n");

  if (exists(path) != 1)
  {
    printf("ERROR: Data file does not exist: %s\n", path);
    return 0;
  }
  
  FILE *data_file = fopen(path, "r");
  fseek(data_file, *buffer_size, SEEK_SET);

  // Read points per sector
  int temp_int;
  *py_npoints=0;

  for (int i=0; i<NSECTORS-1; i++)
  {
    fscanf(data_file, "%d\t", &temp_int);
    points_per_sector[i] = temp_int;
    *py_npoints += temp_int;
  }
  fscanf(data_file, "%d\n", &temp_int);
  points_per_sector[NSECTORS-1] = temp_int;
  *py_npoints += temp_int;

  printf("Data has %ld points \n", *py_npoints);
  *buffer_size = ftell(data_file);
  fclose(data_file);
  return 1;
}

int Load_MCMC_Data_Arrays(const int tic, const int sector, const int run_id, const int secular_drift_flag,
                      const int NPOINTS, long int *buffer_size, double *times, double *fluxes,
                      double *errors, double *magdata, double *magerr)
{
  char path[1024] = "/home/siddhant/scratch/HB_MCMC/data/py_initialize/";;
  Get_Datafile_Name(tic, sector, run_id, secular_drift_flag, path);
  printf("Reading lightcurve and color data \n");

  if (exists(path) != 1)
  {
    printf("ERROR: Data file does not exist: %s\n", path);
    return 0;
  }
  
  FILE *data_file = fopen(path, "r");
  fseek(data_file, *buffer_size, SEEK_SET);

  double time, flux, err;
  for (int i=0; i<NPOINTS; i++)
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


// Taken from https://stackoverflow.com/questions/12733105/c-function-that-counts-lines-in-file
int count_lines(char* fname, long int *line_size)
{
    FILE *file = fopen(fname, "r");
    char buf[BUF_SIZE];
    int counter = 0;
    *line_size = 0;
    for(;;)
    {
        size_t res = fread(buf, 1, BUF_SIZE, file);
        if (ferror(file))
            return -1;

        int i;
        for(i = 0; i < res; i++)
            if (buf[i] == '\n')
                counter++;

        if (feof(file))
            break;
    }

    fclose(file);
    file = fopen(fname, "r");
    for( ;; )
    {
        char c = fgetc( file );
        if( c == EOF || c == '\n' )
            break;
        ++(*line_size);
    }

    fclose(file);
    return counter;
}