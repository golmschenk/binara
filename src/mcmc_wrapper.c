#include "mcmc_wrapper.h"
#include "python_interrupt_handling.h"
#include "random_generator.h"



void Run_MCMC(const int tic, const int sector, const int run_id, const int gmag_flag, const int color_flag, 
              const int secular_drift_flag)
{
  check_for_and_handle_python_interrupt();
  omp_set_num_threads(10);
  // Load the MCMC data
  long int buffer_size;
  int py_niter, py_nchains, py_npars, py_nsectors, py_npast;
  double py_dtemp;

  Load_MCMC_Constants(tic, sector, run_id, secular_drift_flag,
                      &py_niter, &py_nchains, &py_npars, 
                      &py_nsectors, &py_npast, &py_dtemp, &buffer_size);

  const int NCHAINS = py_nchains;
  const int NITER = py_niter;
  const int NPARS = py_npars;
  const int NSECTORS = py_nsectors;
  const int NPAST = py_npast;
  const double dtemp = py_dtemp;

  bounds *limits = malloc(NPARS * sizeof(bounds));
  bounds *limited = malloc(NPARS * sizeof(bounds));
  gauss_bounds *gauss_pars = malloc(NPARS * sizeof(gauss_bounds));
  double *xmap = malloc(NPARS * sizeof(double));
  double *sigma = malloc(NPARS * sizeof(double));

  Load_MCMC_Parameter_Info(tic, sector, run_id, secular_drift_flag, NPARS, &buffer_size,
                          limits, limited, gauss_pars, xmap, sigma);


  long int *points_per_sector = malloc(NSECTORS * sizeof(long int));
  long int py_npoints;
  Load_MCMC_Sector_Points(tic, sector, run_id, secular_drift_flag, NSECTORS, &buffer_size,
                          points_per_sector, &py_npoints);

  const int NPOINTS = py_npoints;
  double *times = malloc(NPOINTS * sizeof(double));
  double *fluxes = malloc(NPOINTS * sizeof(double));
  double *errors = malloc(NPOINTS * sizeof(double));
  double magdata[5], magerr[4];

  Load_MCMC_Data_Arrays(tic, sector, run_id, secular_drift_flag, NPOINTS, &buffer_size,
                        times, fluxes, errors, magdata, magerr);



  const double GAMMA = 2.388/sqrt(2.* (double)NPARS);

  double *logLx = malloc(NCHAINS * sizeof(double));
  double *logPx = malloc(NCHAINS * sizeof(double));
  double *temp = malloc(NCHAINS * sizeof(double));

  double logLy;
  double logPy;
  double logLmap;
  double **x;
  double ***history;

  int *DEtrial_arr = malloc(NCHAINS * sizeof(int));
  int *acc_arr = malloc(NCHAINS * sizeof(int));
  int *DEacc_arr = malloc(NCHAINS * sizeof(int));

  int acc = 0;
  int atrial = 0;
  int DEtrial = 0;
  int DEacc = 0;
  int *index = malloc(NCHAINS * sizeof(int));

  RandomGenerator **random_generators_for_chains = malloc(NCHAINS * sizeof(RandomGenerator *));
  const int NTHREADS = (int)(NCHAINS / 2);

  char chainname[512] = "";
  char outname[512] = "";
  char parname[512] = "";

  Make_Files(tic, sector, run_id, gmag_flag, color_flag, secular_drift_flag, chainname, outname, parname);

  RandomGenerator* random_generator = create_random_generator(0);

  check_for_and_handle_python_interrupt();

  for (int i=0; i<NCHAINS; i++)
  {
      random_generators_for_chains[i] = create_random_generator(i);
    DEtrial_arr[i] = 0;
    DEacc_arr[i] = 0;
    acc_arr[i] = 0;
  }

  // Allocate memory for the mcmc arrays
  x = (double **)malloc(NCHAINS*sizeof(double));
  for(int i=0;i<NCHAINS;i++) 
  {
    x[i]=(double *)malloc(NPARS*sizeof(double));
  }
  
  history = (double ***)malloc(NCHAINS*sizeof(double));
  for(int i=0;i<NCHAINS;i++)
  {

    history[i]=(double **)malloc(NPAST*sizeof(double));
    for(int j=0;j<NPAST;j++)
    {
      history[i][j]=(double *)malloc(NPARS*sizeof(double));
    }
  }


  const double log_LC_PERIOD = xmap[2];
  printf("Loaded period is %f\n", log_LC_PERIOD);
  for (int i=0;i<NPARS;i++) 
  {
    for(int j=0; j<NCHAINS; j++) 
    {
      double tmp1 = get_uniform_random_value(random_generators_for_chains[j]);
      x[j][i] = (limits[i].lo + tmp1*(limits[i].hi - limits[i].lo));
      if (i == 2) 
      {
        x[j][i] = log_LC_PERIOD;
      }
      // 0th chain gets the input parameters
      if (j==0)
      {
        x[j][i] = xmap[i];
      }
    }
  }

  // Initialize parallel tempering
  temp[0] = 1.0;
  index[0] = 0;
  
  for(int i=1; i<NCHAINS; i++) 
  {
    temp[i]  = temp[i-1]*dtemp;
    index[i] = i;
  }

  logLmap = Log_Likelihood(times, fluxes, errors, points_per_sector,
                            NSECTORS, xmap, magdata, magerr, gmag_flag,
                            color_flag, secular_drift_flag);
  printf("Log likelihood is %f\n", logLmap);

  Read_Parameters(parname, x, NPARS, NCHAINS);

  // Main MCMC loop starts here
  for (int iter=0; iter<NITER; iter++) 
  {
    check_for_and_handle_python_interrupt();
    int k = iter - (iter / NPAST) * NPAST;

    #pragma omp parallel for schedule(static) if(ENABLE_OPENMP)
    for(int j=0; j<NCHAINS; j++)
    {
      // Test parameters
      double *y = malloc(NPARS * sizeof(double));

      // Random number
      double alpha = get_uniform_random_value(random_generators_for_chains[j]);

      // Jump scale and steps
      double jscale = pow(10.,-6.+6.*alpha);

      double *dx = malloc(NPARS * sizeof(double));
      double dx_mag = 0;
      int jump = 0;
      int chain_id = index[j];
      int jump_type = 0;

      // Take steps in parameter space
      if((get_uniform_random_value(random_generators_for_chains[j]) < 0.5) && (iter > NPAST))
      {
        jump = 1;
      }

      //gaussian jumps
      if(jump == 0)
      {
          Gaussian_Proposal(x[chain_id], sigma, jscale, temp[j], y, NPARS,
                            random_generators_for_chains, j);
        jump_type = 1;
      }

      //jump along correlations derived from chain history
      if(jump == 1)
      {
        /* DE proposal; happens after 500 cycles */
        if (chain_id == 0)
        {
          DEtrial_arr[j]++;
        }

          Differential_Evolution_Proposal(x[chain_id], history[j], y,
                                          NPARS, NPAST, GAMMA, random_generators_for_chains, j);
        jump_type = 2;

        for (int i=0;i<NPARS;i++) 
        {
          dx_mag += (x[chain_id][i] - y[i]) * (x[chain_id][i] - y[i]);
        }

        if (dx_mag < 1e-6)
        {
            Gaussian_Proposal(x[chain_id], sigma, jscale, temp[j], y, NPARS,
                              random_generators_for_chains, j);
          jump_type = 1;
        }
      }

      // Enforce priors
      for (int i=0;i<NPARS;i++) 
      {

        // Reflecting boundary conditions
        while (((limited[i].lo == 1) && (y[i] < limits[i].lo)) || 
              ((limited[i].hi == 1) && (y[i] > limits[i].hi)))
        {
          if (y[i] < limits[i].lo)
          {
            y[i] = 2.0*limits[i].lo - y[i];
          }

          else
          {
            y[i] = 2.0*limits[i].hi - y[i]; 
          }
        }

        // Periodic boundary conditions
        while ((limited[i].lo == 2) && (y[i] < limits[i].lo))
      {
        y[i] = limits[i].hi + (y[i]-limits[i].lo);
      }

      while ((limited[i].hi == 2) && (y[i] > limits[i].hi))
      {
        y[i] = limits[i].lo + (y[i]-limits[i].hi);
      }

      // Make sure that the mass of the first star is larger than the mass of the second star
      if (y[0] < y[1])
      {
        Swap(&y[0], &y[1]);
      }
    }

    // Fix the period
    y[2] = log_LC_PERIOD;

    // Gaussian priors
    logPx[chain_id] = Log_Prior(NPARS, x[chain_id], gauss_pars);
    logPy = Log_Prior(NPARS, y, gauss_pars);

    //compute current and trial likelihood
    logLx[chain_id] = Log_Likelihood(times, fluxes, errors, points_per_sector,
                            NSECTORS, x[chain_id], magdata, magerr, gmag_flag,
                            color_flag, secular_drift_flag);
    logLy           = Log_Likelihood(times, fluxes, errors, points_per_sector,
                            NSECTORS, y, magdata, magerr, gmag_flag, color_flag, secular_drift_flag);

    /* evaluate new solution */
    alpha = get_uniform_random_value(random_generators_for_chains[j]);

    //Hasting's ratio
    double H = exp((logLy-logLx[chain_id])/temp[j] + (logPy-logPx[chain_id]));

    //conditional acceptance of y
    if (alpha <= H) 
    {

      if (chain_id == 0) 
      {
        acc_arr[j]++;
      }

      for (int i=0;i<NPARS;i++) 
      {
        x[chain_id][i] = y[i];
      }

      logLx[chain_id] = logLy;

      if((jump == 1) && (chain_id == 0)) 
      {
        DEacc_arr[j]++;
      }

    }


  for(int i=0; i<NPARS; i++)
  {
    history[j][k][i] = x[chain_id][i];
  }

  free(y);
  free(dx);
  atrial++;
  }

    /********Chain Loop ends**********/

    for (int i=0; i<NCHAINS; i++)
    {
      acc += acc_arr[i];
      DEacc += DEacc_arr[i];
      DEtrial += DEtrial_arr[i];
      acc_arr[i] = 0;

      /* parallel tempering */
      Ptmcmc(index, temp, logLx, logPx, NCHAINS, random_generator);
    }
    //update progress to screen and write data
    if(iter%100==0) 
    {
      //update best parameters
      if (logLx[index[0]] > logLmap)
      {
        for (int i=0;i<NPARS;i++) 
        {
	        xmap[i]=x[index[0]][i];
        }
        logLmap = logLx[index[0]];
      }

	    printf("%ld/%ld logL=%.10g acc=%.3g DEacc=%.3g",iter,NITER,logLx[index[0]],
	    (double)(acc)/((double)atrial),
	    (double)DEacc/(double)DEtrial);
	    printf("\n");

      //printf("Parameter values: \n");
      // Print first few parameters
      //for (int i=0; i<5; i++)
      //{
      //  printf("%lf\t", x[index[0]][i]);
      //}
      //printf("\n");
    }

    if (iter%100==0)
    {
      Log_Data(chainname, outname, parname, iter, x, logLx, index,
              points_per_sector, times, fluxes, 
              errors, NPARS, NSECTORS, NCHAINS, secular_drift_flag);
    }
  }

  Free_2d(x, NCHAINS);
  Free_3d(history, NCHAINS, NPAST);

  free(limits);
  free(limited);
  free(gauss_pars);
  free(xmap);
  free(sigma);
  free(points_per_sector);
  free(times);
  free(fluxes);
  free(errors);
  free(logLx);
  free(logPx);
  free(temp);
  free(DEtrial_arr);
  free(acc_arr);
  free(DEacc_arr);
  free(index);

  for (int i = 0; i < NCHAINS; i++)
  {
    destroy_random_generator(random_generators_for_chains[i]);
  }
  free(random_generators_for_chains);
  destroy_random_generator(random_generator);
}



int main(int argc, char* argv[])
{
  const int tic = atoi(argv[1]);//461541766;
  const int sector = atoi(argv[2]);//-1;
  const int run_id = atoi(argv[3]);//1;
  const int gmag_flag = atoi(argv[4]);
  const int color_flag = atoi(argv[5]);
  const int secular_drift_flag = atoi(argv[6]);

  Run_MCMC(tic, sector, run_id, gmag_flag, color_flag, secular_drift_flag);
  return 0;
}


void Gaussian_Proposal(double *x, double *sigma, double scale, double temp, double *y, const int NPARS,
                       RandomGenerator **random_generators_for_chains, int chain_number)
{
  int n;
  double gamma;
  double sqtemp;
  double *dx = malloc(NPARS * sizeof(double));

  //scale jumps by temperature
  sqtemp = sqrt(temp);

  //compute size of jumps
  for(n=0; n<NPARS; n++) dx[n] = get_normal_random_value(random_generators_for_chains[chain_number]
          )*sigma[n]*sqtemp*scale;

  //jump in parameter directions scaled by dx
  for(n=0; n<NPARS; n++)
  {
    y[n] = x[n] + dx[n];
  }

  free(dx);
  return;
}


void Differential_Evolution_Proposal(double *x, double **history, double *y, const int NPARS, const int NPAST,
                                     const double GAMMA, RandomGenerator **random_generators_for_chains,
                                     int chain_number)
{
  int n;
  int a;
  int b;
  double c = get_normal_random_value(random_generators_for_chains[chain_number]);
  double *dx = malloc(NPARS * sizeof(double));
  double *epsilon = malloc(NPARS * sizeof(double));

  //choose two samples from chain history
  a = get_uniform_random_value(random_generators_for_chains[chain_number]) * NPAST;
  b = a;
  while(b==a)
  {
    b = get_uniform_random_value(random_generators_for_chains[chain_number]) * NPAST;
  }

  //compute vector connecting two samples
  for(n=0; n<NPARS; n++)
  {
    dx[n] = history[b][n] - history[a][n];
    epsilon[n] = dx[n] * (Gaussian(c, 0, 1.e-4) - 0.5);
  }
  //Blocks?

  //90% of jumps use Gaussian distribution for jump size
  if(get_uniform_random_value(random_generators_for_chains[chain_number]) < 0.9)
  {
    for(n=0; n<NPARS; n++)
    {
      dx[n] *= get_normal_random_value(random_generators_for_chains[chain_number]) * GAMMA;
    }
  }

  //jump along vector w/ appropriate scaling
  for(n=0; n<NPARS; n++)
  {
    dx[n] += epsilon[n];
    y[n] = x[n] + dx[n];
  }

  free(dx);
  free(epsilon);
  return;
}


/* Other functions*/


void Ptmcmc(int *index, double temp[], double logL[], double logP[], const int NCHAINS, RandomGenerator* random_generator)
{
  int a, b;
	int olda, oldb;
	
	double heat1, heat2;
	double logL1, logL2;
  double logP1, logP2;
	double dlogL;
  double dlogP;
  double dlogE;
	double H;
	double alpha;
	double beta;

	/*
	 Possible evidence for over-coupling using this
	 chain swapping scheme?  Neil & Laura & I see that
	 var(logL) < D/2 w/ PT, var(logL) ~ D/2 w/out.  
	 Neil just randomly chooses a pair of chains to 
	 exchange instead of marching down the line in his
	 EMRI code and doesn't see this problem.  But Joey
	 uses this in the eccentric binary code and also doesn't
	 see the problem.  WTF?  In the meantime, I'll switch to
	 randomly proposing a pair, but this is very puzzling.
	 */ 
	
  /* Siddhant: b can be -1, gives seg fault, put bounds on b*/
	b = (int) (get_uniform_random_value(random_generator) * (double)(NCHAINS - 1));
	a = b + 1;
	
	olda = index[a];
	oldb = index[b];
	heat1 = temp[a];
	heat2 = temp[b];
	logL1 = logL[olda];
	logL2 = logL[oldb];
  //logP1 = logP[olda];
  //logP2 = logP[oldb];
	dlogL = logL2 - logL1;
  //dlogP = logP1 - logP2;
	H  = (heat2 - heat1)/(heat2*heat1);
	alpha = exp(dlogL*H);
	beta  = get_uniform_random_value(random_generator);
	if(alpha >= beta)
	{
		index[a] = oldb;
		index[b] = olda;
	}
  return;
}


void Make_Files(const int tic, const int sector, const int run_id, const int gmag_flag, const int color_flag,
                const int secular_drift_flag, char *chainname, char *outname, char *parname)
{

  char prefix[100] = "";
  char suffix[100] = "";
  char tic_num[15] = "";
  char run_num[15] = "";
  char sec_num[15] = "";

  sprintf(prefix, "data");
  sprintf(tic_num, "%d", tic);
  sprintf(run_num, "_%d", run_id);
  if (sector == -1)
  {
    sprintf(sec_num, "_sector_all");
  }
  else
  {
    sprintf(sec_num, "_sector_%d", sector);
  }

  strcat(suffix, tic_num);
  strcat(suffix, sec_num);

  if (gmag_flag)
  {
    strcat(suffix, "_gmag");
  }
  
  if (color_flag)
  {
    strcat(suffix, "_color");
  }

  if (ENABLE_OPENMP)
  {
    strcat(suffix, "_OMP");
  }

  if (secular_drift_flag)
  {
    strcat(suffix, "_drift");
  }

  strcat(suffix, run_num);

  strcat(chainname, prefix);
  strcat(outname, prefix);
  strcat(parname, prefix);
  strcat(chainname, "/chains/chain.");
  strcat(outname,   "/lightcurves/mcmc_lightcurves/");
  strcat(parname, "/pars/pars.");

  strcat(chainname, suffix);
  strcat(outname, suffix);
  strcat(parname, suffix);
  strcat(chainname, ".dat");
  strcat(outname, ".out");
  strcat(parname, ".out");

  printf("Chainfile: %s\n", chainname);
  printf("outfile: %s\n", outname);
  printf("parfile: %s\n", parname);

  FILE *chain_file, *out_file;
  if (!exists(chainname))
  {
    printf("Old chain file not found, creating new file \n");
    chain_file = fopen(chainname, "w");
    fclose(chain_file);
  }

  out_file = fopen(outname, "w");
  fclose(out_file);
  return;
}


// Log the lightcurve file and the data
void Log_Data(char *chainname, char *outname, char *parname, int iter, double **x, double *logLx, int *index,
              long int *points_per_sector, double all_sector_phases[], double all_sector_fluxes[], 
              double all_sector_uncertainties[], const int NPARS, const int NSECTORS, const int NCHAINS, 
              const int secular_drift_flag)
{
  FILE *chain_file, *out_file, *par_file;
  
  if (exists(chainname))
  {
    chain_file = fopen(chainname, "a");
  }
  else
  {
    chain_file = fopen(chainname, "w");
  }

  out_file = fopen(outname, "w");
  par_file = fopen(parname, "w");

  //print parameter chains
  fprintf(chain_file,"%ld %.12g ",iter,logLx[index[0]]);
  for(int i=0; i<NPARS; i++) 
  {
    // write inc not cos(inc)
    if ((i == 4) & (secular_drift_flag != 1))
    {
      fprintf(chain_file,"%.12g ", acos(x[index[0]][i]));
    }
    else
    {
      fprintf(chain_file,"%.12g ",x[index[0]][i]);

    }
    fprintf(par_file, "%.12g ", x[index[0]][i]);

  }

  fprintf(chain_file,"\n");

  // Print the lighcurve for each sector
  const int npars_sector = npars_common + npars_unique;
  double *sector_params = malloc(npars_sector * sizeof(double));
  int skip_samples = 0;

  for (int i=0; i<npars_common; i++)
  {
    sector_params[i] = x[index[0]][i];
  }

  for (int sector_number = 0; sector_number < NSECTORS; sector_number++)
  {
    // Now assign the unique parameters
    for (int i=0; i<npars_unique; i++)
    {
      sector_params[npars_common + i] = x[index[0]][npars_common +
                                            npars_unique * sector_number +
                                            i];
    }

    if (secular_drift_flag == 1)
    {
      double *__temp = malloc(npars_sector * sizeof(double));
      // Current order: logM1, logM2, logP, sigma_r1, sigma_r2, mu_1, tau_1, mu_2, tau_2, alpha_ref_1, alpha_ref_2
      //                alpha_t1, alpha_t2, (e, i, omega, t0, blending, flux_tune, noise_resc)_j
      __temp[0] = sector_params[0];       __temp[1] = sector_params[1];      __temp[2] = sector_params[2];
      __temp[3] = sector_params[15];       __temp[4] = sector_params[16];      __temp[5] = sector_params[17];
      __temp[6] = sector_params[18];       
      for (int ii=7; ii<=18; ii++)
      {
        __temp[ii] = sector_params[ii-4];
      }
      // And now move back to __temp
      for (int ii=0; ii<npars_sector; ii++)
      {
        sector_params[ii] = __temp[ii];
      }
      free(__temp);
    }

    const int Npoints_in_sector = points_per_sector[sector_number];
    double *sector_flux = malloc(Npoints_in_sector * sizeof(double));
    double *sector_phase = malloc(Npoints_in_sector * sizeof(double));
    double *sector_uncetainties = malloc(Npoints_in_sector * sizeof(double));
    double *sector_template = malloc(Npoints_in_sector * sizeof(double));

    for (int index=0; index < Npoints_in_sector; index++)
    {
      sector_flux[index] = all_sector_fluxes[skip_samples + index];
      sector_phase[index] = all_sector_phases[skip_samples + index];
      sector_uncetainties[index] = all_sector_uncertainties[skip_samples + index];
    }

    Calculate_Lightcurve(sector_phase, Npoints_in_sector, sector_params, 
                         sector_template);
    for (int i=0; i<Npoints_in_sector; i++)
    {
      fprintf(out_file,"%12.5e %12.5e %12.5e %12.5e\n",sector_phase[i],sector_flux[i],sector_template[i],
      sector_uncetainties[i]);
    }

    free(sector_flux);
    free(sector_phase);
    free(sector_uncetainties);
    free(sector_template);

    skip_samples += Npoints_in_sector;
  }

  free(sector_params);

  fclose(out_file);
  fclose(chain_file);
  fclose(par_file);
  return;
}

// Read parameters from the existing chain file
void Read_Parameters(char *parname, double **X, const int NPARS, const int NCHAINS)
{
  if (exists(parname))
  {
    printf("Reading pars from existing chain: %s \n", parname);
    FILE *par_file = fopen(parname, "r");
    
    for (int i=0; i<NPARS; i++)
    {
      double temp;
      fscanf(par_file, "%lf\t", &temp);
      for (int j = 0; j<NCHAINS; j++)
      {
        if (j == 4)
        {
          temp = cos(temp);
        }
        X[j][i] = temp;
      }
      printf("Read parameters: %f \t", X[0][i]);
    }
    
  }
  else
  {
    printf("Par file not found; initializing random parameters \n");
  }
  return;
}