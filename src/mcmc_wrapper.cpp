#include "mcmc_wrapper.h"

#include <fstream>
#include <iosfwd>
#include <iomanip>
#include <iostream>

#include "python_interrupt_handling.h"
#include "random_generator.h"
#include "configuration.h"


void Run_MCMC(const int tic_id, const int sector)
{
    check_for_and_handle_python_interrupt();
    initialize_configuration(tic_id, sector);
    omp_set_num_threads(get_configuration().get_number_of_threads());
    // Load the MCMC data
    long int buffer_size;
    int py_niter, py_nchains, py_npars, py_nsectors, py_npast;
    double py_dtemp;

    Load_MCMC_Constants(&py_niter, &py_nchains, &py_npars, &py_nsectors,
                        &py_npast, &py_dtemp, &buffer_size);

    const int NCHAINS = py_nchains;
    const int NITER = py_niter;
    const int NPARS = py_npars;
    const int NSECTORS = py_nsectors;
    const int NPAST = py_npast;
    const double dtemp = py_dtemp;

    bounds* limits = new bounds[NPARS];
    bounds* limited = new bounds[NPARS];
    gauss_bounds* gauss_pars = new gauss_bounds[NPARS];
    double* xmap = new double[NPARS];
    double* sigma = new double[NPARS];

    Load_MCMC_Parameter_Info(NPARS, &buffer_size, limits, limited, gauss_pars, xmap,
                             sigma);


    long int* points_per_sector = new long int[NSECTORS];
    long int py_npoints;
    Load_MCMC_Sector_Points(NSECTORS, &buffer_size, points_per_sector, &py_npoints);

    const int NPOINTS = py_npoints;
    double* times = new double[NPOINTS];
    double* fluxes = new double[NPOINTS];
    double* errors = new double[NPOINTS];
    double magdata[5], magerr[4];

    Load_MCMC_Data_Arrays(NPOINTS, &buffer_size, times, fluxes, errors, magdata,
                          magerr);


    const double GAMMA = 2.388 / sqrt(2. * (double)NPARS);

    double best_recorded_log_likelihood = std::numeric_limits<double>::lowest();
    double* logLx = new double[NCHAINS];
    double* logPx = new double[NCHAINS];
    double* temp = new double[NCHAINS];

    double logLy;
    double logPy;
    double logLmap;
    double** x;
    double*** history;

    int* DEtrial_arr = new int[NCHAINS];
    int* acc_arr = new int[NCHAINS];
    int* DEacc_arr = new int[NCHAINS];

    int acc = 0;
    int atrial = 0;
    int DEtrial = 0;
    int DEacc = 0;
    int* index = new int[NCHAINS];

    auto** random_generators_for_chains = new RandomGenerator*[NCHAINS];

    RandomGenerator* random_generator = create_random_generator(0);

    check_for_and_handle_python_interrupt();

    for (int i = 0; i < NCHAINS; i++)
    {
        random_generators_for_chains[i] = create_random_generator(i);
        DEtrial_arr[i] = 0;
        DEacc_arr[i] = 0;
        acc_arr[i] = 0;
    }

    // Allocate memory for the mcmc arrays
    x = (double**)malloc(NCHAINS * sizeof(double));
    for (int i = 0; i < NCHAINS; i++)
    {
        x[i] = (double*)malloc(NPARS * sizeof(double));
    }

    history = (double***)malloc(NCHAINS * sizeof(double));
    for (int i = 0; i < NCHAINS; i++)
    {
        history[i] = (double**)malloc(NPAST * sizeof(double));
        for (int j = 0; j < NPAST; j++)
        {
            history[i][j] = (double*)malloc(NPARS * sizeof(double));
        }
    }


    const double log_LC_PERIOD = xmap[2];
    printf("Loaded period is %f\n", log_LC_PERIOD);
    for (int i = 0; i < NPARS; i++)
    {
        for (int j = 0; j < NCHAINS; j++)
        {
            double tmp1 = get_uniform_random_value(random_generators_for_chains[j]);
            x[j][i] = (limits[i].lo + tmp1 * (limits[i].hi - limits[i].lo));
            if (i == 2)
            {
                x[j][i] = log_LC_PERIOD;
            }
            // 0th chain gets the input parameters
            if (j == 0)
            {
                x[j][i] = xmap[i];
            }
        }
    }

    // Initialize parallel tempering
    temp[0] = 1.0;
    index[0] = 0;

    for (int i = 1; i < NCHAINS; i++)
    {
        temp[i] = temp[i - 1] * dtemp;
        index[i] = i;
    }

    logLmap = Log_Likelihood(times, fluxes, errors, points_per_sector,
                             NSECTORS, xmap, magdata, magerr);
    printf("Log likelihood is %f\n", logLmap);

    Read_Parameters(x, NPARS, NCHAINS);

    // Main MCMC loop starts here
    for (int iter = 0; iter < NITER; iter++)
    {
        check_for_and_handle_python_interrupt();
        int k = iter - (iter / NPAST) * NPAST;

        #if ENABLE_OPENMP
        #pragma omp parallel for schedule(static)
        #endif
        for (int j = 0; j < NCHAINS; j++)
        {
            // Test parameters
            double* y = new double[NPARS];

            // Random number
            double alpha = get_uniform_random_value(random_generators_for_chains[j]);

            // Jump scale and steps
            double jscale = pow(10., -6. + 6. * alpha);

            double* dx = new double[NPARS];
            double dx_mag = 0;
            int jump = 0;
            int chain_id = index[j];
            int jump_type = 0;

            // Take steps in parameter space
            if ((get_uniform_random_value(random_generators_for_chains[j]) < 0.5) && (iter > NPAST))
            {
                jump = 1;
            }

            //gaussian jumps
            if (jump == 0)
            {
                Gaussian_Proposal(x[chain_id], sigma, jscale, temp[j], y, NPARS,
                                  random_generators_for_chains, j);
                jump_type = 1;
            }

            //jump along correlations derived from chain history
            if (jump == 1)
            {
                /* DE proposal; happens after 500 cycles */
                if (chain_id == 0)
                {
                    DEtrial_arr[j]++;
                }

                Differential_Evolution_Proposal(x[chain_id], history[j], y,
                                                NPARS, NPAST, GAMMA, random_generators_for_chains, j);
                jump_type = 2;

                for (int i = 0; i < NPARS; i++)
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
            for (int i = 0; i < NPARS; i++)
            {
                // Reflecting boundary conditions
                while (((limited[i].lo == 1) && (y[i] < limits[i].lo)) ||
                    ((limited[i].hi == 1) && (y[i] > limits[i].hi)))
                {
                    if (y[i] < limits[i].lo)
                    {
                        y[i] = 2.0 * limits[i].lo - y[i];
                    }

                    else
                    {
                        y[i] = 2.0 * limits[i].hi - y[i];
                    }
                }

                // Periodic boundary conditions
                while ((limited[i].lo == 2) && (y[i] < limits[i].lo))
                {
                    y[i] = limits[i].hi + (y[i] - limits[i].lo);
                }

                while ((limited[i].hi == 2) && (y[i] > limits[i].hi))
                {
                    y[i] = limits[i].lo + (y[i] - limits[i].hi);
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
                                             NSECTORS, x[chain_id], magdata, magerr);
            logLy = Log_Likelihood(times, fluxes, errors, points_per_sector,
                                   NSECTORS, y, magdata, magerr);

            /* evaluate new solution */
            alpha = get_uniform_random_value(random_generators_for_chains[j]);

            //Hasting's ratio
            double H = exp((logLy - logLx[chain_id]) / temp[j] + (logPy - logPx[chain_id]));

            //conditional acceptance of y
            if (alpha <= H)
            {
                if (chain_id == 0)
                {
                    acc_arr[j]++;
                }

                for (int i = 0; i < NPARS; i++)
                {
                    x[chain_id][i] = y[i];
                }

                logLx[chain_id] = logLy;

                if ((jump == 1) && (chain_id == 0))
                {
                    DEacc_arr[j]++;
                }
            }


            for (int i = 0; i < NPARS; i++)
            {
                history[j][k][i] = x[chain_id][i];
            }

            free(y);
            free(dx);
            atrial++;
        }

        /********Chain Loop ends**********/

        for (int i = 0; i < NCHAINS; i++)
        {
            acc += acc_arr[i];
            DEacc += DEacc_arr[i];
            DEtrial += DEtrial_arr[i];
            acc_arr[i] = 0;

            /* parallel tempering */
            Ptmcmc(index, temp, logLx, logPx, NCHAINS, random_generator);
        }
        //update progress to screen and write data
        if (iter % 100 == 0)
        {
            //update best parameters
            if (logLx[index[0]] > logLmap)
            {
                for (int i = 0; i < NPARS; i++)
                {
                    xmap[i] = x[index[0]][i];
                }
                logLmap = logLx[index[0]];
            }

            printf("%ld/%ld logL=%.10g acc=%.3g DEacc=%.3g", iter, NITER, logLx[index[0]],
                   (double)(acc) / ((double)atrial),
                   (double)DEacc / (double)DEtrial);
            printf("\n");

            //printf("Parameter values: \n");
            // Print first few parameters
            //for (int i=0; i<5; i++)
            //{
            //  printf("%lf\t", x[index[0]][i]);
            //}
            //printf("\n");
        }

        if (iter % 100 == 0)
        {
            Log_Data(iter, x, logLx, index,
                     points_per_sector, times, fluxes,
                     errors, NPARS, NSECTORS, NCHAINS, best_recorded_log_likelihood);
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
    const int tic = atoi(argv[1]); //461541766;
    const int sector = atoi(argv[2]); //-1;

    Run_MCMC(tic, sector);
    return 0;
}


void Gaussian_Proposal(double* x, double* sigma, double scale, double temp, double* y, const int NPARS,
                       RandomGenerator** random_generators_for_chains, int chain_number)
{
    int n;
    double gamma;
    double sqtemp;
    double* dx = new double[NPARS];

    //scale jumps by temperature
    sqtemp = sqrt(temp);

    //compute size of jumps
    for (n = 0; n < NPARS; n++)
        dx[n] = get_normal_random_value(random_generators_for_chains[chain_number]
        ) * sigma[n] * sqtemp * scale;

    //jump in parameter directions scaled by dx
    for (n = 0; n < NPARS; n++)
    {
        y[n] = x[n] + dx[n];
    }

    free(dx);
    return;
}


void Differential_Evolution_Proposal(double* x, double** history, double* y, const int NPARS, const int NPAST,
                                     const double GAMMA, RandomGenerator** random_generators_for_chains,
                                     int chain_number)
{
    int n;
    int a;
    int b;
    double c = get_normal_random_value(random_generators_for_chains[chain_number]);
    double* dx = new double[NPARS];
    double* epsilon = new double[NPARS];

    //choose two samples from chain history
    a = get_uniform_random_value(random_generators_for_chains[chain_number]) * NPAST;
    b = a;
    while (b == a)
    {
        b = get_uniform_random_value(random_generators_for_chains[chain_number]) * NPAST;
    }

    //compute vector connecting two samples
    for (n = 0; n < NPARS; n++)
    {
        dx[n] = history[b][n] - history[a][n];
        epsilon[n] = dx[n] * (Gaussian(c, 0, 1.e-4) - 0.5);
    }
    //Blocks?

    //90% of jumps use Gaussian distribution for jump size
    if (get_uniform_random_value(random_generators_for_chains[chain_number]) < 0.9)
    {
        for (n = 0; n < NPARS; n++)
        {
            dx[n] *= get_normal_random_value(random_generators_for_chains[chain_number]) * GAMMA;
        }
    }

    //jump along vector w/ appropriate scaling
    for (n = 0; n < NPARS; n++)
    {
        dx[n] += epsilon[n];
        y[n] = x[n] + dx[n];
    }

    free(dx);
    free(epsilon);
    return;
}


/* Other functions*/


void Ptmcmc(int* index, double temp[], double logL[], double logP[], const int NCHAINS,
            RandomGenerator* random_generator)
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
    b = (int)(get_uniform_random_value(random_generator) * (double)(NCHAINS - 1));
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
    H = (heat2 - heat1) / (heat2 * heat1);
    alpha = exp(dlogL * H);
    beta = get_uniform_random_value(random_generator);
    if (alpha >= beta)
    {
        index[a] = oldb;
        index[b] = olda;
    }
    return;
}

// Log the lightcurve file and the data
void Log_Data(int iter, double** x, double* logLx, int* index,
              long int* points_per_sector, double all_sector_phases[], double all_sector_fluxes[],
              double all_sector_uncertainties[], const int NPARS, const int NSECTORS, const int NCHAINS, double&
              best_recorded_log_likelihood)
{
    std::filesystem::path states_path = get_configuration().get_states_path();
    std::ofstream states_file(states_path, std::ios::app);
    states_file << std::scientific << std::setprecision(12);

    // Append to the states file.
    states_file << iter << " " << logLx[index[0]] << " ";
    for (int i = 0; i < NPARS; i++)
    {
        // write inc not cos(inc)
        if ((i == 4) & (get_configuration().should_use_secular_drift() != 1))
        {
            states_file << acos(x[index[0]][i]) << " ";
        }
        else
        {
            states_file << x[index[0]][i] << " ";
        }
    }
    states_file << "\n";

    if (best_recorded_log_likelihood < logLx[index[0]])
    {
        best_recorded_log_likelihood = logLx[index[0]];

        // Overwrite the parameters file.
        std::filesystem::path parameters_path = get_configuration().get_parameters_path();
        std::ofstream parameters_file(parameters_path);
        parameters_file << std::scientific << std::setprecision(12);
        states_file << iter << " " << logLx[index[0]] << " ";
        for (int i = 0; i < NPARS; i++)
        {
            // write inc not cos(inc)
            if ((i == 4) & (get_configuration().should_use_secular_drift() != 1))
            {
                parameters_file << acos(x[index[0]][i]) << " ";
            }
            else
            {
                parameters_file << x[index[0]][i] << " ";
            }
        }
        states_file << "\n";

        // Overwrite the light curve file.
        std::filesystem::path light_curves_path = get_configuration().get_folded_observed_and_model_light_curves_path();
        std::ofstream light_curves_file(light_curves_path);
        light_curves_file << std::scientific << std::setprecision(5);

        const int npars_sector = npars_common + npars_unique;
        double* sector_params = new double[npars_sector];
        int skip_samples = 0;

        for (int i = 0; i < npars_common; i++)
        {
            sector_params[i] = x[index[0]][i];
        }

        for (int sector_number = 0; sector_number < NSECTORS; sector_number++)
        {
            // Now assign the unique parameters
            for (int i = 0; i < npars_unique; i++)
            {
                sector_params[npars_common + i] = x[index[0]][npars_common +
                    npars_unique * sector_number +
                    i];
            }

            if (get_configuration().should_use_secular_drift() == 1)
            {
                double* __temp = new double[npars_sector];
                // Current order: logM1, logM2, logP, sigma_r1, sigma_r2, mu_1, tau_1, mu_2, tau_2, alpha_ref_1, alpha_ref_2
                //                alpha_t1, alpha_t2, (e, i, omega, t0, blending, flux_tune, noise_resc)_j
                __temp[0] = sector_params[0];
                __temp[1] = sector_params[1];
                __temp[2] = sector_params[2];
                __temp[3] = sector_params[15];
                __temp[4] = sector_params[16];
                __temp[5] = sector_params[17];
                __temp[6] = sector_params[18];
                for (int ii = 7; ii <= 18; ii++)
                {
                    __temp[ii] = sector_params[ii - 4];
                }
                // And now move back to __temp
                for (int ii = 0; ii < npars_sector; ii++)
                {
                    sector_params[ii] = __temp[ii];
                }
                delete[] __temp;
            }

            const int Npoints_in_sector = points_per_sector[sector_number];
            double* sector_flux = new double[Npoints_in_sector];
            double* sector_phase = new double[Npoints_in_sector];
            double* sector_uncetainties = new double[Npoints_in_sector];
            double* sector_template = new double[Npoints_in_sector];

            for (int idx = 0; idx < Npoints_in_sector; idx++)
            {
                sector_flux[idx] = all_sector_fluxes[skip_samples + idx];
                sector_phase[idx] = all_sector_phases[skip_samples + idx];
                sector_uncetainties[idx] = all_sector_uncertainties[skip_samples + idx];
            }

            Calculate_Lightcurve(sector_phase, Npoints_in_sector, sector_params,
                                 sector_template);
            for (int i = 0; i < Npoints_in_sector; i++)
            {
                light_curves_file << sector_phase[i] << " " << sector_flux[i] << " "
                                  << sector_template[i] << " " << sector_uncetainties[i] << "\n";
            }

            delete[] sector_flux;
            delete[] sector_phase;
            delete[] sector_uncetainties;
            delete[] sector_template;

            skip_samples += Npoints_in_sector;
        }

        delete[] sector_params;
    }
}

// Read parameters from the existing paramters file
void Read_Parameters(double** X, const int NPARS, const int NCHAINS)
{
    std::filesystem::path path = get_configuration().get_parameters_path();
    if (!std::filesystem::is_empty(path))
    {
        std::cout << "Reading parameters from existing parameters file: " << path << std::endl;
        std::ifstream par_file(path);

        for (int i = 0; i < NPARS; i++)
        {
            double temp;
            par_file >> temp;
            for (int j = 0; j < NCHAINS; j++)
            {
                if (j == 4)
                {
                    temp = cos(temp);
                }
                X[j][i] = temp;
            }
            std::cout << "Read parameters: " << X[0][i] << "\t";
        }
    }
    else
    {
        std::cout << "Parameters file is empty; initializing random parameters." << std::endl;
    }
}
