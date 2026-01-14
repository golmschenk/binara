# How to run binara
## File structure
### Light curve files
There are two folders:
- #### folded_lightcurves
    - The folder that stores the folded light curves before the mcmc runs on them.
Once you have a light curve from the neural network, you find the period, fold it
based on that period, and then save that information to a .txt file. The file will 
contain a header with the number of data points in the .txt file and the period of 
the light curve in days. Then, the file will have 3 columns that contain the
**phase-folded time, flux, and flux error.**

  - Saving the folded light curve:

    ```python
    t, y, yerr, *_ = load_tic_data(TIC_ID = [tic_id], sector = [sector])

    period = [period]

    phase = (t % period) / period

    sorted_idx = np.argsort(phase)
    phase = phase[sorted_idx]
    y = y[sorted_idx]
    yerr = yerr[sorted_idx]

    out_file = "[tic_id]_sector_[sector].txt"
    np.savetxt(out_file, np.column_stack([phase, y, yerr]),
               fmt="%.6f\t%.6f\t%.6f")

    print(f"Saved folded LC to {out_file}")
    ```
- #### mcmc_lightcurves
    - Once we start fitting the parameters, the best fit light curve that the mcmc
produces is written in this directory. In the template there are two types of files, one with
and one without gmag. It is better to use the gmag files because the mcmc uses gmag as an 
additional constraint to get better estimates. The gmag tag means that the g magnitude
was used in estimating the values, but it is just a single value from the magnitudes folder.
There are four columns in these files: **time in days, flux of original folded
light curve, flux of mcmc light curve, and the error in flux of the original light curve.**

#### magnitudes
Example of file:

<table>
  <tr><td>1000.</td><td></td></tr>
  <tr><td>10.1232</td><td>2.76557439849598</td></tr>
  <tr><td>-0.439</td><td>1.1237103791937538</td></tr>
  <tr><td>0.8318800000000002</td><td>0.8063815510019341</td></tr>
  <tr><td>0.6298999999999993</td><td>0.6608454503128677</td></tr>
</table>

The bottom three rows represent colors, but they are not currently being used. We only
use the top 2 rows. The first row is the distance to the star in parsec, which comes from the tic_id.
The first number in the second row is the gmag, and the second number in the second row is
the uncertainty. You can still run the binara code without providing a magnitude file. If you don't
provide one, the code will still run, but it runs with default values. So, for a more accurate
fit, it is best to include a magnitude file as an input.

#### pars
This file contains the same information that is in the chains directory. The file in the chains 
directory contains 2 additional values in the beginning: the log likelihood and the iteration number.
The file in the pars directory contains the parameter values output by the best fit found by the run.
It goes into the file in the chains directory, find which line has the lowest log likelihood,
and then it copies the parameters of that line into the pars file. If you are only using data from one 
sector, then there are 22 parameters. The 22 parameters written in the file are in this order:
logM1 (mass of primary star 1), logM2 (mass of primary star 2), logP (period), e (eccentricity), 
inc (inclination in rad), omega0 (angle or periastron in rad), T0 (initial time), rr1 (radius scaling factor of primary star),
rr2 (radius scaling factor of secondary star), mu_1, tau_1, mu_2, tau_2, alpha_ref_1, alpha_ref_2,
extra_alpha_beam_1, extra_alpha_beam_2, alpha_T1, alpha_T2, blending, flux_tune, and noise_rescaling.
However, if you are using multiple sectors, 7 additional parameters are added for each additional sector.
In likelihood.c, in the Log_Likelihood function, the variable npars_unique defines the values that are
unique to each sector, and util.c shows that npars_unique is a constant set to 7. These
7 parameters are e, i, omega, t0, blending, flux_tune, and noise_rescaling.

#### py_initialize
This is a combined file created by the write_mcmc_data function. The header is 6 values that are the
number of iterations, the number of chains run in that iteration, number of 
parameters, number of sectors, NPAST which is the number of past samples to remember per chain so 
the code learns from where it has already been, and dTemp, which is the difference in temperature between the chains.
This is followed by 22 lines that represent the same 22 parameters discussed in the pars
section in the same order. Each line has 7 entries, if it is a varying parameter, the initial value, 
the lower bound, the upper bound, the boundary condition, step size, and then the last two entries are the
means and sigmas of each parameter.

## Running binara
First, run
```python
from init_data import load_tic_data
import numpy as np

t, y, yerr, *_ = load_tic_data(TIC_ID, sector= sector)

csv_file = "TESS_data_for_TIC_ID_sector.txt"
np.savetxt(
    csv_file,
    np.column_stack((t, y, yerr)),
    delimiter=" ",
    header="time_days flux flux_err"
)
```
Take the file it outputs, and put it in the 
directory binara/data/original_folded_lightcurves. At the top of the file, erase the headers,
and put the number of lines in the file as the first number, and then the period of the light curve
as the second number. The format of the file should look like this:

<table>
  <tr><td>3358</td><td>5.3475350000</td></tr>
  <tr><td>2229.11878707</td><td>14310.743164</td><td>7.433327</td></tr>
</table>
and it will continue the data for 3358 lines.
Then, the data must be folded and binned. Run the folding_and_binning.py script found
in the src folder. This will create a folded version of the lightcurve that has 150 bins.
The output for the data points of the binned version will be saved in data/original_folded_lightcurves/TIC_ID_sector_binned150.txt.
After that, upload the file and make sure it is in zaratan. You can now run the job using the
sbatch job_script.sh command.

job_script.sh:
```bash
#!/bin/bash
#SBATCH --job-name="binara job"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=40G
#SBATCH --time=0-03:00:00

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Job shell script started."
srun cmake-build-release-zaratan1/binara_exe TIC_ID sector 1 0 0 0
```
The job usually takes a couple of hours (~7) to complete. When it is done,
the output files will be the chains file, the mcmc_lightcurves file, and the pars file.
To see how good the fit is, you can run finding_model_lc_fluxes.c which will create a model_folded.csv
by calling the function Calculate_Lightcurve with the parameters from the pars file.
Then, run plotting_obs_vs_model_fluxes.py to see the observed vs. the model light curves in the same
plot. If you want to see the mode lightcurve alone, you can run example_lightcurve.py, when you
pass in the parameters and the 150 phases, and run the program.