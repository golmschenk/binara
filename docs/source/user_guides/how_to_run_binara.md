# How to run

## File structure

The inputs and outputs to binara primarily live in two directories: `input_data` and `sessions`. Each target has its own subdirectory, named by its identifier (e.g., `tic_id_220052771_sector_6`). All files for a target are stored as a flat list directly within that subdirectory -- there are no nested folders.

### Input data directory

The `input_data` directory contains one subdirectory per target. Each target directory contains the following files:

#### `folded_observed_light_curve.txt`

The folded light curve to be used as input for the MCMC run. The file has a single-line header containing two values: the number of data points and the period of the light curve in days. The remaining lines have 3 columns: **phase-folded time, normalized flux, and flux error**.

Saving the folded light curve from a loaded dataset might look like:

```python
t, y, yerr, *_ = load_tic_data(TIC_ID=[tic_id], sector=[sector])

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

#### `magnitudes_and_colors.txt`

Contains distance and photometric information about the target. The file has 5 rows:

| Row | Contents                                           |
|-----|----------------------------------------------------|
| 1   | Distance to the star in parsecs                    |
| 2   | G magnitude and its uncertainty                    |
| 3   | V-B color and its uncertainty (not currently used) |
| 4   | B-G color and its uncertainty (not currently used) |
| 5   | G-T color and its uncertainty (not currently used) |

Only the distance (row 1) and the G magnitude with its uncertainty (row 2) are used by the code. The color rows are read but not currently used for fitting. You can still run binara without providing a magnitude file; the code will use default values, but providing one produces a more accurate fit.

#### `py_initialize.txt`

A combined initialization file created by the `write_mcmc_data` function. The first line is a header containing 6 values:

1. Number of iterations
2. Number of chains
3. Number of parameters
4. Number of sectors
5. NPAST -- the number of past samples to remember per chain (allows the sampler to learn from where it has already been)
6. dTemp -- the temperature spacing between chains

Following the header, the file contains per-parameter rows with initialization values and prior information, then sector point counts, then the light curve data arrays and magnitude data.

### Sessions directory

The `sessions` directory contains the output from MCMC runs. Each run produces a subdirectory. By default, the subdirectory is named the same as the corresponding `input_data` target directory (e.g., `tic_id_220052771_sector_6`). If the `prefix_session_directory_with_datetime` option is set to `true` in the configuration, the directory name is prefixed with a timestamp (e.g., `2026_02_09_12_03_30_tic_id_220052771_sector_6`).

Each session directory contains the following files:

#### `py_initialize.txt`, `magnitudes_and_colors.txt`, `folded_observed_light_curve.txt`

Copies of the corresponding input files, preserved for reproducibility.

#### `configuration.toml`

A copy of the configuration used for the run.

#### `parameters.txt`

Contains the best-fit parameter values found by the MCMC run. These are the parameters from the chain state with the lowest log likelihood. For a single-sector run, there are 22 parameters, written on a single line in the following order:

1. logM1 -- log mass of primary star
2. logM2 -- log mass of secondary star
3. logP -- log period
4. e -- eccentricity
5. inc -- inclination (radians)
6. omega0 -- argument of periastron (radians)
7. T0 -- initial time
8. rr1 -- radius scaling factor of primary star
9. rr2 -- radius scaling factor of secondary star
10. mu\_1, tau\_1, mu\_2, tau\_2 -- limb darkening parameters
11. alpha\_ref\_1, alpha\_ref\_2 -- reference albedos
12. extra\_alpha\_beam\_1, extra\_alpha\_beam\_2 -- beaming albedo corrections
13. alpha\_T1, alpha\_T2 -- temperature-related albedos
14. blending
15. flux\_tune
16. noise\_rescaling

There are 19 parameters common across all sectors and 3 parameters unique to each sector. For multi-sector runs, 3 additional parameters are added for each additional sector. These per-sector parameters are: blending, flux\_tune, and noise\_rescaling.

#### `states.txt`

Contains the full chain states from the MCMC run. Each line includes the the iteration number, log likelihood, and all parameter values for that state.

#### `folded_observed_and_model_light_curves.txt`

The best-fit model light curve produced by the MCMC alongside the original observed data. This file has four columns: **time in days, flux of the original folded light curve, flux of the MCMC model light curve, and the flux error of the original light curve**.
```