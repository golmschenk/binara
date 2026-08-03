# How to run

## Basic file structure

The inputs and outputs to binara primarily live in two directories: `input_data` and `sessions`. Each target has its own subdirectory, named by its identifier (e.g., `tic_id_220052771_sector_6`). For a complete description of the files and their content structure, see [the file structure document](file_structure.md). The below is a brief overview along with how the files are produced.

The `input_data` directory contains the data which will be used as input to binara MCMC. Inside, is a separate directory for each target and sector (e.g., `tic_id_220052771_sector_6`). Each of these should start with 2 files:
* `folded_observed_light_curve.txt`
* `magnitudes_and_colors.txt`

Currently, there is no built in mechanism in binara to produce these files, and they must be produced independently.

Before running the MCMC, a third file is generated in the `input_data` directory by running the `from binara.init_data import write_mcmc_data` (calling it from the directory that contains the `input_directory`). This produces the `py_initialize.txt` file, such that a given target directory now contains: 
* `folded_observed_light_curve.txt`
* `magnitudes_and_colors.txt`
* `py_initialize.txt`

Lastly, before running the MCMC, you need a `configuration.toml` file in the directory you are running the MCMC from (the same directory that contains `input_data`). These are the global MCMC settings that will be applied to all MCMC runs.

The MCMC code will produce an output directory in the `sessions` directory. Once run, it will contain: 
* `folded_observed_light_curve.txt`
* `magnitudes_and_colors.txt`
* `py_initialize.txt`
* `configuration.toml`
* `parameters.txt`
* `states.txt`
* `folded_observed_and_model_light_curves.txt`

About half of these are simply copies of the input saved for logging purposes. The rest are the output of the MCMC.
