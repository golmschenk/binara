"""
Main python script to run the HB MCMC script that calls
the C MCMC function.
"""
import os
from pathlib import Path

import numpy as np


def load_tic_data(tic_id: int, sector: int = -1):
    """
    Parameters:
        tic_id: int
        sector: int (default = -1)
    if sector == -1 then load all the sector data
    into one array
    Returns:
        A list containing
        tdata: np.ndarray
        ydata: np.ndarray
        yerrdata: np.ndarray
        magdata: np.ndarray
        magerr: np.ndarray
        points_per_sector: np.ndarray (int)
        NSECTORS: int
        Period: float
    """
    # Load lightcurve data
    data_directory = Path('data')
    target_directory = data_directory.joinpath(f'tic_id_{tic_id}_sector_{sector}')
    folded_light_curve_path = target_directory.joinpath('folded_observed_light_curve.txt')
    header = np.genfromtxt(folded_light_curve_path, max_rows=1)
    tdata, ydata, yerrdata = np.genfromtxt(folded_light_curve_path, skip_header=1).T
    # process header
    nsectors = 1
    npoints = len(tdata)
    points_per_sector = []
    if sector == -1:
        nsectors = int(header[0])
        for i in range(nsectors):
            points_per_sector.append(int(header[i + 1]))
    else:
        points_per_sector.append(npoints)
    period = float(header[-1])
    # Load magdata
    magdata, magerr = [], []
    magnitudes_path = target_directory.joinpath('magnitudes.txt')
    # Fill empty gmag data
    with magnitudes_path.open() as magnitudes_file:
        distance = magnitudes_file.readline().split()[0]
        magdata.append(float(distance))
        gmag, gmag_err = magnitudes_file.readline().split()
        magdata.append(float(gmag))
        magerr.append(float(gmag_err))
        for i in range(3):
            color, color_err = magnitudes_file.readline().split()
            magdata.append(float(color))
            magerr.append(float(color_err))
    # Prepare final state
    all_data = [tdata, ydata, yerrdata, magdata, magerr, points_per_sector,
                nsectors, period]
    print("Sectors found: ", nsectors)
    return all_data


def set_mcmc_pars(tic_id: int, sector: int = -1, secular_drift_sources=False):
    """
    Python function to set what parameters to use and what parameters to fix
    Parameters:
            tic_id: int
            sector: int (default = -1)
            run_type: str (default = "plain)
            Choose between "plain", "gmag" and "color"
    """
    all_data = load_tic_data(tic_id=tic_id, sector=sector)
    (tdata, ydata, yerrdata, magdata, magerr, points_per_sector,
     nsectors, period) = all_data
    # set the total number of parameters in the mcmc run
    pars_common = ["logM1", "logM2", "P", "e", "cos_inc", "omega0", "T0", "log rr1", "log rr2",
                   "mu 1", "tau 1", "mu 2", "tau 2", "ref 1", "ref 2", "log beam 1", "log beam 2",
                   "log Teff 1", "log Teff 2"]
    pars_unique = ["blending", "flux tune", "noise resc"]
    if secular_drift_sources:
        pars_common = ["logM1", "logM2", "P", "log rr1", "log rr2", "mu 1", "tau 1", "mu 2", "tau 2",
                       "ref 1", "ref 2", "log beam 1", "log beam 2", "log Teff 1", "log Teff 2"]
        pars_unique = ["e", "cos_inc", "omega0", "T0", "blending", "flux tune", "noise resc"]

    npars_common = len(pars_common)
    npars_unique = len(pars_unique)
    npars = npars_common + nsectors * npars_unique
    # The parameter array values are [fix/vary, value, min, max, 
    #                        boundary condition (ref/per), jump size]
    small_num = 1.e-30
    big_num = 1.e+30
    parameters = {
        "logM1": ["vary", 0., -1.5, 2., "ref", 1.e-2],
        "logM2": ["vary", 0., -1.5, 2., "ref", 1.e-2],
        "P": ["fix", period, -2, 3., "ref", small_num],
        "e": ["vary", 0.1, 0., 1., "ref", 1.e-2],
        "cos_inc": ["vary", 0., 0., 1., "ref", 1.e-3],
        "omega0": ["vary", 0., -np.pi, np.pi, "per", 1.e-3],
        "T0": ["vary", 0., 0., 2. * np.pi, "ref", 1.e-3],
        "log rr1": ["vary", 0., -5., 5., "ref", 1.e-1],
        "log rr2": ["vary", 0., -5., 5., "ref", 1.e-1],
        "mu 1": ["vary", 0.1, 0., 0.4, "ref", 1.e-2],
        "tau 1": ["vary", 0.2, 0.15, 0.55, "ref", 1.e-2],
        "mu 2": ["vary", 0.1, 0., 0.4, "ref", 1.e-2],
        "tau 2": ["vary", 0.2, 0.15, 0.55, "ref", 1.e-2],
        "ref 1": ["vary", 0., 0., 2., "ref", 1.e-2],
        "ref 2": ["vary", 0., 0., 2., "ref", 1.e-2],
        "log beam 1": ["vary", 0., -0.3, 0.3, "ref", 1.e-2],
        "log beam 2": ["vary", 0., -0.3, 0.3, "ref", 1.e-2],
        "log Teff 1": ["vary", 0., -5., 5., "ref", 1.e-1],
        "log Teff 2": ["vary", 0., -5., 5., "ref", 1.e-1],
        "blending": ["vary", 0., 0., 1., "ref", 1.e-3],
        "flux tune": ["vary", 1., 0.99, 1.01, "ref", 1.e-5],
        "noise resc": ["vary", 0., -0.11, 0.11, "ref", 1.e-4]
    }
    # Change the jump size for the parameters that you want to fix
    # change "ref" to 1 and "per" to 2
    for key in parameters.keys():
        if parameters[key][0] == "fix":
            values = parameters[key]
            values[-1] = small_num
            parameters[key] = values
        if parameters[key][4] == "ref":
            values = parameters[key]
            values[4] = 1
            parameters[key] = values
        elif parameters[key][4] == "per":
            values = parameters[key]
            values[4] = 2
            parameters[key] = values
    # change period to log period
    values = parameters["P"]
    values[1] = np.log10(values[1])
    parameters["P"] = values
    # Set the gaussian prior information for the parameters (big number if uniform prior)
    # array values are: prior means and sigmas
    GAUSS_PRIORS = {
        "logM1": [0., big_num],
        "logM2": [0., big_num],
        "P": [0., big_num],
        "e": [0., big_num],
        "cos_inc": [0., big_num],
        "omega0": [0., big_num],
        "T0": [0., big_num],
        "log rr1": [0., 1.],
        "log rr2": [0., 1.],
        "mu 1": [0.16, 0.04],
        "tau 1": [0.34, 0.04],
        "mu 2": [0.16, 0.04],
        "tau 2": [0.34, 0.04],
        "ref 1": [1., 0.2],
        "ref 2": [1., 0.2],
        "log beam 1": [0., 0.1],
        "log beam 2": [0., 0.1],
        "log Teff 1": [0., 1.0],
        "log Teff 2": [0., 1.0],
        "blending": [0., big_num],
        "flux tune": [0., big_num],
        "noise resc": [0., 0.027522]
    }
    # Set the number of chains and the temperature difference between the chains
    niter = 1000000
    nchains = 50
    npast = 500
    d_temp = 1.4
    # Package everything up
    constants = [niter, nchains, npars, nsectors, npast, d_temp]
    parameter_info = [parameters, GAUSS_PRIORS, points_per_sector]
    arrays = [tdata, ydata, yerrdata, magdata, magerr]
    misc = [pars_common, pars_unique]
    return [constants, parameter_info, arrays, misc]


def write_mcmc_data(tic_id: int, sector: int = -1, run_id: int = 1, secular_drift_sources=False):
    """
    Write MCMC data
    Same parameter as set_mcmc_pars
    """
    if set_mcmc_pars(tic_id=tic_id, sector=sector) is None:
        return None
    constants, parameter_info, arrays, misc = set_mcmc_pars(
        tic_id=tic_id, sector=sector, secular_drift_sources=secular_drift_sources)
    data_directory = Path('data')
    target_directory = data_directory.joinpath(f'tic_id_{tic_id}_sector_{sector}')
    py_initialize_path = target_directory.joinpath('py_initialize.txt')
    # .... write the data ...
    with py_initialize_path.open('w') as py_initialize_file:
        # write the constants (NITER, NCHAINS, NPARS, NSECTORS, NPAST, dTemp)
        for const in constants[:-1]:
            py_initialize_file.write("%d\t" % const)
        py_initialize_file.write("%.6f\n" % constants[-1])
        # write the parameter information (repeat the last few pars for multiple sectors)
        PARS_common, PARS_unique = misc
        print(PARS_common, PARS_unique)
        PARAMETERS, GAUSS_PRIORS, points_per_sector = parameter_info
        for par in PARS_common:
            vals = PARAMETERS[par]
            for val in vals[1:]:
                py_initialize_file.write("%.6f\t" % float(val))
            mean, median = GAUSS_PRIORS[par]
            py_initialize_file.write("%.6f\t%.6f\n" % (mean, median))
        NSECTORS = constants[-3]
        for i in range(NSECTORS):
            for par in PARS_unique:
                vals = PARAMETERS[par]
                for val in vals[1:]:
                    py_initialize_file.write("%.6f\t" % float(val))
                mean, median = GAUSS_PRIORS[par]
                py_initialize_file.write("%.6f\t%.6f\n" % (mean, median))
        # Now write the points per sector array
        for pps in points_per_sector[:-1]:
            py_initialize_file.write("%d\t" % pps)
        py_initialize_file.write("%d\n" % points_per_sector[-1])
        # Write the array information
        tdata, ydata, yerrdata, magdata, magerr = arrays
        for (t, f, err) in zip(tdata, ydata, yerrdata):
            py_initialize_file.write("%.6f\t%.6f\t%.6f\n" % (t, f, err))
        # Finally the color data        
        for m in magdata[:-1]:
            py_initialize_file.write("%.6f\t" % m)
        py_initialize_file.write("%.6f\n" % magdata[-1])
        for m in magerr[:-1]:
            py_initialize_file.write("%.6f\t" % m)
        py_initialize_file.write("%.6f\n" % magerr[-1])
    print("Successfully written to ", py_initialize_path)
    return 1


if __name__ == "__main__":
    write_mcmc_data(220052771, sector=6, run_id=1, secular_drift_sources=False)
