'''
    Main python script to run the HB MCMC script that calls
    the C MCMC function.
'''
import numpy as np
import glob, sys, os, re

'''
    Parameters:
        TIC: int
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
'''
def load_tic_data(TIC: int, sector:int = -1):
    # Load lightcurve data
    DATA_DIR = "/Users/golmschenk/Code/heart_beat_generation/data/lightcurves/folded_lightcurves/"
    FNAME = DATA_DIR
    if sector == -1: FNAME += "%d_sector_all.txt" % TIC; print("Loading all lightcurves of TIC %d" %TIC)
    else: FNAME += "%d_sector_%d.txt" % (TIC, sector)
    header = np.genfromtxt(FNAME, max_rows=1)
    tdata, ydata, yerrdata = np.genfromtxt(FNAME, skip_header=1).T
    # process header
    NSECTORS = 1
    NPOINTS = len(tdata)
    points_per_sector = []
    if sector == -1: 
        NSECTORS = int(header[0])
        for i in range(NSECTORS): points_per_sector.append(int(header[i+1]))
    else: points_per_sector.append(NPOINTS)
    period = float(header[-1])
    # Load magdata
    magdata, magerr = [], []
    MAG_DIR = "/Users/golmschenk/Code/heart_beat_generation/data/magnitudes/"
    FNAME = MAG_DIR + "%d.txt" % TIC
    # Fill empty gmag data
    
    if os.path.isfile(FNAME):
        with open(FNAME, "r") as f:
            distance = f.readline().split()[0]
            magdata.append(float(distance))
            gmag, gmag_err = f.readline().split()
            magdata.append(float(gmag)); magerr.append(float(gmag_err))
            for i in range(3):
                color, color_err = f.readline().split()
                magdata.append(float(color)); magerr.append(float(color_err))
    else:
        magdata = [100, 10, 0, 0, 0]
        magerr = [10000., 10000., 10000., 10000.]
    # Prepare final state
    all_data = [tdata, ydata, yerrdata, magdata, magerr, points_per_sector, 
                NSECTORS, period]
    print("Sectors found: ", NSECTORS)
    return all_data

'''
    Python function to set what parameters to use and what parameters to fix
    Parameters:
            TIC: int
            sector: int (default = -1)
            run_type: str (default = "plain)
            Choose between "plain", "gmag" and "color"
'''
def set_mcmc_pars(TIC: int, sector: int = -1, secular_drift_sources = False):
    try: all_data = load_tic_data(TIC=TIC, sector=sector)
    except Exception as error:
        return None
    (tdata, ydata, yerrdata, magdata, magerr, points_per_sector, 
    NSECTORS, period) = all_data
    # set the total number of parameters in the mcmc run
    PARS_common = ["logM1", "logM2", "P", "e", "cos_inc", "omega0", "T0", "log rr1", "log rr2",
                "mu 1", "tau 1", "mu 2", "tau 2", "ref 1", "ref 2", "log beam 1", "log beam 2",
                "log Teff 1", "log Teff 2"]
    PARS_unique = ["blending", "flux tune", "noise resc"]
    if secular_drift_sources:
        PARS_common = ["logM1", "logM2", "P", "log rr1", "log rr2", "mu 1", "tau 1", "mu 2", "tau 2", 
                       "ref 1", "ref 2", "log beam 1", "log beam 2", "log Teff 1", "log Teff 2"]
        PARS_unique = ["e", "cos_inc", "omega0", "T0", "blending", "flux tune", "noise resc"]

    NPARS_common = len(PARS_common)
    NPARS_unique = len(PARS_unique)
    NPARS = NPARS_common + NSECTORS * NPARS_unique
    # The parameter array values are [fix/vary, value, min, max, 
    #                        boundary condition (ref/per), jump size]
    SMALL_NUM = 1.e-30
    BIG_NUM = 1.e+30
    PARAMETERS = {
                "logM1": ["vary", 0., -1.5, 2., "ref", 1.e-2],
                "logM2": ["vary", 0., -1.5, 2., "ref", 1.e-2],
                "P": ["fix", period, -2, 3., "ref", SMALL_NUM],
                "e": ["vary", 0.1, 0., 1., "ref", 1.e-2],
                "cos_inc": ["vary", 0., 0., 1., "ref", 1.e-3],
                "omega0": ["vary", 0., -np.pi, np.pi, "per", 1.e-3],
                "T0": ["vary", 0., 0., 2.*np.pi, "ref", 1.e-3],
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
    for key in PARAMETERS.keys():
        if PARAMETERS[key][0] == "fix": 
            values = PARAMETERS[key]
            values[-1] = SMALL_NUM
            PARAMETERS[key] = values
        if PARAMETERS[key][4] == "ref": 
            values = PARAMETERS[key]
            values[4] = 1
            PARAMETERS[key] = values
        elif PARAMETERS[key][4] == "per":
            values = PARAMETERS[key]
            values[4] = 2
            PARAMETERS[key] = values
    # change period to log period
    values = PARAMETERS["P"]
    values[1] = np.log10(values[1])
    PARAMETERS["P"] = values
    # Set the gaussian prior information for the parameters (big number if uniform prior)
    # array values are: prior means and sigmas
    GAUSS_PRIORS = {
                "logM1": [0., BIG_NUM],
                "logM2": [0., BIG_NUM],
                "P": [0., BIG_NUM],
                "e": [0., BIG_NUM],
                "cos_inc": [0., BIG_NUM],
                "omega0": [0., BIG_NUM],
                "T0": [0., BIG_NUM],
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
                "blending": [0., BIG_NUM],
                "flux tune": [0., BIG_NUM],
                "noise resc": [0., 0.027522]
                }
    # Set the number of chains and the temperature difference between the chains
    NITER = 1000000
    NCHAINS = 50
    NPAST = 500
    dTemp = 1.4
    # Package everything up
    constants = [NITER, NCHAINS, NPARS, NSECTORS, NPAST, dTemp]
    parameter_info = [PARAMETERS, GAUSS_PRIORS, points_per_sector]
    arrays = [tdata, ydata, yerrdata, magdata, magerr]
    misc = [PARS_common, PARS_unique]
    return [constants, parameter_info, arrays, misc]

'''
    Write MCMC data
    Same parameter as set_mcmc_pars
'''
def write_mcmc_data(TIC: int, sector: int = -1, run_id: int = 1, secular_drift_sources = False):
    if set_mcmc_pars(TIC=TIC, sector=sector) is None: return None
    constants, parameter_info, arrays, misc = set_mcmc_pars(
                                                TIC=TIC, sector=sector, secular_drift_sources=secular_drift_sources)
    OUTDIR = "/Users/golmschenk/Code/heart_beat_generation/data/py_initialize/"
    if sector == -1: sector = "all"
    FNAME = OUTDIR + "%d_sector_%s_run_%d.txt" % (TIC, str(sector), run_id)
    if secular_drift_sources: FNAME = OUTDIR + "%d_sector_%s_run_%d_drift.txt" % (TIC, str(sector), run_id)
    # .... write the data ...
    with open(FNAME, "w") as file:
        # write the constants (NITER, NCHAINS, NPARS, NSECTORS, NPAST, dTemp)
        for const in constants[:-1]: file.write("%d\t" % const)
        file.write("%.6f\n" % constants[-1])
        # write the parameter information (repeat the last few pars for multiple sectors)
        PARS_common, PARS_unique = misc
        print(PARS_common, PARS_unique)
        PARAMETERS, GAUSS_PRIORS, points_per_sector = parameter_info
        for par in PARS_common:
            vals = PARAMETERS[par]
            for val in vals[1:]: file.write("%.6f\t" % float(val))
            mean, median = GAUSS_PRIORS[par]
            file.write("%.6f\t%.6f\n" % (mean, median))
        NSECTORS = constants[-3]
        for i in range(NSECTORS):
            for par in PARS_unique:
                vals = PARAMETERS[par]
                for val in vals[1:]: file.write("%.6f\t" % float(val))
                mean, median = GAUSS_PRIORS[par]
                file.write("%.6f\t%.6f\n" % (mean, median))
        # Now write the points per sector array
        for pps in points_per_sector[:-1]: file.write("%d\t" % pps)
        file.write("%d\n" % points_per_sector[-1])
        # Write the array information
        tdata, ydata, yerrdata, magdata, magerr = arrays
        for (t, f, err) in zip(tdata, ydata, yerrdata):
            file.write("%.6f\t%.6f\t%.6f\n" % (t, f, err))
        # Finally the color data        
        for m in magdata[:-1]: file.write("%.6f\t" % m)
        file.write("%.6f\n" % magdata[-1])
        for m in magerr[:-1]: file.write("%.6f\t" % m)
        file.write("%.6f\n" % magerr[-1])
    print("Successfully written to ", FNAME)
    return 1

def get_all_tics():
    TICS_LOC = glob.glob("/Users/golmschenk/Code/heart_beat_generation/data/lightcurves/folded_lightcurves/*.txt")
    TICS_LIST = []

    for loc in TICS_LOC:
        pattern = r'\d+'  # Regular expression pattern to match one or more digits
        numbers = re.findall(pattern, loc)
        
        if numbers:
            TICS_LIST.append(int(numbers[0]))

    TICS_LIST = list(set(TICS_LIST))
    return TICS_LIST

def make_par_files(tic: int):
    DATA_DIR = "/Users/golmschenk/Code/heart_beat_generation/data/"
    chainfname = DATA_DIR + "chains/chain.%d_sector_all_gmag_OMP_1.dat" % tic
    data = np.loadtxt(chainfname,  skiprows=5000)
    pars = data[-110, 2:]

    parfname = DATA_DIR + "pars/pars.%d_sector_all_gmag_OMP_1.out" % tic
    with open(parfname, "w") as file:
        for i in range(pars.shape[0]): file.write("%.7f\t" % pars[i])
    return


import subprocess
def move_teo_tics(tic: int, sector: int):
    data_dir = r'/Users/golmschenk/Code/heart_beat_generation/data/chains/'
    data_dir += 'chain.%d_sector_%d_gmag_OMP_1.dat' % (tic, sector)
    subprocess.call(["cp", data_dir, "teo_tics"])
    return

def main():
    TICS = get_all_tics()
    HB_EB_tics = [53194798, 59090149, 75751455, 77372151,
    102289966, 145974291, 146820000, 150284425,
    169449165, 188258837, 189333345, 191457397,
    219224399, 219340003, 220052771, 220399390,
    233841767, 237957506, 240918551, 249033993,
    261724521, 268599292, 269323587, 270859445,
    272822324, 283539216, 283839379, 285128834,
    293525651, 293618358, 296633772, 296633787,
    305252841, 305454334, 306107122, 306507325,
    316119373, 316119373, 317050129, 326374705,
    326484902, 330943024, 332533061, 334074735,
    336538437, 338282749, 342520115, 343173397,
    343626774, 344586348, 356169556, 363679519,
    371966579, 375977934, 376499580, 378275980,
    390337472, 391244527, 396201191, 408618201,
    444555685, 448860246, 451129217, 451708707]

    HB_ECL_tics = [90547242, 118305806, 219224398, 
    219707463, 255876795, 269692669, 277236190, 
    282876586, 284144129, 293618358, 293950421, 
    302828770, 303577635, 316119373, 326484902,
    327725463, 330605074, 336823975, 367944808,
    369033532, 395413286, 427377458, 441626681,
    462940910]

    Vanilla_tics = [169398679, 431305729, 2021686529, 2021686530, 337376644, 371584261, 
    363674500, 461541766, 312344965, 312344969, 137810570, 154225294, 336168208, 365831185, 
    444446481, 470847250, 336091799, 356245399, 374211609, 316668953, 457264281, 275703451, 
    124412957, 251972126, 316621853, 53824793, 275586333, 468721314, 78379043, 421924260, 
    383518759, 78379048, 286324138, 271554988, 79686189, 156175661, 174094640, 362081271, 
    427657141, 427377463, 154348089, 272822330, 274686141, 256308031, 236881602, 358370373, 
    469289800, 415968973, 64290127, 420388688, 45895631, 405319631, 352770261, 462292185, 
    352835929, 286221659, 746425948, 10602978, 418183908, 269511526, 343878759, 2046878952, 
    153043302, 456042860, 390334189, 305448302, 272924271, 252588526, 441512942, 405320687, 
    133664118, 172985206, 271652854, 82091902]
    
    Shporer_kics = [4659476, 5017127, 5090937, 5790807,
                5818706, 5877364, 5960989, 6370558, 
                6775034, 8027591, 8164262, 9016693,
                9965691, 10334122, 11071278, 11403032,
                11649962, 11923629, 12255108]

    TICS = [378275980, 391244527, 444555685, 283839379, 451708707]
    sectors = [[14, 15, 41, 54, 55], [10,11,37,38,63,64], [18,24,25,58], 
               [17,18,24,58], [10,11,37,64]]
    new_drifitng_TICS = [240918551, 305454334, 115394297, 239714064, 283539216]
    new_drifting_sectors = [[17, 18, 58], [10, 11, 63, 64], [19, 43, 44, 45, 59], [19, 43, 44, 45, 59], [16, 17, 56, 57]]
    paper_drifting = [115394297, 178739533, 239714064, 240918551, 283539216, 283839379, 305454334, 336538437, 378275980, 391244527,444555685, 451708707]
    #sectors_tic_336538437 = [7, 33]
    #for sector in sectors_tic_336538437: write_mcmc_data(178739533, sector=sector, run_id=1)
    # for tic in paper_drifting: write_mcmc_data(tic, sector=-1, run_id=10, secular_drift_sources=True)

    write_mcmc_data(220052771, sector=6, run_id=1, secular_drift_sources=False)
    #write_mcmc_data(240918551, sector=-1, run_id=2)
    #write_mcmc_data(240918551, sector=17, run_id=2)
    #write_mcmc_data(240918551, sector=18, run_id=2)
    #write_mcmc_data(240918551, sector=58, run_id=2)
    
    #write_mcmc_data(1111111, 1, run_id=1)
    #write_mcmc_data(2222222, 1, run_id=1)
    #write_mcmc_data(3333333, 1, run_id=1)

    #diff = list((set(TICS) - set(HB_EB_tics)) - set(HB_ECL_tics))
    #print(len(diff), diff)
    #for sector in [16, 17, 18, 24, 57, 58]:
    #    write_mcmc_data(336538437, sector, run_id=1)

    all_tics = []
    #with open("../tic_classification.txt", "r") as file:
    #    for line in file.readlines():
    #        tic = int(line.split()[0])
    #        good_flag = int(line.split()[1])
    #        if good_flag == 1: all_tics.append(tic)

    #for i, tic in enumerate(TICS):
    #    for sector in sectors[i]:
    #        move_teo_tics(tic, sector)
            #write_mcmc_data(tic, sector)
    #for TIC in all_tics:
    #    write_mcmc_data(TIC, sector=-1)
        #make_par_files(TIC)
        #if int(TIC) == 283539216: continue
        #else: 
        #    if(write_mcmc_data(TIC, sector=-1)): good_tics.append(TIC)
    return

if __name__ == "__main__":
    main()