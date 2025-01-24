from pathlib import Path


def enforce_expected_data_directory_tree(data_directory_parent=Path('.')):
    data_directory = data_directory_parent.joinpath('data')
    data_directory.mkdir(exist_ok=True, parents=True)
    data_directory.joinpath('chains').mkdir(exist_ok=True, parents=True)
    data_directory.joinpath('lightcurves').mkdir(exist_ok=True, parents=True)
    data_directory.joinpath('lightcurves/folded_lightcurves').mkdir(exist_ok=True, parents=True)
    data_directory.joinpath('lightcurves/mcmc_lightcurves').mkdir(exist_ok=True, parents=True)
    data_directory.joinpath('magnitudes').mkdir(exist_ok=True, parents=True)
    data_directory.joinpath('pars').mkdir(exist_ok=True, parents=True)
    data_directory.joinpath('py_initialize').mkdir(exist_ok=True, parents=True)
