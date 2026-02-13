from pathlib import Path

import pandas as pd
from bokeh.io import show
from bokeh.plotting import figure

output_path = Path('data/lightcurves/mcmc_lightcurves/150284425_sector_all_gmag_OMP_1.out')
output_path = Path('data/lightcurves/mcmc_lightcurves/118305806_sector_37_gmag_OMP_1.out')
output_path = Path('data/lightcurves/mcmc_lightcurves/220052771_sector_6_gmag_OMP_1.out')

folded_path = Path('data/lightcurves/folded_lightcurves/150284425_sector_all.txt')
folded_path = Path('data/lightcurves/folded_lightcurves/118305806_sector_37.txt')
folded_path = Path('data/lightcurves/folded_lightcurves/220052771_sector_6.txt')

data_frame = pd.read_csv(output_path, sep='\s+', header=None, index_col=False)
folded_df = pd.read_csv(folded_path, sep=r'\s+', header=None, skiprows=1)

model_phases = data_frame[0]
observed_fluxes = data_frame[1]
model_fluxes = data_frame[2]

flux_errors = folded_df[2]

figure_ = figure()
figure_.scatter(x=model_phases, y=model_fluxes, legend_label="Model", color='firebrick')
figure_.scatter(x=model_phases, y=observed_fluxes, legend_label="Observed", color='mediumblue')

figure_.segment(
    x0=model_phases,
    y0=observed_fluxes - flux_errors,
    x1=model_phases,
    y1=observed_fluxes + flux_errors,
    color='mediumblue'
)

show(figure_)
