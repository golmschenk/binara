from pathlib import Path

import pandas as pd
from bokeh.io import show
from bokeh.plotting import figure

output_path = Path('sessions/2026_02_13_11_31_00_tic_id_118305806_sector_37/folded_observed_and_model_light_curves.txt')
data_frame = pd.read_csv(output_path, sep='\s+', header=None, index_col=False)

model_phases = data_frame[0]
observed_fluxes = data_frame[1]
model_fluxes = data_frame[2]

figure_ = figure()
figure_.scatter(x=model_phases, y=model_fluxes, legend_label="Model", color='firebrick')
figure_.scatter(x=model_phases, y=observed_fluxes, legend_label="Observed", color='mediumblue')
show(figure_)