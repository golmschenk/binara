from pathlib import Path

import pandas as pd
import numpy.typing as npt
from bokeh.colors import Color
from bokeh.io import show, save
from bokeh.plotting import figure

from gobo.corner_plot import create_corner_plot
from gobo.internal.corner_plot import add_2d_histogram_to_figure
from gobo.internal.palette import default_discrete_palette
from bokeh.io import show, output_file

def create_2d_histogram_figure(
        array0: npt.NDArray,
        array1: npt.NDArray,
        *,
        color: Color = default_discrete_palette.blue
) -> figure:
    figure_ = figure()
    add_2d_histogram_to_figure(figure_, array0, array1, color=color)
    return figure_

chain_output_path = Path('data/chains/chain.118305806_sector_37_gmag_OMP_1.dat')
chain_output_path = Path('data/chains/chain.150284425_sector_all_gmag_OMP_1.dat')
chain_output_path = Path('data/chains/chain.220052771_sector_6_gmag_OMP_1.dat')

output_chain_data_frame = pd.read_csv(chain_output_path, delimiter=r'\s+', header=None)
column_indexes_to_compare = list(range(1, 23))  # Set the column indexes you want to compare.

iteration_column_index = 0
burn_in_iteration = 1_000_000
filtered_rows = output_chain_data_frame[
    output_chain_data_frame[iteration_column_index] >= burn_in_iteration]

filtered_data_frame = filtered_rows.loc[:, column_indexes_to_compare]
parameter_array = filtered_data_frame.values

print("Number of points shown in corner plot:", parameter_array.shape[0])
corner_plot = create_corner_plot(parameter_array, marginal_2d_figure_function=create_2d_histogram_figure)
save(corner_plot,"220052771_corner_plot_burnin_1_000_000.html")
show(corner_plot)
