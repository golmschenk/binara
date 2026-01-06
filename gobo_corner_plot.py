from pathlib import Path

import pandas as pd
from bokeh.io import show
from gobo.corner_plot import create_corner_plot

chain_output_path = Path('data/chains/chain.110602878_sector_34_gmag_OMP_1.dat')
example_output_chain_data_frame = pd.read_csv(chain_output_path, delimiter=r'\s+', header=None)
column_indexes_to_compare = [2, 3, 5, 7]  # Set the column indexes you want to compare.

iteration_column_index = 0
burn_in_iteration = 2_200_000
filtered_rows = example_output_chain_data_frame[
    example_output_chain_data_frame[iteration_column_index] >= burn_in_iteration]

filtered_data_frame = filtered_rows.loc[:, column_indexes_to_compare]
parameter_array = filtered_data_frame.values

#burn_in_states_to_remove = 1000
# filtered_data_frame = example_output_chain_data_frame.iloc[burn_in_states_to_remove:, column_indexes_to_compare]
# parameter_array = filtered_data_frame.values
#
print("Number of points shown in corner plot:", parameter_array.shape[0])
corner_plot = create_corner_plot(parameter_array)
show(corner_plot)