from pathlib import Path
import numpy as np
from bokeh.io import output_file, save, show

from gobo.corner_plot import create_multi_distribution_corner_plot


def load_states_params(states_path, burn_in=500000, param_indexes=None):
    data = np.loadtxt(states_path)
    data = np.atleast_2d(data)

    iteration = data[:, 0]
    samples = data[:, 2:24]

    if burn_in is not None:
        samples = samples[iteration >= burn_in]

    if param_indexes is not None:
        samples = samples[:, param_indexes]

    return samples


#states_path_1 = Path("old_sessions/2026_06_18_13_17_50_tic_id_302795956_sector_37/states.txt")
states_path_1 = Path("old_sessions/2026_06_13_21_30_00_tic_id_302795956_sector_37/states.txt")
states_path_2 = Path("old_sessions/2026_06_22_14_07_01_tic_id_302795956_sector_37/states.txt")

param_indexes_to_plot = list(range(22))
param_indexes_to_plot = [i for i in range(22) if i not in (2, 6, 21)]
#param_indexes_to_plot = [0,1,3,4,5,6]

array1 = load_states_params(states_path_1, burn_in=500000, param_indexes=param_indexes_to_plot)
array2 = load_states_params(states_path_2, burn_in=500000, param_indexes=param_indexes_to_plot)

print("Array 1 shape:", array1.shape)
print("Array 2 shape:", array2.shape)

label_map = {
    0: '$$logM1$$',
    1: '$$logM2$$',
    3: '$$e$$',
    4: '$$inclination$$',
    5: '$$omega0$$',
    7: '$$radius\\_resc\\_factor\\_of\\_star1$$',
    8: '$$radius\\_resc\\_factor\\_of\\_star2$$',
    9: '$$mu\\_1$$',
    10: '$$tau\\_1$$',
    11: '$$mu\\_2$$',
    12: '$$tau\\_2$$',
    13: '$$alpha\\_ref\\_1$$',
    14: '$$alpha\\_ref\\_2$$',
    15: '$$extra\\_alpha\\_beam\\_1$$',
    16: '$$extra\\_alpha\\_beam\\_2$$',
    17: '$$alpha\\_T1$$',
    18: '$$alpha\\_T2$$',
    19: '$$blending$$',
    20: '$$flux\\_tune$$',
}

dimension_labels = [
    label_map[i] for i in param_indexes_to_plot
]

corner_plot = create_multi_distribution_corner_plot(
    [array1, array2],
    dimension_labels=dimension_labels,
)

output_file("two_lc_corner_plot.html")
save(corner_plot)
show(corner_plot)