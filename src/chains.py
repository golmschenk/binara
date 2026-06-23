import numpy as np
from pathlib import Path

from bokeh.io import output_file, save, show
from bokeh.layouts import column
from bokeh.plotting import figure

chains = [
    (
        "122224576_sector_54",
        Path("scripts/2026_03_01_00_00_15_tic_id_122224576_sector_54/states.txt"),
    )
]

plots = []

for title, chain_path in chains:
    data = np.loadtxt(chain_path)
    data = data[data[:, 0] != 0]

    iters = data[:, 0]
    vals  = data[:, 1]

    p = figure(
        title=f"{title} Iteration vs Log Likelihood",
        x_axis_label="Iteration",
        y_axis_label="Log Likelihood",
        height=350,
        sizing_mode="stretch_width",
        tools="pan,wheel_zoom,box_zoom,reset,save",
    )

    p.line(iters, vals, line_width=1, color="navy", alpha=0.8)
    plots.append(p)

layout = column(*plots, sizing_mode="stretch_width")

save(layout, "chains.html")
show(layout)
