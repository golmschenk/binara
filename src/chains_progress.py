import numpy as np
import matplotlib.pyplot as plt
from bokeh.io import output_file, save
from bokeh.layouts import column
from bokeh.plotting import figure
from bokeh.io import show, save

data = np.loadtxt("data/chains/chain.150284425_sector_all_gmag_OMP_1.dat")
#data = np.loadtxt("data/chains/chain.118305806_sector_37_gmag_OMP_1.dat")
#data = np.loadtxt("data/chains/chain.220052771_sector_6_gmag_OMP_1.dat")

data = data[data[:, 0] != 0]
# data = data[data[:, 0] <= 175_000]
iters = data[:, 0]
vals  = data[:, 1]

p = figure(
    title="Iteration vs. Log Likelihood",
    x_axis_label="Iteration",
    y_axis_label="Log Likelihood",
    width=900,
    height=500,
    tools="pan,wheel_zoom,box_zoom,reset,save"
)

p.line(iters, vals, line_width=1, color="navy", alpha=0.8)

output_file("loglikelihood_chain_110602878_sector_34.html")
save(p)
show(p)