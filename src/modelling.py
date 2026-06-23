from pathlib import Path
import numpy as np
from bokeh.layouts import column
from bokeh.plotting import figure, show

def plot_mcmc(path):
    arr = np.loadtxt(path)
    phase, data, model = arr[:,0], arr[:,1], arr[:,2]
    err = arr[:,3] if arr.shape[1] >= 4 else None

    # clean + sort by phase for a nice line
    m = np.isfinite(phase) & np.isfinite(data) & np.isfinite(model)
    if err is not None: m &= np.isfinite(err)
    phase, data, model = phase[m], data[m], model[m]
    if err is not None: err = err[m]
    s = np.argsort(phase)
    phase, data, model = phase[s], data[s], model[s]
    if err is not None: err = err[s]

    # model is ~1 ⇒ scale to data median
    model *= np.median(data)

    p = figure(title=Path(path).name, x_axis_label="Phase", y_axis_label="Flux",
               height=350, sizing_mode="stretch_width")
    p.scatter(phase, data, size=3, alpha=0.6, legend_label="Data")
    p.line(phase, model, alpha=0.9, legend_label="Model")
    p.legend.location = "top_left"

    # residuals
    r = data - model
    p2 = figure(x_axis_label="Phase", y_axis_label="Residual", height=200,
                sizing_mode="stretch_width", x_range=p.x_range)
    p2.scatter(phase, r, size=3, alpha=0.6)

    show(column(p, p2))

plot_mcmc("data/lightcurves/mcmc_lightcurves/110602878_sector_34_gmag_OMP_1.out")
