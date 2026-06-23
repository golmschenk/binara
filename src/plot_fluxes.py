from pathlib import Path

import pandas as pd
from bokeh.io import output_file, save, show
from bokeh.layouts import column
from bokeh.plotting import figure

pairs = [
    (
        "150284425 sector_all",
        Path('corrected_original_input_data/lightcurves/mcmc_lightcurves/150284425_sector_all_gmag_OMP_1.out'),
        Path('corrected_original_input_data/lightcurves/folded_lightcurves/150284425_sector_all.txt'),
    ),
    (
        "118305806 sector_37",
        Path('corrected_original_input_data/lightcurves/mcmc_lightcurves/118305806_sector_37_gmag_OMP_1.out'),
        Path('corrected_original_input_data/lightcurves/folded_lightcurves/118305806_sector_37.txt'),
    ),
    (
        "220052771 sector_6",
        Path('corrected_original_input_data/lightcurves/mcmc_lightcurves/220052771_sector_6_gmag_OMP_1.out'),
        Path('corrected_original_input_data/lightcurves/folded_lightcurves/220052771_sector_6.txt'),
    ),
]

plots = []

for title, output_path, folded_path in pairs:
    df = pd.read_csv(output_path, sep=r'\s+', header=None)
    folded_df = pd.read_csv(folded_path, sep=r'\s+', header=None, skiprows=1)

    model_phases = df[0]
    observed_fluxes = df[1]
    model_fluxes = df[2]
    flux_errors = folded_df[2]

    n = min(len(model_phases), len(observed_fluxes), len(model_fluxes), len(flux_errors))
    model_phases = model_phases.iloc[:n]
    observed_fluxes = observed_fluxes.iloc[:n]
    model_fluxes = model_fluxes.iloc[:n]
    flux_errors = flux_errors.iloc[:n]

    p = figure(title=title, height=350, sizing_mode="stretch_width")
    p.scatter(model_phases, model_fluxes, legend_label="Model", color="firebrick")
    p.scatter(model_phases, observed_fluxes, legend_label="Observed", color="mediumblue")
    p.segment(
        x0=model_phases,
        y0=observed_fluxes - flux_errors,
        x1=model_phases,
        y1=observed_fluxes + flux_errors,
        color="mediumblue",
    )
    p.legend.location = "top_left"
    plots.append(p)

layout = column(*plots, sizing_mode="stretch_width")

save(layout,"fluxes.html")
show(layout)
