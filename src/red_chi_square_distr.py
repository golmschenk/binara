from pathlib import Path
import numpy as np
import pandas as pd
from bokeh.layouts import column

from bokeh.plotting import figure, output_file, save
from bokeh.models import HoverTool, ColumnDataSource


SESSIONS_DIR = Path("sessions")
SESSIONS_DIR = Path("old_sessions")

LC_FILENAME = "folded_observed_and_model_light_curves.txt"

N_FIT_PARAMS = 22


def reduced_chi_square(path: Path, n_fit_params: int = 0):
    data = np.loadtxt(path)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    time = data[:, 0]
    observed = data[:, 1]
    model = data[:, 2]
    error = data[:, 3]

    good = (
            np.isfinite(time)
            & np.isfinite(observed)
            & np.isfinite(model)
            & np.isfinite(error)
            & (error > 0)
    )

    observed = observed[good]
    model = model[good]
    error = error[good]

    n_points = len(observed)

    if n_points == 0:
        return np.nan, np.nan, 0

    chi2 = np.sum(((observed - model) / error) ** 2)

    dof = n_points - n_fit_params
    red_chi2 = chi2 / dof

    return chi2, red_chi2, n_points

def main():
    session_dirs = sorted([
        p for p in SESSIONS_DIR.iterdir()
        if p.is_dir() and "06_07" in p.name
    ])
    rows = []

    print(f"Found {len(session_dirs)} session folders to check.")

    for session_dir in session_dirs:
        lc_path = session_dir / LC_FILENAME

        if not lc_path.exists():
            print(f"Skipping {session_dir.name}: missing {LC_FILENAME}")
            continue

        try:
            chi2, red_chi2, n_points = reduced_chi_square(
                lc_path,
                n_fit_params=N_FIT_PARAMS
            )

            rows.append({
                "session": session_dir.name,
                "path": str(lc_path),
                "chi2": chi2,
                "reduced_chi2": red_chi2,
                "n_points": n_points,
                "n_fit_params": N_FIT_PARAMS,
            })

            print(
                f"{session_dir.name}: "
                f"chi2={chi2:.6g}, red_chi2={red_chi2:.6g}, N={n_points}"
            )

        except Exception as e:
            print(f"Failed for {session_dir.name}: {e}")

    df = pd.DataFrame(rows)

    if df.empty:
        print("No valid light curve files found.")
        return

    csv_out = Path("red_chi2.csv")
    df.to_csv(csv_out, index=False)
    print(f"Wrote {csv_out}")

    values = df["reduced_chi2"].replace([np.inf, -np.inf], np.nan).dropna()

    hist, edges = np.histogram(values, bins=40)

    source = ColumnDataSource(data={
        "left": edges[:-1],
        "right": edges[1:],
        "top": hist,
        "bottom": np.zeros_like(hist),
    })

    p = figure(
        title=f"Reduced chi-square distribution,",
        x_axis_label="Reduced chi-square",
        y_axis_label="Count",
        width=900,
        height=500,
        tools="pan,wheel_zoom,box_zoom,reset,save"
    )

    p.quad(
        source=source,
        left="left",
        right="right",
        top="top",
        bottom="bottom",
        line_color="black",
        fill_alpha=0.7
    )

    p.add_tools(HoverTool(
        tooltips=[
            ("bin left", "@left{0.000}"),
            ("bin right", "@right{0.000}"),
            ("count", "@top"),
        ]
    ))

    # Make log-scale x-axis histogram too
    positive_values = values[values > 0]

    if len(positive_values) > 0:
        log_min = positive_values.min()
        log_max = positive_values.max()

        if log_min == log_max:
            log_min = log_min / 10
            log_max = log_max * 10

        log_edges = np.logspace(
            np.log10(log_min),
            np.log10(log_max),
            40
        )

        log_hist, log_edges = np.histogram(positive_values, bins=log_edges)

        log_source = ColumnDataSource(data={
            "left": log_edges[:-1],
            "right": log_edges[1:],
            "top": log_hist,
            "bottom": np.zeros_like(log_hist),
        })

        p_log = figure(
            title="Reduced chi-square distribution, log x-axis",
            x_axis_label="Reduced chi-square",
            y_axis_label="Count",
            width=900,
            height=500,
            x_axis_type="log",
            tools="pan,wheel_zoom,box_zoom,reset,save"
        )

        p_log.quad(
            source=log_source,
            left="left",
            right="right",
            top="top",
            bottom="bottom",
            line_color="black",
            fill_alpha=0.7
        )

        p_log.add_tools(HoverTool(
            tooltips=[
                ("bin left", "@left{0.000}"),
                ("bin right", "@right{0.000}"),
                ("count", "@top"),
            ]
        ))

        layout = column(p, p_log)

    else:
        print("No positive reduced chi-square values, so skipping log plot.")
        layout = p


    html_out = "red_chi2.html"
    output_file(html_out)
    save(layout)

    print(f"Wrote {html_out}")
    print()
    print("Summary:")
    print(df["reduced_chi2"].describe())


if __name__ == "__main__":
    main()