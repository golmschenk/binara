from pathlib import Path
import re
import csv
import numpy as np

from bokeh.plotting import figure, output_file, save
from bokeh.models import ColumnDataSource, HoverTool


SESSIONS_DIR = Path("sessions")
filename = Path("input_data/best_periods.csv")
OUT_HTML = Path("period_vs_ecc.html")
common_date = "2026_06_08"

session_reg = re.compile(r"tic_id_(\d+)_sector_(\d+)", re.IGNORECASE)


def get_best_periods():
    best_periods = {}

    with open(filename, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tic = int(float(row["TIC ID"]))
            sector = int(float(row["Sector"]))
            period = float(row["Best Period (days)"])
            best_periods[(tic, sector)] = period

    return best_periods


def parse_tic_sector(folder_name):
    m = session_reg.search(folder_name)
    if not m:
        return None
    tic = int(m.group(1))
    sector = int(m.group(2))
    return tic, sector


def main():
    best_periods = get_best_periods()

    xs = []
    ys = []
    tics = []
    sectors = []
    folders = []

    for folder in sorted(SESSIONS_DIR.iterdir()):
        if not folder.is_dir():
            continue

        if common_date not in folder.name:
            continue

        parsed = parse_tic_sector(folder.name)
        if parsed is None:
            continue

        tic, sector = parsed

        if (tic, sector) not in best_periods:
            continue

        params_path = folder / "parameters.txt"
        if not params_path.exists():
            continue

        vals = np.loadtxt(params_path).flatten()
        if len(vals) < 4:
            continue

        ecc = float(vals[3])
        period = best_periods[(tic, sector)]

        xs.append(period)
        ys.append(ecc)
        tics.append(tic)
        sectors.append(sector)
        folders.append(folder.name)

    source = ColumnDataSource(data=dict(
        period=xs,
        ecc=ys,
        tic=tics,
        sector=sectors,
        folder=folders,
    ))

    p = figure(
        width=900,
        height=600,
        title="Period vs Eccentricity",
        x_axis_label="Best Period (days)",
        y_axis_label="Eccentricity",
        tools="pan,wheel_zoom,box_zoom,reset,save"
    )

    r = p.scatter("period", "ecc", source=source, size=7, alpha=0.8)

    hover = HoverTool(
        renderers=[r],
        tooltips=[
            ("TIC", "@tic"),
            ("Sector", "@sector"),
            ("Period", "@period{0.000000}"),
            ("Ecc", "@ecc{0.000000}"),
            ("Folder", "@folder"),
        ]
    )
    p.add_tools(hover)

    output_file(str(OUT_HTML))
    save(p)

    print(f"Saved {OUT_HTML}")
    print(f"Plotted {len(xs)} points")


if __name__ == "__main__":
    main()