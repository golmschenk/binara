from pathlib import Path
from collections import Counter
import re

import numpy as np
import lightkurve as lk
from astropy.timeseries import BoxLeastSquares

from bokeh.plotting import figure, output_file, save
from bokeh.layouts import column
from bokeh.models import Div


REPEAT_PATH = Path("src/repeat_heartbeats.txt")
INPUT_ROOT = Path("input_data")
OUT_HTML = Path("repeat_tics_single_vs_multi.html")

# old corrected/default Lightkurve behavior
FLUX_COLUMN = None

# for SAP instead, use:
# FLUX_COLUMN = "sap_flux"

SINGLE_SECTOR_NBINS = 150

reg = re.compile(r"phot_(\d+)-s(\d{4})", re.IGNORECASE)


def extract_tic_sector_from_line(line):
    filename = line.strip().replace("\\", "/").split("/")[-1]
    m = reg.search(filename)

    if not m:
        return None

    tic_id = int(m.group(1))
    sector = int(m.group(2))

    return tic_id, sector, filename


def get_repeat_tic_entries():
    entries = []

    with open(REPEAT_PATH, "r") as f:
        for line in f:
            if not line.strip():
                continue

            result = extract_tic_sector_from_line(line)
            if result is not None:
                entries.append(result)

    counts = Counter(tic_id for tic_id, sector, filename in entries)

    repeat_entries = [
        (tic_id, sector, filename)
        for tic_id, sector, filename in entries
        if counts[tic_id] > 1
    ]

    seen = set()
    repeat_entries_unique = []

    for tic_id, sector, filename in repeat_entries:
        if tic_id in seen:
            continue
        seen.add(tic_id)
        repeat_entries_unique.append((tic_id, sector, filename))

    return repeat_entries_unique


def read_existing_multi_folded(tic_id, sector):
    folded_path = INPUT_ROOT / f"tic_id_{tic_id}_sector_{sector}" / "folded_observed_light_curve.txt"

    if not folded_path.exists():
        print(f"Missing input folded file: {folded_path}", flush=True)
        return None

    with open(folded_path, "r") as f:
        header = f.readline().strip().split()

    n_header = int(float(header[0]))
    period = float(header[1])

    data = np.loadtxt(folded_path, skiprows=1)

    if data.ndim == 1:
        data = data.reshape(1, -1)

    return folded_path, n_header, period, data


def load_and_normalize_one_sector(tic, sector):
    search = lk.search_lightcurve(f"TIC {tic}", mission="TESS", sector=sector)

    if len(search) == 0:
        print(f"TIC {tic}, sector {sector}: no light curve found", flush=True)
        return None

    lc = search[0].download()

    if lc is None:
        print(f"TIC {tic}, sector {sector}: download failed", flush=True)
        return None

    lc = lc.remove_nans()

    t = np.asarray(lc.time.value, dtype=np.float64)
    y = np.asarray(lc.flux.value, dtype=np.float64)
    yerr = np.asarray(lc.flux_err.value, dtype=np.float64)

    good = np.isfinite(t) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)

    t = t[good]
    y = y[good]
    yerr = yerr[good]

    if len(t) == 0:
        return None

    med = np.nanmedian(y)

    if not np.isfinite(med) or med == 0:
        return None

    y_norm = y / med
    yerr_norm = yerr / abs(med)

    return np.column_stack([t, y_norm, yerr_norm])


def find_best_period_from_array(arr, tic_id):
    t = arr[:, 0]
    y = arr[:, 1]
    yerr = arr[:, 2]

    bls = BoxLeastSquares(t, y, dy=yerr)

    periods = np.linspace(1, 30, 10000)
    bls_results = bls.power(periods, 0.1)

    best_idx = np.nanargmax(bls_results.power)
    best_period = float(bls_results.period[best_idx])

    print(f"TIC {tic_id}: single-sector best period = {best_period:.6f} days", flush=True)

    return best_period


def bin_error_from_linear_residuals(x, y):
    if len(x) < 3:
        return np.std(y, ddof=1) if len(y) > 1 else np.nan

    a, b = np.polyfit(x, y, 1)
    residuals = y - (a * x + b)

    return np.std(residuals, ddof=1)


def fold_array_same_way(combined_arr, period, nbins):
    time = combined_arr[:, 0]
    flux = combined_arr[:, 1]
    err = combined_arr[:, 2]

    t0 = time.min()
    phase = ((time - t0) / period) % 1.0

    bin_edges = np.linspace(0.0, 1.0, nbins + 1)
    bin_centers_phase = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    binned_phase_days = bin_centers_phase * period

    rows = []

    for i in range(nbins):
        in_bin = (phase >= bin_edges[i]) & (phase < bin_edges[i + 1])
        idx = np.where(in_bin)[0]
        n = len(idx)

        if n == 0:
            continue

        fl_bin = flux[idx]
        ph_bin = phase[idx]

        f_mean = np.mean(fl_bin)

        if n >= 3:
            f_err = bin_error_from_linear_residuals(ph_bin, fl_bin)
        elif n == 2:
            f_err = np.std(fl_bin, ddof=1)
        else:
            f_err = err[idx][0]

        if np.isfinite(f_mean) and np.isfinite(f_err):
            rows.append((binned_phase_days[i], f_mean, f_err))

    if len(rows) == 0:
        return None

    return np.asarray(rows, dtype=float)


def make_plot(tic_id, sector, multi_period, single_period, multi_data, single_data):
    p = figure(
        width=1000,
        height=350,
        title=f"TIC {tic_id} sector {sector} | multi P={multi_period:.6f} | single P={single_period:.6f}",
        x_axis_label="Folded time [days]",
        y_axis_label="Normalized flux",
        tools="pan,wheel_zoom,box_zoom,reset,save"
    )

    p.circle(
        multi_data[:, 0],
        multi_data[:, 1],
        size=4,
        alpha=0.65,
        legend_label="multi"
    )

    p.circle(
        single_data[:, 0],
        single_data[:, 1],
        size=4,
        alpha=0.65,
        color="red",
        legend_label="single"
    )

    p.legend.location = "top_right"
    p.legend.click_policy = "hide"

    return p

def main():
    repeat_entries = get_repeat_tic_entries()
    plots = []
    made = 0
    skipped = 0

    for tic_id, sector, filename in repeat_entries:
        print(f"\nWorking on TIC {tic_id}, sector {sector}", flush=True)

        multi = read_existing_multi_folded(tic_id, sector)

        if multi is None:
            skipped += 1
            continue

        multi_path, n_header, multi_period, multi_data = multi

        single_arr = load_and_normalize_one_sector(tic_id, sector)

        if single_arr is None:
            skipped += 1
            continue

        single_period = find_best_period_from_array(single_arr, tic_id)

        single_data = fold_array_same_way(
            combined_arr=single_arr,
            period=single_period,
            nbins=SINGLE_SECTOR_NBINS
        )

        if single_data is None:
            skipped += 1
            continue

        plots.append(
            make_plot(
                tic_id=tic_id,
                sector=sector,
                multi_period=multi_period,
                single_period=single_period,
                multi_data=multi_data,
                single_data=single_data
            )
        )

        made += 1

    output_file(str(OUT_HTML))
    save(column(*plots))

    print(f"\nDone. Saved HTML: {OUT_HTML}", flush=True)
    print(f"Made {made} plots. Skipped {skipped} entries.", flush=True)


if __name__ == "__main__":
    main()