from pathlib import Path
import re
import math
import csv

import numpy as np
from astropy.io import fits

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, HoverTool
from bokeh.palettes import Category20, Turbo256
from bokeh.resources import INLINE
from bokeh.embed import file_html
from bokeh.layouts import column


TIC_SECTOR_RE = re.compile(r"phot_(\d+)-s(\d{4})", re.IGNORECASE)
SECTOR_RE = re.compile(r"-s(\d{4})", re.IGNORECASE)

ROOT = Path("mastDownload/HLSP")
REPEAT_FILE = Path("src/repeat_heartbeats.txt")

FLUX_TYPE = "SAP"
PLOTS_PER_PAGE = 50

OUTPUT_ROOT = Path("stitched_all_tics")
TXT_DIR = OUTPUT_ROOT / "stitched_txt"
HTML_DIR = OUTPUT_ROOT / "html_pages"

# This only down-samples the HTML display so the files don't become insane.
# The saved .txt still contains all stitched points.
MAX_POINTS_PER_SECTOR_FOR_HTML = 2500


def extract_tic_sector_from_name(filename):
    m = TIC_SECTOR_RE.search(filename)
    if not m:
        return None
    tic_id = int(m.group(1))
    sector = int(m.group(2))
    return tic_id, sector


def extract_sector(filename):
    m = SECTOR_RE.search(filename)
    if not m:
        return -1
    return int(m.group(1))


def read_unique_tics_from_repeat_file(path):
    unique_tics = []
    seen = set()

    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            filename = line.replace("\\", "/").split("/")[-1]
            parsed = extract_tic_sector_from_name(filename)

            if parsed is None:
                continue

            tic_id, sector = parsed

            if tic_id not in seen:
                unique_tics.append(tic_id)
                seen.add(tic_id)

    print(f"Found {len(unique_tics)} unique TICs in {path}", flush=True)
    return unique_tics


def build_file_index_for_tics(root, tic_ids):
    tic_set = set(int(t) for t in tic_ids)
    files_by_tic = {int(t): [] for t in tic_ids}

    print(f"Scanning {root} once for matching TIC FITS files...", flush=True)

    n_scanned = 0
    n_matched = 0

    for path in root.rglob("*_lc.fits"):
        n_scanned += 1

        if n_scanned % 10000 == 0:
            print(f"Scanned {n_scanned} FITS files...", flush=True)

        parsed = extract_tic_sector_from_name(path.name)
        if parsed is None:
            continue

        tic_id, sector = parsed

        if tic_id in tic_set:
            files_by_tic[tic_id].append(path)
            n_matched += 1

    for tic_id in files_by_tic:
        files_by_tic[tic_id] = sorted(
            files_by_tic[tic_id],
            key=lambda p: extract_sector(p.name)
        )

    print(f"Done scanning. Scanned {n_scanned} files, matched {n_matched}.", flush=True)

    return files_by_tic


def read_normalized_sector(path, flux_type="SAP"):
    if flux_type.upper() == "SAP":
        flux_col = "SAP_FLUX"
        err_col = "SAP_FLUX_ERR"
    elif flux_type.upper() == "PDCSAP":
        flux_col = "PDCSAP_FLUX"
        err_col = "PDCSAP_FLUX_ERR"
    else:
        raise ValueError("flux_type must be SAP or PDCSAP")

    sector = extract_sector(path.name)

    with fits.open(path) as hdul:
        data = hdul[1].data

        time = np.asarray(data["TIME"], dtype=float)
        flux = np.asarray(data[flux_col], dtype=float)
        flux_err = np.asarray(data[err_col], dtype=float)

        if "QUALITY" in data.columns.names:
            quality = np.asarray(data["QUALITY"], dtype=int)
            good_quality = quality == 0
        else:
            good_quality = np.ones(len(time), dtype=bool)

    good = (
            np.isfinite(time)
            & np.isfinite(flux)
            & np.isfinite(flux_err)
            & (flux_err > 0)
            & good_quality
    )

    time = time[good]
    flux = flux[good]
    flux_err = flux_err[good]

    if len(time) == 0:
        print(f"Sector {sector}: no good points after filtering", flush=True)
        return None

    med = np.nanmedian(flux)

    if not np.isfinite(med) or med == 0:
        print(f"Sector {sector}: bad median flux", flush=True)
        return None

    flux_norm = flux / med
    flux_err_norm = flux_err / abs(med)

    arr = np.column_stack([
        time,
        flux_norm,
        flux_err_norm,
        np.full(len(time), sector)
    ])

    return arr


def stitch_tic_from_files(tic_id, files, flux_type="PDCSAP"):
    arrays = []

    for path in files:
        arr = read_normalized_sector(path, flux_type=flux_type)
        if arr is not None:
            arrays.append(arr)

    if len(arrays) == 0:
        return None, []

    stitched = np.vstack(arrays)
    stitched = stitched[np.argsort(stitched[:, 0])]

    sectors = sorted(set(stitched[:, 3].astype(int)))

    return stitched, sectors


def get_sector_color(i, n):
    if n <= 20:
        return Category20[20][i % 20]
    return Turbo256[int(i * 255 / max(n - 1, 1))]


def downsample_for_html(time, flux, err, sector_col):
    if MAX_POINTS_PER_SECTOR_FOR_HTML is None:
        return time, flux, err, sector_col

    if len(time) <= MAX_POINTS_PER_SECTOR_FOR_HTML:
        return time, flux, err, sector_col

    idx = np.linspace(
        0,
        len(time) - 1,
        MAX_POINTS_PER_SECTOR_FOR_HTML,
        dtype=int
    )

    return time[idx], flux[idx], err[idx], sector_col[idx]


def make_stitched_figure(tic_id, stitched, sectors, flux_type="SAP"):
    time = stitched[:, 0]
    flux = stitched[:, 1]
    err = stitched[:, 2]
    sector_col = stitched[:, 3].astype(int)

    p = figure(
        title=f"TIC {tic_id} stitched {flux_type}; sectors {sectors}",
        x_axis_label="Time [BTJD]",
        y_axis_label="Normalized flux",
        width=1100,
        height=260,
        tools="pan,wheel_zoom,box_zoom,reset,save"
    )

    renderers = []

    for i, sector in enumerate(sectors):
        mask = sector_col == sector

        t_s = time[mask]
        f_s = flux[mask]
        e_s = err[mask]
        sec_s = sector_col[mask]

        t_s, f_s, e_s, sec_s = downsample_for_html(t_s, f_s, e_s, sec_s)

        source = ColumnDataSource({
            "time": t_s,
            "flux": f_s,
            "err": e_s,
            "sector": sec_s,
        })

        r = p.scatter(
            x="time",
            y="flux",
            source=source,
            size=2,
            alpha=0.6,
            color=get_sector_color(i, len(sectors)),
            legend_label=f"S{sector}"
        )

        renderers.append(r)

    hover = HoverTool(
        renderers=renderers,
        tooltips=[
            ("sector", "@sector"),
            ("time", "@time{0.000000}"),
            ("flux", "@flux{0.000000}"),
            ("err", "@err{0.000000}"),
        ],
    )
    p.add_tools(hover)

    p.legend.location = "top_left"
    p.legend.click_policy = "hide"

    return p


def save_stitched_txt(tic_id, stitched, flux_type="SAP"):
    TXT_DIR.mkdir(parents=True, exist_ok=True)

    txt_path = TXT_DIR / f"tic_{tic_id}_{flux_type.lower()}_stitched.txt"

    np.savetxt(
        txt_path,
        stitched,
        fmt="%.10f %.10f %.10f %d",
        header="time normalized_flux normalized_flux_err sector"
    )

    return txt_path


def save_html_page(figures, page_num, flux_type="SAP"):
    HTML_DIR.mkdir(parents=True, exist_ok=True)

    page_path = HTML_DIR / f"stitched_{flux_type.lower()}_page_{page_num}.html"

    layout = column(*figures)

    html = file_html(
        layout,
        INLINE,
        title=f"Stitched {flux_type} light curves page {page_num}"
    )

    page_path.write_text(html)

    print(f"Wrote HTML page: {page_path}", flush=True)

    return page_path


def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)

    unique_tics = read_unique_tics_from_repeat_file(REPEAT_FILE)

    files_by_tic = build_file_index_for_tics(
        root=ROOT,
        tic_ids=unique_tics
    )

    figures = []
    page_num = 1
    summary_rows = []

    for i, tic_id in enumerate(unique_tics, start=1):
        files = files_by_tic.get(tic_id, [])

        print(
            f"\n[{i}/{len(unique_tics)}] TIC {tic_id}: found {len(files)} files",
            flush=True
        )

        if len(files) == 0:
            print(f"Skipping TIC {tic_id}: no files found", flush=True)
            summary_rows.append([tic_id, 0, "", "NO_FILES"])
            continue

        stitched, sectors = stitch_tic_from_files(
            tic_id=tic_id,
            files=files,
            flux_type=FLUX_TYPE
        )

        if stitched is None:
            print(f"Skipping TIC {tic_id}: no usable data", flush=True)
            summary_rows.append([tic_id, len(files), "", "NO_USABLE_DATA"])
            continue

        txt_path = save_stitched_txt(
            tic_id=tic_id,
            stitched=stitched,
            flux_type=FLUX_TYPE
        )

        fig = make_stitched_figure(
            tic_id=tic_id,
            stitched=stitched,
            sectors=sectors,
            flux_type=FLUX_TYPE
        )

        figures.append(fig)

        summary_rows.append([
            tic_id,
            len(files),
            ",".join(str(s) for s in sectors),
            "OK"
        ])

        print(f"Saved stitched txt: {txt_path}", flush=True)

        if len(figures) == PLOTS_PER_PAGE:
            save_html_page(figures, page_num, flux_type=FLUX_TYPE)
            page_num += 1
            figures = []

    if len(figures) > 0:
        save_html_page(figures, page_num, flux_type=FLUX_TYPE)

    summary_path = OUTPUT_ROOT / f"stitched_{FLUX_TYPE.lower()}_summary.csv"

    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["tic_id", "n_files", "sectors", "status"])
        writer.writerows(summary_rows)

    print(f"\nDone. Summary: {summary_path}", flush=True)
    print(f"HTML pages in: {HTML_DIR}", flush=True)
    print(f"Stitched txt files in: {TXT_DIR}", flush=True)


if __name__ == "__main__":
    main()