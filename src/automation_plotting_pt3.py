import numpy as np
from pathlib import Path
from bokeh.io import show, save, output_file
from bokeh.layouts import column
from bokeh.plotting import figure
from binara.binara_ext import internal_calculate_light_curve, EclipseMethod
import tomllib
from bokeh.models import Div, Button, CustomJS
import lightkurve as lk

SESSIONS_DIR = Path("old_sessions")
PAGE_SIZE = 50

PARAM_NAMES = [
    "logM1",
    "logM2",
    "logPdays",
    "e",
    "cos_inc",
    "omega0",
    "T0",
    "rr1",
    "rr2",
    "mu_1",
    "tau_1",
    "mu_2",
    "tau_2",
    "alpha_ref_1",
    "alpha_ref_2",
    "extra_alpha_beam_1",
    "extra_alpha_beam_2",
    "alpha_Teff_1",
    "alpha_Teff_2",
    "blending",
    "flux_tune",
]

def load_eclipse_method(session_dir: Path):
    config_file = session_dir / "configuration.toml"

    if not config_file.exists():
        print(f"{session_dir.name}: no configuration.toml, using REGULAR")
        return EclipseMethod.REGULAR

    with open(config_file, "rb") as f:
        config = tomllib.load(f)

    method = config.get("modeling", {}).get("eclipse_method", "regular")
    method = method.lower().strip()

    if method == "regular":
        return EclipseMethod.REGULAR
    elif method == "limb_darkening":
        return EclipseMethod.LIMB_DARKENING

def load_params(session_dir: Path):
    param_file = session_dir / "parameters.txt"

    if not param_file.exists():
        return None

    for line in param_file.read_text().splitlines():
        line = line.strip()

        if line == "" or line.startswith("#"):
            continue

        nums = line.split()

        if len(nums) >= 22:
            params = np.array(nums[:22], dtype=np.float64)
            params[4] = np.cos(params[4])

            return np.ascontiguousarray(params, dtype=np.float64)

        if len(nums) >= 21:
            params = np.array(nums[:21], dtype=np.float64)
            params[4] = np.cos(params[4])

            return np.ascontiguousarray(params, dtype=np.float64)

    return None

def make_param_table(session_dir: Path):
    param_file = session_dir / "parameters.txt"
    pars = None

    for line in param_file.read_text().splitlines():
        line = line.strip()

        if line == "" or line.startswith("#"):
            continue

        nums = line.split()

        if len(nums) >= 21:
            pars = nums[:21]
            break

    if pars is None:
        return Div(
            text="<p><b>No 21-parameter line found in parameters.txt.</b></p>",
            width=900
        )

    rows = ""

    for name, value in zip(PARAM_NAMES, pars):
        rows += f"""
        <tr>
            <td style="border: 1px solid black; padding: 6px;">{name}</td>
            <td style="border: 1px solid black; padding: 6px;">{value}</td>
        </tr>
        """

    table_html = f"""
    <h3>Parameters</h3>
    <table style="border-collapse: collapse; width: 500px;">
        <tr>
            <th style="border: 1px solid black; padding: 6px;">Parameter</th>
            <th style="border: 1px solid black; padding: 6px;">Value</th>
        </tr>
        {rows}
    </table>
    """

    return Div(text=table_html, width=900)

def produce_residual_plot_from_file(lc_file: Path, linked_plot):
    data = np.loadtxt(lc_file)

    if data.size == 0:
        print(f"Skipping {lc_file.parent.name}: no data in {lc_file.name}")
        return None

    data = np.atleast_2d(data)
    time = data[:, 0]
    flux1 = data[:, 1]
    flux2 = data[:, 2]
    err1 = data[:, 3]
    residuals = flux1 - flux2

    p = figure(
        title=f"Residuals: Observed - Model ({lc_file.parent.name})",
        x_axis_label="Time / Phase",
        y_axis_label="Residual",
        width=900,
        height=250,
        tools="pan,wheel_zoom,box_zoom,reset,save",
        x_range=linked_plot.x_range
    )
    p.scatter(time, residuals, size=5, color="darkgreen", alpha=0.8)
    p.segment(
        x0=time,
        y0=residuals - err1,
        x1=time,
        y1=residuals + err1,
        color="darkgreen",
        alpha=0.7
    )
    p.line([time.min(), time.max()], [0, 0], line_width=1, color="black", alpha=0.8)
    p.grid.grid_line_alpha = 0.3
    return p
def produce_light_curve_from_file(lc_file: Path, session_dir: Path):
    data = np.loadtxt(lc_file)

    if data.size == 0:
        print(f"Skipping {lc_file.parent.name}: no data in {lc_file.name}")
        return None

    data = np.atleast_2d(data)

    time = data[:, 0]
    flux1 = data[:, 1]
    flux2 = data[:, 2]
    err1 = data[:, 3]

    p = figure(
        title=f"Folded Light Curve: Observed vs Model ({lc_file.parent.name})",
        x_axis_label="Time / Phase",
        y_axis_label="Normalized Flux",
        width=900,
        height=450,
        tools="pan,wheel_zoom,box_zoom,reset,save"
    )

    params = load_params(session_dir)
    eclipse_method = load_eclipse_method(session_dir)

    if 1:
        params = np.asarray(params, dtype=np.float64)

        smooth_x = np.linspace(time.min(), time.max(), 1000, dtype=np.float64)

        smooth_flux = internal_calculate_light_curve(
            smooth_x,
            params,
            eclipse_method
        )

        p.line(
            smooth_x,
            smooth_flux,
            line_width=2,
            color="mediumblue",
            alpha=0.9,
            legend_label="Model"
        )
    p.scatter(time, flux1, size=5, color="firebrick", alpha=0.8, legend_label="Observed")

    p.segment(
        x0=time,
        y0=flux1 - err1,
        x1=time,
        y1=flux1 + err1,
        color="firebrick",
        alpha=0.7
    )

    p.legend.click_policy = "hide"
    p.grid.grid_line_alpha = 0.3
    return p


def produce_chains_plot(chain_file: Path):
    try:
        data = np.loadtxt(chain_file)
    except Exception as e:
        print(f"Could not load {chain_file}: {e}")
        return None

    if data.size == 0:
        print(f"Skipping {chain_file.parent.name}: no data in {chain_file.name}")
        return None

    data = np.atleast_2d(data)

    data = data[data[:, 0] != 0]
    if len(data) == 0:
        print(f"Skipping {chain_file.parent.name}: no nonzero iteration rows in {chain_file.name}")
        return None

    iters = data[:, 0]
    vals  = data[:, 1]

    p = figure(
        title="Iteration vs Log Likelihood",
        x_axis_label="Iteration",
        y_axis_label="Log Likelihood",
        width=900,
        height=350,
        tools="pan,wheel_zoom,box_zoom,reset,save"
    )

    p.line(iters, vals, line_width=1, color="navy", alpha=0.8)
    p.grid.grid_line_alpha = 0.3
    return p


def extract_tic_sector(session_name: str):
    parts = session_name.split("_")

    tic_id = None
    sector = None

    for i, part in enumerate(parts):
        if part == "id" and i + 1 < len(parts):
            tic_id = parts[i + 1]
        if part == "sector" and i + 1 < len(parts):
            sector = parts[i + 1]

    return tic_id, sector

def make_original_lc_dropdown(tic_id, sector):
    original_plot = produce_original_tess_lc_plot(tic_id, sector)

    if original_plot is None:
        return None

    original_plot.visible = False

    button = Button(
        label="▶ Original TESS light curve",
        width=900,
        button_type="default"
    )

    button.js_on_click(CustomJS(
        args=dict(plot=original_plot, button=button),
        code="""
        plot.visible = !plot.visible;

        if (plot.visible) {
            button.label = "▼ Original TESS light curve";
        } else {
            button.label = "▶ Original TESS light curve";
        }
        """
    ))

    return column(button, original_plot)

def produce_original_tess_lc_plot(tic_id, sector, max_points=5000):
    try:
        tic_id = int(tic_id)
        sector = int(sector)

        print(f"Downloading original TESS LC for TIC {tic_id}, sector {sector}", flush=True)

        search = lk.search_lightcurve(
            f"TIC {tic_id}",
            mission="TESS",
            sector=sector
        )

        if len(search) == 0:
            return Div(
                text=f"<p><b>No original TESS light curve found for TIC {tic_id}, sector {sector}.</b></p>",
                width=900
            )

        lc = search[0].download()

        if lc is None:
            return Div(
                text=f"<p><b>Download failed for TIC {tic_id}, sector {sector}.</b></p>",
                width=900
            )

        lc = lc.remove_nans()

        time = np.asarray(lc.time.value, dtype=float)
        flux = np.asarray(lc.flux.value, dtype=float)

        good = np.isfinite(time) & np.isfinite(flux)
        time = time[good]
        flux = flux[good]

        if len(time) == 0:
            return Div(
                text=f"<p><b>No finite original light curve data for TIC {tic_id}, sector {sector}.</b></p>",
                width=900
            )

        if max_points is not None and len(time) > max_points:
            step = int(np.ceil(len(time) / max_points))
            time = time[::step]
            flux = flux[::step]

        p = figure(
            title=f"Original TESS Light Curve: TIC {tic_id}, Sector {sector}",
            x_axis_label="Time",
            y_axis_label="Flux",
            width=900,
            height=350,
            tools="pan,wheel_zoom,box_zoom,reset,save"
        )

        p.scatter(time, flux, size=3, alpha=0.7)
        p.grid.grid_line_alpha = 0.3

        return p

    except Exception as e:
        return Div(
            text=f"<p><b>Could not load original TESS LC for TIC {tic_id}, sector {sector}:</b> {e}</p>",
            width=900
        )

def make_one_session_plot(session_dir: Path):
    lc_file = session_dir / "folded_observed_and_model_light_curves.txt"
    chain_file = session_dir / "states.txt"

    if not lc_file.exists():
        print(f"Skipping {session_dir.name}: missing folded_observed_and_model_light_curves.txt")
        return None

    if not chain_file.exists():
        print(f"Skipping {session_dir.name}: missing states.txt")
        return None

    p_chain = produce_chains_plot(chain_file)
    p_flux = produce_light_curve_from_file(lc_file, session_dir)
    if p_chain is None or p_flux is None:
        return None

    p_resid = produce_residual_plot_from_file(lc_file, p_flux)
    tic_id, sector = extract_tic_sector(session_dir.name)

    header = Div(
        text=f"""
        <hr>
        <h1>TIC {tic_id} Sector {sector}</h1>
        <h3>Folded Observed vs Model Light Curve + MCMC Chain</h3>
        <p>{session_dir.name}</p>
        """,
        width=900
    )

    param_table = make_param_table(session_dir)
    # original_lc_dropdown = make_original_lc_dropdown(tic_id, sector)
    #
    # return column(
    #     header,
    #     param_table,
    #     original_lc_dropdown,
    #     p_chain,
    #     p_flux,
    #     p_resid
    # )
    return column(
        header,
        param_table,
        p_chain,
        p_flux,
        p_resid
    )


def make_nav_div(page_num, total_pages):
    links = ""

    if page_num > 1:
        prev_file = f"all_sessions_lc_and_chain_page_{page_num - 1}.html"
        links += f"""
        <a href="{prev_file}" style="
            display: inline-block;
            padding: 10px 15px;
            border: 1px solid black;
            margin-right: 10px;
            text-decoration: none;
            color: black;
        ">Previous</a>
        """

    if page_num < total_pages:
        next_file = f"all_sessions_lc_and_chain_page_{page_num + 1}.html"
        links += f"""
        <a href="{next_file}" style="
            display: inline-block;
            padding: 10px 15px;
            border: 1px solid black;
            margin-right: 10px;
            text-decoration: none;
            color: black;
        ">Next</a>
        """

    return Div(
        text=f"""
        <div style="margin: 20px 0;">
            <h3>Page {page_num} of {total_pages}</h3>
            {links}
        </div>
        """,
        width=900
    )


def read_bad_models_txt(txt_path):
    targets = set()

    with open(txt_path, "r") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            parts = line.replace(",", " ").split()

            if len(parts) < 2:
                continue

            tic_id = str(int(parts[0]))
            sector = str(int(parts[1]))

            targets.add((tic_id, sector))

    return targets

def main():
    session_dirs = sorted([
        p for p in SESSIONS_DIR.iterdir()
        if p.is_dir() and "2026_06_07" in p.name
    ])
    # session_dirs = sorted([
    #     p for p in SESSIONS_DIR.iterdir()
    #     if p.is_dir()
    # ])[-300:]

    total_sessions = len(session_dirs)
    total_pages = (total_sessions + PAGE_SIZE - 1) // PAGE_SIZE

    print(f"Found {total_sessions} sessions")
    print(f"Making {total_pages} pages with {PAGE_SIZE} light curves per page")

    for page_index in range(total_pages):
        page_num = page_index + 1

        start = page_index * PAGE_SIZE
        end = start + PAGE_SIZE
        page_session_dirs = session_dirs[start:end]

        sections = []

        title_div = Div(
            text=f"""
            <h1>All Session Light Curves and MCMC Chains</h1>
            <h2>Page {page_num} of {total_pages}</h2>
            <p>Showing sessions {start + 1} to {min(end, total_sessions)} out of {total_sessions}</p>
            """,
            width=900
        )

        sections.append(title_div)
        sections.append(make_nav_div(page_num, total_pages))

        for session_dir in page_session_dirs:
            section = make_one_session_plot(session_dir)
            if section is not None:
                sections.append(section)

        sections.append(make_nav_div(page_num, total_pages))

        layout = column(*sections)

        out_file = f"all_sessions_lc_and_chain_page_{page_num}_old.html"
        output_file(out_file)
        save(layout)

        print(f"Saved {out_file}")

if __name__ == "__main__":
    main()