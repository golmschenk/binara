from pathlib import Path
from astropy.timeseries import BoxLeastSquares
import pandas as pd
import lightkurve as lk
from astroquery.mast import Observations
from astropy.io import fits
import numpy as np
import re
from binara.init_data import write_mcmc_data # correct one
#from init_data import write_mcmc_data
from astroquery.mast import Catalogs

reg = re.compile(r"phot_(\d+)-s(\d{4})", re.IGNORECASE)

from pathlib import Path
import subprocess
import shutil

def submit_binara_job(tic_id: int, sector: int):
    subprocess.run(
        ["sbatch", "job3_script.sh", str(tic_id), str(sector)],
        check=True
    )

    print(f"Submitted job for TIC {tic_id}, sector {sector}")

def load_and_normalize_one_sector(tic, sector):
    search = lk.search_lightcurve(f"TIC {tic}", mission="TESS", sector=sector)

    if len(search) == 0:
        print(f"TIC {tic}, sector {sector}: no light curve found")
        return None

    lc = search[0].download()

    if lc is None:
        print(f"TIC {tic}, sector {sector}: download failed")
        return None

    lc = lc.remove_nans()

    t = np.asarray(lc.time.value, dtype=np.float64)
    y = np.asarray(lc.flux.value, dtype=np.float64)
    yerr = np.asarray(lc.flux_err.value, dtype=np.float64)

    good = np.isfinite(t) & np.isfinite(y) & np.isfinite(yerr) & (yerr > 0)
    t, y, yerr = t[good], y[good], yerr[good]

    if len(t) == 0:
        return None

    med = np.nanmedian(y)
    if not np.isfinite(med) or med == 0:
        return None

    y_norm = y / med
    yerr_norm = yerr / abs(med)
    return np.column_stack([t, y_norm, yerr_norm])

def extract_tic_sector(path: str):
    m = reg.search(path)
    if not m:
        return None
    tic_id = int(m.group(1))
    sector = int(m.group(2))
    return tic_id, sector

def write_period_to_csv(tic_id, sector, best_period):
    best_path = Path("input_data/best_periods.csv")
    columns = ["TIC ID", "Sector", "Best Period (days)"]

    if best_path.exists():
        best_df = pd.read_csv(best_path)
    else:
        best_df = pd.DataFrame(columns=columns)

    mask = (
            (best_df["TIC ID"].astype(float).astype(int) == int(tic_id)) &
            (best_df["Sector"].astype(float).astype(int) == int(sector))
    )

    if mask.any():
        first_idx = best_df[mask].index[0]

        best_df.loc[first_idx, "TIC ID"] = int(tic_id)
        best_df.loc[first_idx, "Sector"] = int(sector)
        best_df.loc[first_idx, "Best Period (days)"] = round(float(best_period), 6)

        duplicate_indices = best_df[mask].index[1:]
        best_df = best_df.drop(index=duplicate_indices)

    else:
        new_row = pd.DataFrame([{
            "TIC ID": int(tic_id),
            "Sector": int(sector),
            "Best Period (days)": round(float(best_period), 6)
        }])

        best_df = pd.concat([best_df, new_row], ignore_index=True)

    best_df.to_csv(best_path, index=False)

def find_best_period_from_array(arr, tic_id):
    t = arr[:, 0]
    y = arr[:, 1]
    yerr = arr[:, 2]

    bls = BoxLeastSquares(t, y, dy=yerr)

    periods = np.linspace(1, 30, 10000)
    bls_results = bls.power(periods, 0.1)

    best_idx = np.nanargmax(bls_results.power)
    best_period = float(bls_results.period[best_idx])

    print(f"TIC {tic_id}: best period = {best_period:.6f} days")

    return best_period

def find_best_period_all_available_sectors(tic, sector_min, sector_max):
    arrays = []
    sectors_used = []

    for s in range(sector_min, sector_max + 1):
        arr = load_and_normalize_one_sector(tic, s)

        if arr is not None:
            arrays.append(arr)
            sectors_used.append(s)
            print(f"Using sector {s}, points = {len(arr)}")

    if len(arrays) == 0:
        raise ValueError(f"No usable sectors found for TIC {tic}")
    big_arr = np.vstack(arrays)
    big_arr = big_arr[np.argsort(big_arr[:, 0])]

    print(f"\nCombined sectors for TIC {tic}: {sectors_used}")
    print(f"Total points: {len(big_arr)}")

    best_period = find_best_period_from_array(big_arr, tic)
    return best_period, sectors_used, big_arr

def get_period_from_csv(tic_id, sector, csv_path="input_data/best_periods.csv"):
    df = pd.read_csv(csv_path)

    row = df[
        (df["TIC ID"].astype(int) == int(tic_id)) &
        (df["Sector"].astype(int) == int(sector))
        ]

    if row.empty:
        raise ValueError(f"No period found for TIC {tic_id}, sector {sector}")

    return float(row["Best Period (days)"].iloc[0])

def create_input_directory(tic_id, sector_num, parent_dir):
    parent_dir = Path(parent_dir)
    folder_name = f"tic_id_{tic_id}_sector_{sector_num}"
    folder_path = parent_dir / folder_name
    folder_path.mkdir(parents=True, exist_ok=True)
    file_names = ["folded_observed_light_curve.txt", "magnitudes_and_colors.txt", "py_initialize.txt"]
    for file_name in file_names:
        file_path = folder_path / file_name
        file_path.touch(exist_ok=True)

    return folder_path

# def save_tess_lc_from_tic_sector(tic_id: int, sector: int):
#     obs = Observations.query_criteria(
#         provenance_name="TESS-SPOC",
#         target_name=str(tic_id),
#         sequence_number=sector
#     )
#
#     if len(obs) == 0:
#         raise FileNotFoundError(
#             f"No TESS-SPOC product found for TIC {tic_id}, sector {sector}"
#         )
#
#     products = Observations.get_product_list(obs)
#
#     mask = (
#             (products["productType"] == "SCIENCE") &
#             np.char.endswith(products["productFilename"].astype(str), "lc.fits")
#     )
#     lc_products = products[mask]
#
#     if len(lc_products) == 0:
#         raise FileNotFoundError(
#             f"No lc.fits file found for TIC {tic_id}, sector {sector}"
#         )
#
#     manifest = Observations.download_products(lc_products[:1], mrp_only=False)
#     local_path = manifest["Local Path"][0]
#
#     with fits.open(local_path) as hdul:
#         data = hdul[1].data
#         time = data["TIME"]
#         flux = data["PDCSAP_FLUX"]
#         flux_err = data["PDCSAP_FLUX_ERR"]
#
#     arr = np.column_stack((time, flux, flux_err))
#     arr = arr[np.isfinite(arr).all(axis=1)]
#
#     Path("input_data/raw_data").mkdir(parents=True, exist_ok=True)
#     output_file = f"input_data/raw_data/tic_id_{tic_id}_sector_{sector}_raw_lightcurve.txt"
#     period = get_period_from_csv(tic_id, sector)
#     with open(output_file, "w") as f:
#         f.write(f"{len(arr)+2} {period}\n")
#         np.savetxt(f, arr, fmt="%.10f")

# def get_tess_lc_array(tic_id: int, sector: int):
#     obs = Observations.query_criteria(
#         provenance_name="TESS-SPOC",
#         target_name=str(tic_id),
#         sequence_number=sector
#     )
#
#     if len(obs) == 0:
#         raise FileNotFoundError(
#             f"No TESS-SPOC product found for TIC {tic_id}, sector {sector}"
#         )
#
#     products = Observations.get_product_list(obs)
#
#     mask = (
#             (products["productType"] == "SCIENCE") &
#             np.char.endswith(products["productFilename"].astype(str), "lc.fits")
#     )
#
#     lc_products = products[mask]
#
#     if len(lc_products) == 0:
#         raise FileNotFoundError(
#             f"No lc.fits file found for TIC {tic_id}, sector {sector}"
#         )
#
#     manifest = Observations.download_products(lc_products[:1], mrp_only=False)
#     local_path = manifest["Local Path"][0]
#
#     with fits.open(local_path) as hdul:
#         data = hdul[1].data
#         time = data["TIME"]
#         flux = data["PDCSAP_FLUX"]
#         flux_err = data["PDCSAP_FLUX_ERR"]
#
#     arr = np.column_stack((time, flux, flux_err))
#     arr = arr[np.isfinite(arr).all(axis=1)]
#
#     return arr
def bin_error_from_linear_residuals(x, y):
    if len(x) < 3:
        return np.std(y, ddof=1) if len(y) > 1 else np.nan

    a, b = np.polyfit(x, y, 1)
    residuals = y - (a * x + b)
    return np.std(residuals, ddof=1)

# def fold_and_bin2(input_file: str, output_file: str, nbins: int = 150):
#     input_file = Path(input_file)
#     output_file = Path(output_file)
#
#     with open(input_file, "r") as f:
#         header = f.readline().strip()
#
#     parts = header.split()
#     if len(parts) < 2:
#         raise ValueError(f"Bad header: {header}")
#
#     P_days = float(parts[1])
#
#     data = np.loadtxt(input_file, skiprows=1)
#     time = data[:, 0]
#     flux = data[:, 1]
#     err = data[:, 2]
#
#     median_flux = np.median(flux)
#     flux_norm = flux / median_flux
#
#     t0 = time.min()
#     phase = ((time - t0) / P_days) % 1.0
#
#     bin_edges = np.linspace(0.0, 1.0, nbins + 1)
#     bin_centers_phase = 0.5 * (bin_edges[:-1] + bin_edges[1:])
#     binned_phase_days = bin_centers_phase * P_days
#
#     binned_flux_mean = np.full(nbins, np.nan)
#     binned_flux_err  = np.full(nbins, np.nan)
#
#     for i in range(nbins):
#         in_bin = (phase >= bin_edges[i]) & (phase < bin_edges[i+1])
#         idx = np.where(in_bin)[0]
#         n = len(idx)
#
#         if n == 0:
#             continue
#
#         fl_bin = flux_norm[idx]
#         ph_bin = phase[idx]
#
#         binned_flux_mean[i] = np.mean(fl_bin)
#
#         if n >= 3:
#             binned_flux_err[i] = bin_error_from_linear_residuals(ph_bin, fl_bin)
#         elif n == 2:
#             binned_flux_err[i] = np.std(fl_bin, ddof=1)
#         else:
#             binned_flux_err[i] = err[idx][0] / median_flux
#
#     with open(output_file, "w") as f:
#         f.write(f"{nbins} {P_days:.10f}\n")
#         for t_bin, f_mean, f_err in zip(
#                 binned_phase_days,
#                 binned_flux_mean,
#                 binned_flux_err
#         ):
#             if np.isnan(f_mean):
#                 continue
#             f.write(f"{t_bin:.10f} {f_mean:.10f} {f_err:.10f}\n")

def fold_and_write_combined_light_curve(output_file, combined_arr, period, nbins):
    output_file = Path(output_file)

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

        # f_err = f_err / 10.0
        if np.isfinite(f_mean) and np.isfinite(f_err):
            rows.append((binned_phase_days[i], f_mean, f_err))

    with open(output_file, "w") as f:
        f.write(f"{len(rows)} {period:.10f}\n")

        for t_bin, f_mean, f_err in rows:
            f.write(f"{t_bin:.10f} {f_mean:.10f} {f_err:.10f}\n")

    print(f"Wrote {len(rows)} binned points to {output_file}")

def read_lightcurves_from_txt(txt_path="repeat_heartbeats.txt"):
    lightcurves = []

    with open(txt_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            filename = line.replace("\\", "/").split("/")[-1]
            if filename.endswith(".fits"):
                lightcurves.append(filename)

    return lightcurves

def main():
    #df = pd.read_csv("input_data/best_periods.csv")
    input_root = Path("input_data")

    lightcurves = read_lightcurves_from_txt("src/repeat_heartbeats.txt")
    #lightcurves = [
    #    f"hlsp_tess-spoc_tess_phot_{int(tic_id):016d}-s{int(sector):04d}_tess_v1_lc.fits"
    #    for tic_id, sector in zip(df["TIC ID"], df["Sector"])
    #]
    #lightcurves = ["hlsp_tess-spoc_tess_phot_0000000302795956-s0037_tess_v1_lc.fits"]
    seen_tics = set()
    for lc in lightcurves:
        # extract the tic id and sector number for lc
        tic_id, sector = extract_tic_sector(lc)
        print(f"Working on TIC {tic_id}, sector {sector}", flush=True)
        print(f"Working on TIC {tic_id}, sector {sector}")
        if tic_id in seen_tics:
            print(f"Skipping repeated TIC {tic_id}, sector {sector}", flush=True)
            continue
        # create an input data folder for lc
        par_dir = "input_data"

        try:
            best_period, sectors_used, combined_arr = find_best_period_all_available_sectors(
                tic_id,
                sector_min=27,
                sector_max=55
            )

            if len(sectors_used) > 1:
                nbins = 300
            else:
                nbins = 150

            print(f"TIC {tic_id}: sectors used = {sectors_used}")
            print(f"TIC {tic_id}: using nbins = {nbins}")

            write_period_to_csv(tic_id, sector, best_period)
            create_input_directory(tic_id, sector, par_dir)

            folded_file = f"input_data/tic_id_{tic_id}_sector_{sector}/folded_observed_light_curve.txt"

            fold_and_write_combined_light_curve(
                output_file=folded_file,
                combined_arr=combined_arr,
                period=best_period,
                nbins=nbins
            )

            catalog_data = Catalogs.query_object(
                f"TIC {tic_id}",
                radius=0.0001,
                catalog="TIC"
            )

            try:
                row = catalog_data[0]
                distance = row["d"]
                gmag = row["GAIAmag"]
                gmag_err = row["e_GAIAmag"]

                input_mag_file = f"input_data/tic_id_{tic_id}_sector_{sector}/magnitudes_and_colors.txt"

                with open(input_mag_file, "w") as f:
                    f.write(f"{distance}\n")
                    f.write(f"{gmag} {gmag_err}\n")
                    f.write("-0.439\t1.123710371937538\n")
                    f.write("0.8318000000000002\t0.8003815310019341\n")
                    f.write("0.6298999999999993\t0.6600445403128677\n")
            except:
                continue

            write_mcmc_data(tic_id, sector)
            submit_binara_job(tic_id, sector)
            print(f"Submitted job for TIC {tic_id}, sector {sector}", flush=True)

        except Exception as e:
            print(f"Failed for TIC {tic_id}, sector {sector}: {e}")
            continue

if __name__ == '__main__':
    main()