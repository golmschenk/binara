import numpy as np
import matplotlib.pyplot as plt

# taking the original TESS light curve, folding it, and then binning it using
# 150 bins and using Siddhant's method for the more accurate error bars

input_file  = "data/lightcurves/original_folded_lightcurves/110602878_sector_34.txt"
output_file = "data/lightcurves/original_folded_lightcurves/110602878_sector_34_binned150.txt"

nbins = 150

def bin_error_from_linear_residuals(x, y):
    a, b = np.polyfit(x, y, 1)
    r = y - (a * x + b)
    return np.std(r, ddof=1)

with open(input_file, "r") as f:
    Nt, P_days = f.readline().split()
    Nt = int(Nt)
    P_days = float(P_days)

data = np.loadtxt(input_file, skiprows=1)
time = data[:, 0]
flux = data[:, 1]
err  = data[:, 2]

median_flux = np.median(flux)
flux_norm = flux / median_flux

t0 = time.min()
phase = ((time - t0) / P_days) % 1.0

bin_edges = np.linspace(0.0, 1.0, nbins + 1)
bin_centers_phase = 0.5 * (bin_edges[:-1] + bin_edges[1:])
binned_phase_days = bin_centers_phase * P_days

binned_flux_mean = np.full(nbins, np.nan)
binned_flux_err  = np.full(nbins, np.nan)

for i in range(nbins):
    in_bin = (phase >= bin_edges[i]) & (phase < bin_edges[i+1])
    idx = np.where(in_bin)[0]
    n = len(idx)

    if n == 0:
        continue

    fl_bin = flux_norm[idx]
    ph_bin = phase[idx]

    binned_flux_mean[i] = np.mean(fl_bin)

    if n >= 3:
        binned_flux_err[i] = bin_error_from_linear_residuals(ph_bin, fl_bin)
    elif n == 2:
        binned_flux_err[i] = np.std(fl_bin, ddof=1)
    else:
        binned_flux_err[i] = err[idx][0] / median_flux

with open(output_file, "w") as f:
    f.write(f"{nbins} {P_days:.10f}\n")
    for t_bin, f_mean, f_err in zip(binned_phase_days, binned_flux_mean, binned_flux_err):
        f.write(f"{t_bin:.10f} {f_mean:.10f} {f_err:.10f}\n")

print(f"Wrote folded + binned LC to: {output_file}")

err_file = "data/lightcurves/original_folded_lightcurves/110602878_sector_34_binned150_err.txt"
in_file = "data/lightcurves/original_folded_lightcurves/110602878_sector_34_binned150.txt"
data = np.loadtxt(in_file, skiprows=1)
old_err= data[:, 2]
with open(err_file, "w") as f:
    f.write("original lc error bars  model lc error bars\n")
    for t, e in zip(old_err, binned_flux_err):
        f.write(f"{t:.10f} {e:.10e}\n")

plt.figure(figsize=(10, 5))

plt.scatter(
    phase * P_days,
    flux_norm,
    s=6,
    color="gray",
    alpha=0.35,
    label="Folded data (raw)"
)

plt.errorbar(
    binned_phase_days,
    binned_flux_mean,
    yerr=binned_flux_err,
    fmt="o",
    color="firebrick",
    ecolor="firebrick",
    elinewidth=1,
    capsize=2,
    markersize=4,
    label="Binned (150 bins)"
)

plt.xlabel("Phase (days)")
plt.ylabel("Normalized Flux")
plt.title("Folded & Binned Light Curve with Error Bars")
plt.legend()
plt.tight_layout()
plt.show()
