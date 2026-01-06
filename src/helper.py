import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

show_error_bars = True
df = pd.read_csv("model_folded_limbdarkening.csv", header=None, skiprows=1,
                 names=["time_days", "flux_model", "flux_other"])

df = pd.read_csv("model_folded.csv", header=None, skiprows=1,
                names=["time_days", "flux_model", "flux_other"])

time = df["time_days"]
flux_model = df["flux_model"]
flux_other = df["flux_other"]

# plt.figure(figsize=(8,4))
# plt.plot(time, flux_model, '.', ms=2)
# plt.xlabel("Time (days)")
# plt.ylabel("Flux Model")
# plt.title("Time vs Model Flux")
# plt.tight_layout()
# plt.show()
err_file = "data/lightcurves/original_folded_lightcurves/110602878_sector_34_binned150_err.txt"
err_data = np.loadtxt(err_file, skiprows=1)
other_flux_err = err_data[:, 0]
model_flux_err = err_data[:, 1]

plt.figure(figsize=(8,4))
if show_error_bars:
    plt.errorbar(
        time,
        flux_other,
        yerr=other_flux_err,
        fmt='.',
        ms=2,
        alpha=0.6,
        label="Observed Flux"
    )
else:
    plt.plot(
        time,
        flux_other,
        '.',
        ms=2,
        alpha=0.6,
        label="Observed Flux"
    )

plt.plot(
    time,
    flux_model,
    '.',
    ms=2,
    label="Model Flux"
)

plt.xlabel("Time")
plt.ylabel("Flux")
plt.legend(['Model Flux', 'Observed Flux'])
plt.title("Time vs Flux")
plt.tight_layout()

plt.show()