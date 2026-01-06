import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("data/chains/chain.110602878_sector_34_gmag_OMP_1.dat")
# data = np.loadtxt("data/chains/chain.110602878_3mil_nolimbdarkening.dat")
#
# data = np.loadtxt("data/chains/chain.110602878_3mil_limbdarkening.dat")

data = data[data[:, 0] != 0]
data = data[data[:, 1] < 1e6]
# data = data[data[:, 0] <= 175_000]
iters = data[:, 0]
vals  = data[:, 1]

plt.figure()
plt.plot(iters, vals)
plt.xlabel("Iteration")
plt.ylabel("Log Likelihood")
plt.title("Iteration vs. Log Likelihood")
plt.grid(True)
plt.tight_layout()
plt.show()