import numpy as np
import matplotlib.pyplot as plt

datafile = "p-pdot_search_best_profile.dat"

profile = np.loadtxt(datafile)

plt.figure("P, Pdot search")
plt.plot(profile)
plt.xlabel("bin")
plt.ylabel("Mean flux")

plt.show()
