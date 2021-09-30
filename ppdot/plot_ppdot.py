import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage

datafile = "p-pdot_search.dat"

P0    = 1091.1558
Pf    = 1091.1858
Pstep = 0.0005
Ps    = np.arange(P0, Pf, Pstep)

Pdot0    = -4e-9
Pdotf    = 4e-9
Pdotstep = 0.1e-9
Pdots    = np.arange(Pdot0, Pdotf, Pdotstep)

SNRs = np.loadtxt(datafile)
# For smooth contours
data = ndimage.zoom(SNRs, 3)


fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(111)
#plt.figure("P, Pdot search")
#plt.imshow(SNRs, extent=(Pdots[0]-Pdotstep/2, Pdots[-1]+Pdotstep/2, Ps[0]-Pstep/2, Ps[-1]+Pstep/2), origin='lower', aspect='auto')
#plt.colorbar()
ax.contour(data, extent=(Pdots[0]-Pdotstep/2, Pdots[-1]+Pdotstep/2, Ps[0]-Pstep/2, Ps[-1]+Pstep/2), origin='lower', aspect='auto', levels=[10,11,12,13,14,15], cmap="Blues")
ax.axvline(0, ls=":", color="black")
# Values from p-pdot_search.dat
ax.scatter(6.000000000000029e-10, 1091.1688000000026, marker="+", color="darkred", s=300, zorder=100)
#ax.colorbar()
ax.set_xlabel(r"$\dot{P}$ / s s$^{-1}$")
ax.set_ylabel(r"$P$ /s")
fig.savefig("P_Pdot.pdf", bbox_inches="tight")
fig.savefig("P_Pdot.png", bbox_inches="tight")


