import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib.ticker import ScalarFormatter

# Nature requires sans-serif fonts
plt.rcParams.update({
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})
#    "font.sans-serif": ["Arial"]})
#    "font.family": "sans-serif",

cm = 1/2.54  # centimeters in inches
#    "text.usetex": False,

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


fig = plt.figure(figsize=(8*cm,8*cm))
ax = fig.add_subplot(111)
#plt.figure("P, Pdot search")
#plt.imshow(SNRs, extent=(Pdots[0]-Pdotstep/2, Pdots[-1]+Pdotstep/2, Ps[0]-Pstep/2, Ps[-1]+Pstep/2), origin='lower', aspect='auto')
#plt.colorbar()
ax.contour(data, extent=(Pdots[0]-Pdotstep/2, Pdots[-1]+Pdotstep/2, Ps[0]-Pstep/2, Ps[-1]+Pstep/2), origin='lower', aspect='auto', levels=[10,11,12,13,14,15], cmap="Blues", linewidths=0.5)
ax.axvline(0, ls=":", color="black", lw=0.5)
# Values from p-pdot_search.dat
ax.scatter(6.000000000000029e-10, 1091.1688000000026, marker="+", color="darkred", s=300, zorder=100, linewidth=0.5)
#ax.colorbar()
ax.set_xlabel(r"$\dot{P}$ / s s$^{-1}$")
ax.set_ylabel(r"$P$ /s")
ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True, useOffset=False))
#ax.ticklabel_format(axis="y", style="plain", useOffset=True)
fig.savefig("P_Pdot.pdf", bbox_inches="tight", dpi=300)
fig.savefig("P_Pdot.png", bbox_inches="tight", dpi=300)
fig.savefig("P_Pdot.eps", bbox_inches="tight", dpi=300)


