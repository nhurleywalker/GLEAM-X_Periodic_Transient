import numpy as np
import matplotlib.pyplot as plt

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

plt.figure("P, Pdot search")
plt.imshow(SNRs, extent=(Pdots[0]-Pdotstep/2, Pdots[-1]+Pdotstep/2, Ps[0]-Pstep/2, Ps[-1]+Pstep/2), origin='lower', aspect='auto')
plt.colorbar()
plt.xlabel("Pdot (s/s)")
plt.ylabel("P (s)")

plt.show()
