import sys
import numpy as np
from astropy.io import fits
from astropy.time import Time
from matplotlib import pyplot as plt
import matplotlib.pylab as pl
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
from astropy.visualization import PercentileInterval
from glob import glob
import matplotlib.colors as mcol

def wordMonth(intMonth):
    if intMonth == 3:
        cmonth = "Mar"
    elif intMonth == 2:
        cmonth = "Feb"
    elif intMonth == 1:
        cmonth = "Jan"
    else:
        cmonth = None
    return cmonth

# I want to try red to blue for the frequencies
cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["c","m"])

# Sans-serif fonts for Nature
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54  # centimeters in inches

# This is the best value of the period for Pdot = 0, which results in a stack that lines up good enough for the figure
P = 1091.1708 # s

metas = sorted(glob("*.metafits"), reverse=False)

fig = plt.figure(figsize=(10*cm,8*cm))
# One for each of the profiles, [and one for the colorbar -- if it were working]
#gs = gridspec.GridSpec(len(metas)+1, 1)
gs = gridspec.GridSpec(len(metas), 1)
#c = np.arange(1, 6) # we have five frequencies
# Make dummy mappable so we can add the colorbar later
#ax = fig.add_subplot(111)
#dummie_cax = ax.scatter(c, c, c=c, cmap=cm1)
# Clear axis
#ax.cla()
ax = None
month = None
for i in np.arange(0, len(metas)):
    meta = metas[i]
    obsid = meta[0:10]
    ti = Time(obsid, format="gps")
    ti.format="ymdhms"
    t = ti.value
    cmonth = wordMonth(t[1])
    hdu = fits.open(meta)
    freqcent = hdu[0].header["FREQCENT"]
    bandwidth = hdu[0].header["BANDWDTH"]
    hdu.close()
    dat = np.loadtxt("{0}_I_57.0_flagged_profile.dat".format(obsid))
    dat = dat.T
    if ax:
        ax = plt.subplot(gs[i], sharex = ax)
    else:
        ax = plt.subplot(gs[i])
    x = np.mod(dat[0] + float(obsid) + P/2, P) - 0.48*P
    ind = np.where(x!=0)
#    line1, = ax.plot(x, dat[1], lw=1, color=pl.cm.jet((freqcent - 88.)/(215.-88.)))
    line1, = ax.plot(x[ind], dat[1][ind], lw=0.75, color=cm1((freqcent - 88.)/(215.-88.)))
    ax.axes.get_xaxis().set_visible(False)
#    ax.set_frame_on(False)
    ax.set_yticks([])
    # Means each month is only printed once
    if cmonth != month:
        ylabel = f"{cmonth} {t[2]:02d} {t[3]:02d}:{t[4]:02d}"
        ax.yaxis.set_label_coords(0.13,0.2)
        month = cmonth
    else:
        ylabel = f"{t[2]:02d} {t[3]:02d}:{t[4]:02d}"
        ax.yaxis.set_label_coords(0.17,0.2)
    ax.set_ylabel(ylabel, rotation=0)#, loc="bottom")
    ax.set_xlim([-100,100])
    ax.axvline(-30, ls="--", lw=1, color="black")
    ax.axvline(26, ls="--", lw=1, color="black")
ax.axes.get_xaxis().set_visible(True)
ax.set_xlabel("Time (seconds)")
ax.xaxis.get_major_formatter()._usetex = False
ax.yaxis.get_major_formatter()._usetex = False
plt.subplots_adjust(hspace=.0)
fig.savefig("subpulse_comparison.png", bbox_inches="tight", dpi=300)
fig.savefig("subpulse_comparison.pdf", bbox_inches="tight", dpi=300)
fig.savefig("subpulse_comparison.eps", bbox_inches="tight", dpi=300)
