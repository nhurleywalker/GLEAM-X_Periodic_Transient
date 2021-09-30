#!/usr/bin/env python
import numpy as np

from astropy.io import fits
from astropy.time import Time
from astropy.visualization import PercentileInterval

from glob import glob
import os
import sys

from matplotlib import pyplot as plt
#import matplotlib.pylab as pl
from matplotlib.ticker import FormatStrFormatter
from matplotlib import gridspec
import matplotlib.colors as mcol
import matplotlib.cm as cm

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
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

# This is the best value of the period for Pdot = 0, which results in a stack that lines up good enough for the figure
P = 1091.1708 # s
# Obsids where you can't really see the pulse in the stack
excludes = [1205955280, 1205955160, 1204998416, 1204837904, 1204831488, 1204225936, 1204225816]

metas = sorted(glob("../metafits/*.metafits"), reverse=False)

for ex in excludes:
    print(f"Excluding {ex}")
    metas.remove(f"../metafits/{ex}.metafits")

fig = plt.figure(figsize=(5,10))
# One for each of the profiles
gs = gridspec.GridSpec(int(len(metas)/2), 2)
ax = None
month = None
j = 0
m = 0
for i in np.arange(0, len(metas)):
    meta = metas[i]
    obsid = os.path.basename(meta)[0:10]
    ti = Time(obsid, format="gps")
    ti.format="ymdhms"
    t = ti.value
    cmonth = wordMonth(t[1])
    hdu = fits.open(meta)
    freqcent = hdu[0].header["FREQCENT"]
    bandwidth = hdu[0].header["BANDWDTH"]
    hdu.close()
    dat = np.loadtxt("../dedispersed_profiles/{0}_I_57.0_flagged_profile.dat".format(obsid))
    dat = dat.T
    # Figure out which column we're indexing in gs
    if j == 0:
        m = i
    else:
        m = i - 32
    # If it's not the first time you've created a profile plot, share the x axis
    if ax:
        ax = plt.subplot(gs[m, j], sharex = ax)
    # If it's the first time you've created a profile plot, make a new axis
    else:
        ax = plt.subplot(gs[m, j])
    x = np.mod(dat[0] + float(obsid) + P/2, P) - 0.48*P
    ind = np.where(x!=0)
    line1, = ax.plot(x[ind], dat[1][ind], lw=1, color=cm1((freqcent - 88.)/(215.-88.)), zorder=1)
    ax.axes.get_xaxis().set_visible(False)
#    ax.set_frame_on(False)
    ax.set_yticks([])
    # So each month is only printed once
    if cmonth != month:
        ylabel = f"{cmonth} {t[2]:02d} {t[3]:02d}:{t[4]:02d}"
        ax.yaxis.set_label_coords(0,0.3)
        month = cmonth
    else:
        ylabel = f"{t[2]:02d} {t[3]:02d}:{t[4]:02d}"
        ax.yaxis.set_label_coords(0.075,0.3)
    ax.set_ylabel(ylabel, rotation=0, zorder = 4)#, loc="bottom")
    ax.set_xlim([-100,100])
    # Sans-serif fonts
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    # Remove unnecessary lines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    # Reduce padding
    ax.margins(0.0)
    # Sans-serif fonts for Nature
    ax.xaxis.get_major_formatter()._usetex = False
    ax.yaxis.get_major_formatter()._usetex = False
    if i == int(len(metas)/2) - 1 or i == int(len(metas)) - 1:
    # Increment the index for the gridspec so it starts a new column
        j += 1
        ax.spines["bottom"].set_visible(True)
        ax.axes.get_xaxis().set_visible(True)
        ax.set_xlabel("Time / seconds")

plt.subplots_adjust(hspace=.0)
fig.savefig("pulse_stack.pdf", bbox_inches="tight")
fig.savefig("pulse_stack.png", bbox_inches="tight", dpi=300)
