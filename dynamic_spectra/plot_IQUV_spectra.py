import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.font_manager
from matplotlib import rc
# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

try:
    obsid = sys.argv[1]
except IndexError:
    obsid = 1199591976
DM = 57.0
stokes = ["I", "Q", "U", "V"]
nstokes = len(stokes)
dt = 0.5 # seconds per pixel
t0 = 50

freq_lo = 170.24
BW = 30.72

# Minimum and maximum time range
xmin = -30.25
xmax = 20.75

# Minimum and maximum image range
vmin = -20.
vmax = 25.

def get_dat_filename(obsid, stokes):
    return "{}_{}_{:.1f}_dedispersed_dynspec.dat".format(obsid, stokes, DM)

def plot_dynspec(obsid, stokes, ax):
    filename = get_dat_filename(obsid, stokes)
    dat = np.loadtxt(filename).T
    t = np.arange(dat.shape[1])*dt - t0
    extent = np.array([t[0] - 0.5*dt, t[-1] + 0.5*dt, freq_lo, freq_lo + BW])
    image = ax.imshow(dat, cmap="plasma", origin='lower', interpolation='none', aspect='auto', extent=extent, vmin=vmin, vmax=vmax)
    return image, dat, t

if __name__ == "__main__":
    fig, axs = plt.subplots(nrows=2, ncols=nstokes, sharex=True)
    for i in range(nstokes):
        print(i)
        ax = axs[-1][i]
        S = stokes[i]
        image, dat, t = plot_dynspec(obsid, S, ax)
        ax.set_xlim([xmin, xmax])
        ax.set_xlabel("Time (s)")
        ax.text(xmin+4, 198., S, backgroundcolor="white")
            # Sans-serif fonts for Nature
        ax.xaxis.get_major_formatter()._usetex = False
        ax.yaxis.get_major_formatter()._usetex = False
        if i > 0:
            ax.set_yticks([])
#        else:
#            axs[0][0].plot(t, profile)
#            axs[0][0].set_ylabel("$S$ (Jy\,beam$^{-1}$)")
        if i == 0:
            profile = np.nanmean(dat, axis=0)
            ax.set_ylabel("Frequency (MHz)")
        fig.delaxes(axs[0][i])
    miniax = fig.add_subplot(445, sharex=axs[0][1])
    miniax.set_ylabel("Frequency (MHz)")
    miniax.set_ylim([-3, 28])
    # Only this axis; https://stackoverflow.com/questions/4209467/matplotlib-share-x-axis-but-dont-show-x-axis-tick-labels-for-both-just-one
    plt.setp(miniax.get_xticklabels(), visible=False)
    miniax.plot(t, profile, color="black", lw=0.5)
    miniax.set_ylabel("$S$ (Jy\,beam$^{-1}$)")
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.savefig("{0}_dynspec.pdf".format(obsid), bbox_inches="tight")
    fig.savefig("{0}_dynspec.png".format(obsid), bbox_inches="tight", dpi=1000)
