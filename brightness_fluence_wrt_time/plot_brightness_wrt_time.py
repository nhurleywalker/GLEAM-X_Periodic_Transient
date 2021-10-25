import matplotlib as mpl
mpl.use("Agg")
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

import numpy as np
from astropy.io import fits
from astropy.time import Time

from glob import glob

import os
import sys

import matplotlib.font_manager
from matplotlib import rc
# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

from astropy.time import Time

import matplotlib.patches as mpatches
from matplotlib import markers
from matplotlib.path import Path

cm = 1/2.54  # centimeters in inches

# https://stackoverflow.com/questions/26686722/align-matplotlib-scatter-marker-left-and-or-right
# To align the limits with their real values
def align_marker(marker, halign='center', valign='middle',):
    """
    create markers with specified alignment.

    Parameters
    ----------

    marker : a valid marker specification.
      See mpl.markers

    halign : string, float {'left', 'center', 'right'}
      Specifies the horizontal alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'center',
      -1 is 'right', 1 is 'left').

    valign : string, float {'top', 'middle', 'bottom'}
      Specifies the vertical alignment of the marker. *float* values
      specify the alignment in units of the markersize/2 (0 is 'middle',
      -1 is 'top', 1 is 'bottom').

    Returns
    -------

    marker_array : numpy.ndarray
      A Nx2 array that specifies the marker path relative to the
      plot target point at (0, 0).

    Notes
    -----
    The mark_array can be passed directly to ax.plot and ax.scatter, e.g.::

        ax.plot(1, 1, marker=align_marker('>', 'left'))

    """

    if isinstance(halign, str): #, unicode)):
        halign = {'right': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'left': 1.,
                  }[halign]

    if isinstance(valign, str): #, unicode)):
        valign = {'top': -1.,
                  'middle': 0.,
                  'center': 0.,
                  'bottom': 1.,
                  }[valign]

    # Define the base marker
    bm = markers.MarkerStyle(marker)

    # Get the marker path and apply the marker transform to get the
    # actual marker vertices (they should all be in a unit-square
    # centered at (0, 0))
    m_arr = bm.get_path().transformed(bm.get_transform()).vertices

    # Shift the marker vertices for the specified alignment.
    m_arr[:, 0] += halign / 2
    m_arr[:, 1] += valign / 2

    return Path(m_arr, bm.get_path().codes)

# Delete the radio pseudo-luminosity width-frequency text file before starting
if os.path.exists("luminosity_nuW.txt"):
    os.remove("luminosity_nuW.txt")

# Previously found spectral index to be -1.16
alpha = -1.16

# Previously found distance to be 1.3 kpc
dist = 1.3 #kpc

# Get the central channels and pulse widths of each obsid
chans = {}
widths = {}
with open("det_freqs_widths.txt") as f:
    for line in f:
        (key, chan, width) = line.split()
        chans[int(key)] = chan
        widths[int(key)] = width

limits = np.loadtxt("non-detections.txt")
limits = limits.T

files = sorted(glob("../dedispersed_profiles/1?????????_I_57.0_flagged_profile.dat"))

# Six pulses are so close to the edge of the observation that the pulse does not appear at all frequencies
# So the profile is flagged after dedispersion. We will also not include these in the plotting.
# One pulse is split between the final two observations (1205955160, 1205955280): treat this separately below
to_remove = [1199510192, 1204225816, 1204225936, 1204831488, 1204837904, 1205955160, 1205955280, 1205956360]

for f in to_remove:
    files.remove("../dedispersed_profiles/{0}_I_57.0_flagged_profile.dat".format(f))

fig = plt.figure(figsize=(18.3*cm,7*cm))
ax = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
peaks = []
fluences = []
times = []
trunc_peaks = []
trunc_fluences = []
trunc_times = []
for f in files:
    obsid = int(os.path.basename(f)[0:10])
    vals = np.loadtxt(f).T[1]
# Remove very negative values -- they are artefacts of the dedispersion
    vals[np.where(vals < -0.1)] = np.nan
    peak = np.nanmax(vals)
    fluence = np.nansum(vals*0.5)
    ind = np.squeeze(np.where(vals==peak))
# Central frequency is MWA channel number x 1.28
    nu = 1.28 * float(chans[obsid])
# Normalise to 154 MHz
    peak = peak * (154. / nu)**alpha
    fluence = fluence * (154. / (1.28 * float(chans[obsid])))**alpha
    t = float(obsid) + (float(ind) / 2.)
# If the peak is close to the start or end of the observation, then it is truncated and only a lower limit
# Note that these files all have 30 timesteps of padding at the beginning and end, so the constraint below
# is really within 15s (30 timesteps + padding) of the start or end of the observation
    if ind < 60 or ind > (len(vals) - 60):
        trunc_times.append(t)
        trunc_peaks.append(peak)
        trunc_fluences.append(fluence)
    else:
        times.append(t)
        peaks.append(peak)
        fluences.append(fluence)
    # For well-measured pulses, accumulate values for luminosity vs nuW plot
        with open("luminosity_nuW.txt", "a") as f:
            f.write("{0} {1}\n".format(peak*(dist**2), float(widths[obsid])*nu/1.e3))

# Read in the pair of observations which split a pulse
pair_peak = 0.0
pair_fluence = 0.0
files = ["../dedispersed_profiles/1205955160_I_57.0_flagged_profile.dat", "../dedispersed_profiles/1205955280_I_57.0_flagged_profile.dat"]
for f in files:
    obsid = int(os.path.basename(f)[0:10])
    vals = np.loadtxt(f)
    peak = np.max(vals)
    fluence = np.sum(vals*0.5)
    ind = np.squeeze(np.where(vals==peak))
# Normalise to 154 MHz
    peak = peak * (154. / (1.28 * float(chans[obsid])))**alpha
    fluence = fluence * (154. / (1.28 * float(chans[obsid])))**alpha
# Use whichever peak is greatest
    pair_peak = np.max([pair_peak, peak])
# Sum the fluences
    pair_fluence += fluence
# Set the peak to just between the two obsids
pair_t = obsid - 4.
trunc_times.append(pair_t)
trunc_peaks.append(pair_peak)
trunc_fluences.append(pair_fluence)

times = np.array(times)
peaks = np.array(peaks)
fluences = np.array(fluences)

trunc_times = np.array(trunc_times)
trunc_peaks = np.array(trunc_peaks)
trunc_fluences = np.array(trunc_fluences)

ax.scatter((times-trunc_times[0])/(24.*3600.), peaks, color="blue", marker=".", label="Detections", s=20)
ax.scatter((trunc_times-trunc_times[0])/(24.*3600.), trunc_peaks, color="purple", marker=align_marker(r"$\uparrow$", valign="bottom"), label="Truncated detections", s=20)
ax.scatter((limits[0]-trunc_times[0])/(24.*3600), limits[1], color="red", marker=align_marker(r"$\downarrow$", valign="top"), label="Non-detections", s=20)

ax2.scatter((times-trunc_times[0])/(24.*3600.), fluences, color="blue", marker=".", label="detections", s=20)
ax2.scatter((trunc_times-trunc_times[0])/(24.*3600.), trunc_fluences, color="purple", marker=align_marker(r"$\uparrow$", valign="bottom"), label="truncated detections", s=20)
ax2.scatter((limits[0]-trunc_times[0])/(24.*3600), limits[1], marker=align_marker(r"$\downarrow$", valign="top"), color="red", label="limits", s=20)

before_gap = np.where(times < 1203374696)
after_gap = np.where(times > 1203374696)

ax.plot((np.array(times[before_gap])-trunc_times[0])/(24.*3600.), peaks[before_gap], color="blue", lw=0.5)
ax.plot((np.array(times[after_gap])-trunc_times[0])/(24.*3600.), peaks[after_gap], color="blue", lw=0.5)

ax2.plot((np.array(times[before_gap])-trunc_times[0])/(24.*3600.), fluences[before_gap], color="blue", lw=0.5)
ax2.plot((np.array(times[after_gap])-trunc_times[0])/(24.*3600.), fluences[after_gap], color="blue", lw=0.5)

# Remove unnecessary labels from top axis, and put ticks inside so we can see them
ax.set_xticks([])
ax_dup = ax2.twiny()
ax_dup.set_xlim(ax.get_xlim())
ax_dup.tick_params(axis="x", direction="in")
ax_dup.set_xticklabels([])
plt.subplots_adjust(hspace=0.0)

ax.legend(loc=2, frameon=False, handletextpad=0.1, borderaxespad=0.1, borderpad=0.2)
ax.axhline(0, alpha=1.0, color="black", lw=0.5, ls=":")
ax.set_ylabel("Flux density (Jy)")

ax2.axhline(0, alpha=1.0, color="black", lw=0.5, ls=":")
ax2.set_ylabel("Fluence (Jy s)")
t = Time(int(trunc_times[0]), format="gps")
t.format = "ymdhms"
t = t.value
ax2.set_xlabel(f"Day since first pulse (at {t[0]}-Jan-{t[2]:02d} {t[3]}:{t[4]}:{t[5]:02.0f})")

fig.savefig("brightness_fluence_wrt_time.pdf", bbox_inches="tight")
fig.savefig("brightness_fluence_wrt_time.png", bbox_inches="tight", dpi=1200)
#ax.set_yscale("log")


