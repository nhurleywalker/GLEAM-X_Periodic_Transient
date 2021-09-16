import sys
import numpy as np
from astropy.io import fits
# A least-squares estimator that gives a useable error estimate
from scipy.optimize import leastsq
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter

# Nature requires sans-serif fonts
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica"]})

# http://scipy-cookbook.readthedocs.org/items/FittingData.html
# Define function for calculating a power law
powerlaw = lambda x, amp, index: amp * (x**index)
# define our (line) fitting function
fitfunc = lambda p, x: p[0] + p[1] * x
errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err

def fit_spectrum(freq_array,flux_array,flux_errors):
    pinit = [10.0, -0.7]
    fit = leastsq(errfunc, pinit, args=(freq_array, flux_array, flux_errors), full_output=1)
    covar = fit[1]
    P = fit[0]
    residual = errfunc(P,freq_array, flux_array, flux_errors)
    chi2red = sum(np.power(residual,2))/(len(freqs)-len(pinit))
    alpha=P[1]
    amp = np.exp(P[0])
# Errors
    if covar is not None:
        err_alpha = np.sqrt(covar[1][1])
        err_amp = np.sqrt(covar[0][0])
    else:
        err_alpha = None
        err_amp = None
    return alpha, err_alpha, amp, err_amp, chi2red

freq = 85
linefreqs = []
freqs = np.loadtxt("{0}MHz_freq.txt".format(freq))
scales = np.loadtxt("{0}MHz_scale.txt".format(freq))
x = np.arange(int(len(freqs)/(len(scales)*2)), 9*int(len(freqs)/(len(scales)*2)), int(len(freqs)/(len(scales))))
errs = 0.01 * np.ones(len(scales))
m, err_m, b, err_b, chi2red = fit_spectrum(x, scales, errs)
dat = np.loadtxt("{0}MHz_dedisper.txt".format(freq))
correction = m * np.arange(0, dat.shape[0]) + np.log(b)
correction_2D = np.tile(correction, [dat.shape[1],1]).T
dat *= correction_2D
# Before alignment
dat_b4 = np.loadtxt("{0}MHz_aligned.txt".format(freq))
dat_b4*= correction_2D
for freq in [119, 154, 185, 215]:
    raw_freqs = np.loadtxt("{0}MHz_freq.txt".format(freq))
    linefreqs.append(raw_freqs[0])
    scales = np.loadtxt("{0}MHz_scale.txt".format(freq))
    x = np.arange(int(len(raw_freqs)/(len(scales)*2)), 9*int(len(raw_freqs)/(len(scales)*2)), int(len(raw_freqs)/(len(scales))))
    errs = 0.01 * np.ones(len(scales))
    m, err_m, b, err_b, chi2red = fit_spectrum(x, scales, errs)
# After alignment
    raw_dat = np.loadtxt("{0}MHz_dedisper.txt".format(freq))
    correction = m * np.arange(0, raw_dat.shape[0]) + np.log(b)
    correction_2D = np.tile(correction, [raw_dat.shape[1],1]).T
    raw_dat *= correction_2D
    freqs = np.append(freqs, raw_freqs, 0)
    dat = np.append(dat, raw_dat, 0)
# Before alignment
    raw_b4 = np.loadtxt("{0}MHz_aligned.txt".format(freq))
    raw_b4*= correction_2D
    dat_b4 = np.append(dat_b4, raw_b4, 0)
# Orbcomm satellites yield frequencies between 134 and 139 MHz unuseable: need to zero out
    if freq == 119:
        padfreqs = np.arange((133595000.0+320000.0),(139035000.0+2*320000.0), 320000.0)
        padvals = np.zeros((len(padfreqs), dat.shape[1]))
        freqs = np.append(freqs, padfreqs, 0)
        dat = np.append(dat, padvals, 0)
        dat_b4 = np.append(dat_b4, padvals, 0)
        linefreqs.append(padfreqs[0])

# Making this plot is quite slow so make it optional
makeDyn = True
if makeDyn is True:
    # Time range to display
    tstart = 0
    tend1 = 200
    tend2 = 150
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(3,5), gridspec_kw={'width_ratios': [tend1, tend2]})
    # Before alignment
    ax1.set_ylabel("Frequency (MHz)")
    ax1.pcolormesh(np.arange(tstart, tend1), freqs/1.e6, dat_b4[:,tstart:tend1], vmin=-2, vmax=35., cmap="plasma", edgecolors="face", shading="nearest")
    # After alignment
    ax2.pcolormesh(np.arange(tstart, tend2), freqs/1.e6, dat[:,tstart:tend2], vmin=-2, vmax=35., cmap="plasma", edgecolors="face", shading="nearest")
    ax2.set_yticks([])
    for ax in [ax1, ax2]:
        ax.set_xlabel("Time (seconds)")
        ax.invert_yaxis()
        for l in linefreqs:
            ax.axhline(l/1.e6, color="black", ls="-", lw=1)
    plt.subplots_adjust(wspace=0.05, hspace=0.05)
    fig.savefig("data_2D.png", bbox_inches="tight", dpi=1000)
    fig.savefig("data_2D.pdf", bbox_inches="tight", dpi=1000)

# Transpose
dat = dat.T

# The averaged dedispersed spectrum
avg = np.average(dat, axis=1)

# Need to define RMS from some region without any signal in it
rms = np.nanstd(avg[125:200])

# We're only going to use points where the source is on, i.e. above 3 x the RMS level
cutoff = 3
profile = avg[np.where(avg > cutoff*rms)]
total_profile = np.sum(profile)

mask_1D = avg > cutoff*rms
r = np.arange(0,len(avg))

# Helper plot to show what data is therefore included
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlabel("timestep")
ax.set_ylabel("flux density (Jy/beam)")
ax.fill_between(range(0, len(avg)), -rms, rms, color="red", alpha=0.3)
plt.axhline(cutoff*rms, color="red")
ax.plot(r, avg)
ax.plot(r[mask_1D],avg[mask_1D], color="green")
fig.savefig("average_profile.png", bbox_inches="tight")
fig.savefig("average_profile.pdf", bbox_inches="tight")

# 1D weight profile 
weight = avg / total_profile
weight_tile = np.tile(weight, [dat.shape[1],1]).T
weighted_dat = weight_tile * dat
# But we only want to use the points where the source is on, so mask everything else
# The mask works the opposite way around here
mask_1D = avg < cutoff*rms
mask_2D = np.tile(mask_1D, [dat.shape[1],1]).T
weighted_fluxes = np.ma.sum(np.ma.array(weighted_dat, mask=mask_2D), axis=0)

# Hardcode freqcent to 154MHz
freqcent = 154.e6

# Scipy can silently fail on 32-bit arrays: convert all to float64
# Frequencies -- w.r.t. central frequency
freq_array = np.log(np.array([x/freqcent for x in freqs]),dtype="float64")
flux_array = np.log(np.array(weighted_fluxes, dtype="float64"))
# Make error array rms of each slice, divided by sqrt(number of samples in profile)
err_array = np.std(dat[125:200]/np.sqrt(np.count_nonzero(r[mask_1D])), axis=0)
# Convert to percentage
err_array = err_array / flux_array
err_array = np.sqrt(err_array**2 + (0.05*np.ones(len(err_array)))**2)

# Data with indices < 70 are affected by ionospheric scintillation;
# Data with indices > 382 appear to have a different spectral index
# Not enough data to support fitting a curved spectrum
# Fit will fail if we feed it any NaNs
ind = np.where(weighted_fluxes[70:382] > 0)
alpha, err_alpha, amp, err_amp, chi2red = fit_spectrum(freq_array[70:382][ind],flux_array[70:382][ind],err_array[70:382][ind])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.errorbar(freqs/1.e6, weighted_fluxes, yerr=err_array*weighted_fluxes, lw=0, fmt=".", elinewidth=1, markerfacecolor="gray", markeredgecolor="black", ecolor="gray", zorder=1)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Frequency / MHz")
ax.set_xticks([70, 80, 90, 100, 120, 140, 160, 180, 200, 220])
ax.set_ylabel("Average weighted flux density (Jy/beam)")
ax.yaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
ax.xaxis.set_major_formatter(FormatStrFormatter('%3.0f'))
ax.yaxis.set_minor_formatter(FormatStrFormatter('%3.0f'))
ax.xaxis.set_minor_formatter(FormatStrFormatter('%3.0f'))
ax.set_ylim([1,40])
ax.plot(freqs/1.e6, powerlaw(freqs/freqcent, amp, alpha), color="blue", zorder=10) # Fit
for l in linefreqs:
    ax.axvline(l/1.e6, color="black", ls="-", lw=1)
fig.savefig("spectrum.pdf", bbox_inches="tight")
if err_alpha is None:
    ax.set_title(r"$\alpha = {0:4.2f}$, $\chi^2_r = {1:4.2f}$".format(alpha, chi2red))
else:
    ax.set_title(r"$\alpha = {0:4.2f}\pm{1:4.2f}$, $\chi^2_r = {2:4.2f}$".format(alpha, err_alpha, chi2red))
fig.savefig("spectrum.png", bbox_inches="tight")
