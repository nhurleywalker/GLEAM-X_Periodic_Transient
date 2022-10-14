from scipy.io import wavfile as wav
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from glob import glob

arr = np.loadtxt("../../85MHz_aligned.txt")
files = sorted(glob("../../???MHz_aligned.txt"))

for f in files:
    arr = np.vstack([arr, np.loadtxt(f)])

ntimebins = arr.shape[1]

# Print dynamic spectrum, if desired
#plt.imshow(arr, aspect='auto', origin='lower', interpolation='none')
#plt.show()

# To support crossfade, each time sample (worth 0.5 seconds of real data) will be mapped to a wave form that lasts 1 whole second:
# |       |       |       |  - sample boundaries
# |       *-------*       |  - original sample
# |   *---+-------+---*   |  - sample lasting 1 second

iffted_time_per_sample = 1*u.second

# For "ordinary" waves, sample rate is 44100 Hz
sample_rate = 44100 * u.hertz
nyquist_rate = sample_rate/2

iffted_bins_per_sample = int(np.round(iffted_time_per_sample*sample_rate))
Nt = iffted_bins_per_sample # A convenient alias
dt = (1/sample_rate).to(u.second)

# So, the number of frequency bins in each column needs to be the same as the nyquist rate (for rfft).
# Each frequency bin is therefore
df = (1/iffted_time_per_sample).to(u.hertz)

# We want the pulse to be in the audible range, so we want to place the lower frequency such that the 478
# frequency bins of the data don't go *too* high
freq_lo = 300*u.hertz
freq_lo_bin = int(np.round(freq_lo/df))

freq_hi_bin = freq_lo_bin + arr.shape[0]
freq_hi = freq_hi_bin*df

audible_range = [freq_lo, freq_hi]
print("Mapping to frequency range: ", audible_range)

# So now, pad up to the Nyquist rate and pad down to 0 frequency
Nf = int(nyquist_rate/df)
Npad_bins_lo = freq_lo_bin
Npad_bins_hi = Nf - freq_hi_bin
#print(Npad_bins_lo + Npad_bins_hi + arr.shape[0])

padded_arr = np.vstack([np.zeros((Npad_bins_lo, ntimebins)), arr, np.zeros((Npad_bins_hi, ntimebins))])
#plt.imshow(padded_arr, aspect='auto', origin='lower', interpolation='none')
#plt.show()

# Now ifft each column to recover a timeseries
iffted = np.fft.irfft(padded_arr, n=Nt, axis=0)
#plt.imshow(iffted, aspect='auto', origin='lower', interpolation='none')
#plt.show()

# First guess: reshape (without crossfade) and see what it sounds like
waveform = iffted.flatten(order='F')
t = np.arange(Nt*iffted.shape[1])*dt
#plt.plot(t, waveform)
#plt.show()

# Normalise the wave for writing (see https://docs.scipy.org/doc/scipy/reference/generated/scipy.io.wavfile.write.html)
wavmin = np.min(waveform)
wavmax = np.max(waveform)
normalised_waveform = (waveform - wavmin)/(wavmax - wavmin) * 2.0 - 1.0 # Now in range -1 < x < 1

wav.write("test.wav", int(sample_rate.value), waveform)
