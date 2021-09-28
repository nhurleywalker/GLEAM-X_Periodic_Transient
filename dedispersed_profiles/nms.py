##########################################
# nms.py
# (Natasha's Mystery Source)
#
# Sam McSweeney, 2021
# sam.mcsweeney@curtin.edu.au
##########################################

from astropy.io import fits
import numpy as np
import re


# The dedispersion constant
D = 4.148008e3 # MHz^2 s (pc/cm^3)^-1

def fitsname(path, obsid, time, chan, stokes):
    '''
    Return the fits filename corresponding to the given obsic, timestep,
    channel number, and stokes parameter.
    '''
    return "{}/{}_time-t{:04d}-{:04d}-{}-image-pb.fits".format(path, obsid, time, chan, stokes)

def MFSfitsname(path, obsid, time, stokes):
    '''
    Return the fits filename corresponding to the given obsic, timestep,
    channel number, and stokes parameter.
    '''
    return "{}/{}_time-t{:04d}-MFS-{}-image-pb.fits".format(path, obsid, time, stokes)

def readfits(fitsfile):
    '''
    Open the fits file and read in the image and frequency information
    '''
    #print("  Opening file {}".format(fitsfile))
    hdul = fits.open(fitsfile)
    image    = np.squeeze(hdul[0].data)
    freq_Hz  = hdul[0].header['CRVAL3']

    return image, freq_Hz

def form_dynspec(path, obsid, starttime, endtime, startchan, endchan, stokes, criterion):
    '''
    Extract the brightest pixel from each image, and use them to form a
    dynamic spectrum.
    '''
    # Construct an empty numpy array to hold the dynamic spectrum
    nchan = endchan - startchan + 1
    ntime = endtime - starttime + 1

    dynspec = np.zeros((ntime, nchan))

    # Loop through the times and chans
    for t in range(ntime):
        # Convert time "idx" to time "number"
        time = t + starttime
        print("Reading time t = {}...".format(time))

        for c in range(nchan):
            # Convert chan "idx" to chan "number"
            chan = c + startchan
            print("  Reading chan c = {}...".format(chan))

            # Get the fits file name corresponding to this time and chan
            fitsfile = fitsname(path, obsid, time, chan, stokes)
            #print("    fitsfile = {}".format(fitsfile))

            # Open the fits file and read in the image
            try:
                image, freq = readfits(fitsfile)
            except:
                print("Could not open {}. Corrupt?".format(fitsfile))
                sys.exit(1)

            # Get the max pixel and put it in the dynamic spectrum
            dynspec[t, c], _ = choose_pixel(image, criterion)

    return dynspec

def choose_pixel(image, criterion="max"):
    '''
    Pick out the value of one single pixel from the given image, according
    to the selected criterion.
    
    criterion can be
      - "max" (default): the value with the largest value is selected
      - "min": the value with the minimum value is selected
      - (i,j): return the (i,j)th pixel
    '''
    if criterion == "max":
        pos = np.unravel_index(image.argmax(), image.shape)
    elif criterion == "min":
        pos = np.unravel_index(image.argmin(), image.shape)
    else:
        try:
            pos = criterion
        except:
            raise("criterion for choosing pixel not one of recognised types")

    val = image[pos[0], pos[1]]

    return val, pos

def time_shift(dynspec, delays, dt):
    # FFT the dynamic spectrum in the time direction
    ffted = np.fft.rfft(dynspec, axis=0)
    nsamples = dynspec.shape[0]

    delays_samples = delays/dt                           # The same delays converted to sample units
    nbins    = ffted.shape[0]                            # The number of bins in the fft
    dph      = 2*np.pi*delays_samples/nsamples           # The delay of the first frequency bin converted to phase (radians)
    dphs, bins = np.meshgrid(dph, np.arange(nbins))
    ramp     = np.exp(1j*dphs*bins)
    ramped   = ramp*ffted
    shifted  = np.fft.irfft(ramped, axis=0)

    return shifted

def incoh_dedisperse(dynspec, dm, flo, bw, dt, startchan_no, zeropad=None):
    # Calculate the phase ramp to apply to each channel
    nsample_orig = dynspec.shape[0]
    nchan  = dynspec.shape[1]
    freqs  = bw*np.arange(startchan_no, startchan_no+nchan) + flo  # A list of the frequencies
    fmid   = np.mean(freqs)                            # The mid (reference) frequency
    delays = D*dm*(1/freqs**2 - 1/fmid**2)             # The delays (sec), relative to the mid frequency

    # Straddle the dynamic spectrum with zeros in the time dimension
    if zeropad is not None:
        if zeropad >= 1:
            nsample = nsample_orig + zeropad*2
            zeropadded = np.zeros((nsample, nchan))
            zeropadded[zeropad:zeropad+nsample_orig,:] = dynspec # Nestle dynspec between the zeros
            dynspec = zeropadded # Rename dynspec

    dedispersed = time_shift(dynspec, delays, dt)

    # Replace any dedispersed pixel that effectively started out as a padded zero with a nan
    dedispersed_clean = np.copy(dedispersed)
    if zeropad is not None:
        if zeropad >= 1:
            for c in range(nchan):
                first = zeropad-int(np.round(delays[c]/dt))
                dedispersed[:first,c] = np.nan;
                dedispersed[first+nsample_orig:,c] = np.nan;

            # Create an alternative spectrum where samples are nan if ANY channel has dedispersed channel
            dedispersed_clean[:zeropad-int(np.round(delays[-1]/dt)),:] = np.nan
            dedispersed_clean[zeropad-int(np.round(delays[0]/dt))+nsample_orig:,:] = np.nan

    return dedispersed, dedispersed_clean


def get_processing_meta_info(obsid):
    '''
    Search for the obsid (a string) in processing_meta_info.txt and return
    the parsed line that matches.
    '''
    with open("processing_meta_info.txt", "r") as processing_file:
        for line in processing_file:
            words = line.split()
            if words[0] == obsid:
                path = words[1]
                starttime = int(words[2])
                endtime   = int(words[3])
                startchan = int(words[4])
                endchan   = int(words[5])
                vmin      = int(words[6])
                vmax      = int(words[7])
                onstart   = int(words[8])
                onstop    = int(words[9])
                ctr_freq  = int(words[10])
                return path, starttime, endtime, startchan, endchan, vmin, vmax, onstart, onstop, ctr_freq

    # If we've got this far, then no obsid was found
    raise("Obsid {} was not found in processing_meta_info.txt".format(obsid))
