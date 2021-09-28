import numpy as np
import sys
import os
import nms
from gps2mjd import gps2mjd

PATH_TO_BC_CORR="../barycentric_correction_utility/bc_corr"
PATH_TO_EPHEMERIS="../barycentric_correction_utility/de430.bsp"

def usage():
    print("python {} OUTPUT_DAT_FILE BW_MHZ DT\n".format(sys.argv[0]))
    print("OUTPUT_DAT_FILE must have the following form:")
    print("  OBSID_STOKES_DM_dedispersed_dynspec.dat")
    print("where")
    print("  OBSID  is 10-digit numeric")
    print("  STOKES is one of (I, Q, U, V)")
    print("  DM     is the desired dedispersion measure in pc/cm^3")
    print("BW_MHZ       the single channel bandwidth in MHz")
    print("DT           the length of one sample (sec)")
    print("")
    print("This script will look for a file named OBSID_flagged_chans.txt, which is expected")
    print("to be a file containing list of channel numbers to flag for that OBSID. If the file")
    print("exists, those channels will be flagged. If it doesn't, no channels are flagged")
    #print("\nNB: This script also does a fractional sample barycentric time correction as well")
 
if __name__ == "__main__":
    # Do some basic parsing of the command line arguments
    if len(sys.argv) < 4:
        usage()
        exit()

    output = sys.argv[1]
    parsed_output = output.split("_")
    obsid  = parsed_output[0] # (Keep as string -- no point converting to int)
    stokes = parsed_output[1]
    dm     = float(parsed_output[2])
    #flo    = float(sys.argv[2])
    bw     = float(sys.argv[2])
    dt     = float(sys.argv[3])

    # Look for the file containing the original dynamic spectrum
    dynspec_filename = '_'.join([obsid, stokes, "dynspec.dat"])
    dynspec = np.loadtxt(dynspec_filename)
    _, _, _, startchan, _, _, _, _, _, ctr_freq = nms.get_processing_meta_info(obsid)

    # Cut off the last 5 seconds' worth of data, since they're often crap.
    dynspec = dynspec[:-10,:]

    # Create a frequency axis
    flo         = (ctr_freq - 30720000/2 + bw/2)*1e-6
    nchans      = dynspec.shape[1]
    freqs       = bw*np.arange(startchan, startchan + nchans) + flo

    # Zero pad the dynamic spectrum to avoid wrapping an "edge" pulse around in time
    zeropad     = 30
    dedispersed, dedispersed_clean = nms.incoh_dedisperse(dynspec, dm, flo, bw, dt, startchan, zeropad=zeropad)
    nsamples    = dedispersed.shape[0]

    # Create a time axis
    times = np.arange(nsamples)*dt - zeropad*dt

    # Apply the barycentric correction, relative to OBSID 1199496944
    MJD = gps2mjd([obsid])
    RA   = 16.4665305555556      # From Natasha's pinned Slack message
    DEC  = -52.5845305555556     # ditto
    cmd = "{} {} {} {} {}".format(PATH_TO_BC_CORR, RA, DEC, MJD[0], PATH_TO_EPHEMERIS)
    #print(cmd)
    stream  = os.popen(cmd)
    bc_corr = float(stream.read())
    times += bc_corr

    # Apply the (bulk) DM correction
    fmid    = np.mean(freqs)
    D = 4.148008e3 # MHz^2 s (pc/cm^3)^-1
    dm_corr = -D*dm/fmid**2
    times += dm_corr


    # Shift every channel by the (same) fractional amount of a sample that
    # corresponds to the barycentric and DM corrections
    #delay        = np.mod(bc_corr + dm_corr, dt)
    #delays       = delay*np.ones((1, nchans))
    #bc_corrected = nms.time_shift(dedispersed, delays, dt)

    # Read in list of channels to flag
    flagged_chans_file = "{}_flagged_chans.txt".format(obsid)
    try:
        flagged_chans = np.loadtxt(flagged_chans_file).astype(int)
        if flagged_chans.size == 1:
            flagged_chans = [flagged_chans,]
        print("Flagging channels {}".format(flagged_chans))
        for k in flagged_chans:
            if k >= 0 and k < nchans:
                dedispersed[:,k]       = np.nan
                dedispersed_clean[:,k] = np.nan
    except:
        pass

    # Write out the dedispersed dynamic spectrum
    header = "Created with\n"
    header += " ".join(sys.argv)
    np.savetxt(output, dedispersed, header=header)

    # Average over frequency (i.e. fscrunch) and write out the profile
    profile = np.array([times, np.nanmean(dedispersed_clean, axis=1)])
    profile_outfile = "{}_{}_{}_flagged_profile.dat".format(parsed_output[0], parsed_output[1], parsed_output[2])
    header += "\n Barycentric correction = {} s".format(bc_corr)
    header += "\n DM correction at centre frequency {} MHz = {} s".format(fmid, dm_corr)
    header += "\n(time - OBSID)(sec)  flux_density(?)"
    np.savetxt(profile_outfile, profile.T, header=header)
