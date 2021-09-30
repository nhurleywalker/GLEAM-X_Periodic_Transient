import numpy as np
import matplotlib.pyplot as plt
import sys

# Dispersion constant
D = 4.148008e3 # MHz^2 s (pc/cm^3)^-1

def dm_delay(DM, freq_MHz, reffreq=None):
    delay = D*DM/freq_MHz**2
    if reffreq is not None:
        delay -= D*DM/reffreq**2

    return delay

def get_profile(obsid, DM, sampling_rate=0.5):
    filename = "../dedispersed_profiles/{}_I_{:.1f}_flagged_profile.dat".format(obsid, DM)
    data = np.loadtxt(filename)
    fluxes = data[:,1]
    times_rel_to_obs = data[:,0]
    return fluxes, times_rel_to_obs

def create_default_time_axis(profile, sampling_rate=0.5):
    npoints = len(profile)
    t_axis = np.arange(npoints)*sampling_rate
    return t_axis

def get_metainfo(filename):
    '''
    Search for the obsid (a string) in processing_meta_info.txt and return
    the parsed line that matches.
    '''
    with open(filename, "r") as processing_file:

        freqs_MHz = {} # A dictionary where key = obsid and value = frequency in MHz
        for line in processing_file:

            if line[0].isnumeric() == False:
                continue

            words = line.split()
            obsid = int(words[0])
            freqs_MHz[obsid] = float(words[10])/1e6

    return freqs_MHz

def apply_dm_corrections(t_axes, DM, freqs_MHz, reffreq=None):

    for obsid in t_axes:
        t_axes[obsid] -= dm_delay(DM, freqs_MHz[obsid], reffreq)

def get_bc_corrections(filename):
    file_contents = np.loadtxt(filename)
    bc_corrections = {int(file_contents[i][0]):file_contents[i][2] for i in range(file_contents.shape[0])}
    return bc_corrections

def apply_bc_corrections(bc_corrections, t_axes):
    for obsid in t_axes.keys():
        t_axes[obsid] += bc_corrections[obsid]

def convert_to_phases(t_axes, P, Pdot=0.0):
    phases = {}
    # Use the first obsid as a reference
    first_obsid = list(t_axes.keys())[0]

    for obsid in t_axes.keys():
        # Get the absolute time (relative to the first obsid)
        t = t_axes[obsid] + (obsid - first_obsid)

        # Model the period and period-derivative to get the "phase" time
        phases[obsid] = np.mod(-0.5*Pdot*(t/P)**2 + (t/P), 1.0)

    return phases

def make_mean_profile(t_axes, fluxes, dt=0.5):

    # Construct a time axis
    first_obsid = list(t_axes.keys())[0]
    tmin = t_axes[first_obsid][0]
    tmax = t_axes[first_obsid][-1]
    for obsid in t_axes.keys():
        if np.min(t_axes[obsid]) < tmin:
            tmin = np.min(t_axes[obsid])
        if np.max(t_axes[obsid]) > tmax:
            tmax = np.max(t_axes[obsid])

    nbins = int(np.round((tmax - tmin)/dt)) + 1
    t = np.arange(nbins)*dt

    # Now bin everything together
    summed = np.zeros(t.shape)
    counts = np.zeros(t.shape)
    for obsid in t_axes.keys():
        bins = np.round((t_axes[obsid] - tmin)/dt).astype(int)
        for i in range(len(bins)):
            if not np.isnan(fluxes[obsid][i]):
                summed[bins[i]] += fluxes[obsid][i]
                counts[bins[i]] += 1
    mean_profile = np.array([summed[i] / counts[i] if counts[i] > 0 else 0.0 for i in range(len(summed))])

    return mean_profile, t

def get_mean_peak(t_axes, fluxes, dt=0.5):
    mean_profile, _ = make_mean_profile(t_axes, fluxes, dt)
    return np.max(mean_profile), mean_profile

if __name__ == '__main__':

    # Get the list of obsids and frequencies from processing_meta_info.txt
    freqs_MHz = get_metainfo("../dedispersed_profiles/processing_meta_info.txt")

    # Get all the profiles
    fluxes = {}
    t_axes = {}
    for obsid in freqs_MHz:
        fluxes[obsid], t_axes[obsid] = get_profile(obsid, 57.0)

    # Apply a weighting due to spectral index
    ref_freq = 154 # MHz
    spec_idx = -1.16
    weighted_fluxes = {obsid:fluxes[obsid]*(ref_freq/freqs_MHz[obsid])**spec_idx for obsid in freqs_MHz}

    # Calculate S/N in P-Pdot parameter space
    P0    = 1091.1558
    Pf    = 1091.1858
    Pstep = 0.0005
    Ps    = np.arange(P0, Pf, Pstep)

    Pdot0    = -4e-9
    Pdotf    = 4e-9
    Pdotstep = 0.1e-9
    Pdots    = np.arange(Pdot0, Pdotf, Pdotstep)

    SNRs = np.zeros((len(Ps),len(Pdots)))
    maxSNR = 0.0
    best_P = 0.0
    best_Pdot = 0.0
    for i in range(len(Ps)):
        P = Ps[i]
        print("Doing period {} ({}/{})".format(P, i+1, len(Ps)))
        for j in range(len(Pdots)):
            Pdot = Pdots[j]
            phases = convert_to_phases(t_axes, P, Pdot=Pdot)
            SNRs[i,j], mean_profile = get_mean_peak(phases, weighted_fluxes, dt=4.5822e-4)
            #print("  SNR = {}".format(SNRs[i,j]))

            if SNRs[i,j] > maxSNR:
                maxSNR = SNRs[i,j]
                best_P = P
                best_Pdot = Pdot
                best_profile = mean_profile

    header = "Created with\n"
    header += " ".join(sys.argv)
    header += "\nBest (P, Pdot) = ({}, {})\n".format(best_P, best_Pdot)

    np.savetxt("p-pdot_search_best_profile.dat", best_profile, header=header)

    header += "Rows = periods (P_start = {}, P_step = {}, Num_Ps = {})\n".format(P0, Pstep, len(Ps))
    header += "Cols = Pdots (Pd_start = {}, Pd_step = {}, Num_Pds = {})\n".format(Pdot0, Pdotstep, len(Pdots))

    np.savetxt("p-pdot_search.dat", SNRs, header=header)

