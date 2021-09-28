# Dedispersed dynamic spectra and profiles

The scripts in this folder will convert a dynamic spectrum (see the folder `dynamic_spectra` for a few examples) into both dediserpsed dynamic spectra (`[OBSID]_[STOKES]_[DM]_dedispersed_dynspec.dat`), as well as dedispersed profiles (`[OBSID]_[STOKES]_[DM]_flagged_profile.dat`). Due to file size constraints, only a few examples of the former are included here (the ones that correspond to the examples given in `dynamic_spectra`). However, the full set of Stokes I dedispersed profiles are here.

In some observations, individual channels containing RFI (or otherwise contaminated) were manually removed. The list of flagged channels are found in the files `[OBSID]_flagged_chans.txt`. If no such file exists for a given observation, then no channels were flagged.

## Usage

The main script to generate the dedispersed products is `dedispersed_dynspec.py`. Running it without arguments gives the following usage information:

```
python dedispersed_dynspec.py OUTPUT_DAT_FILE BW_MHZ DT

OUTPUT_DAT_FILE must have the following form:
  OBSID_STOKES_DM_dedispersed_dynspec.dat
where
  OBSID  is 10-digit numeric
  STOKES is one of (I, Q, U, V)
  DM     is the desired dedispersion measure in pc/cm^3
BW_MHZ       the single channel bandwidth in MHz
DT           the length of one sample (sec)

This script will look for a file named OBSID_flagged_chans.txt, which is expected
to be a file containing list of channel numbers to flag for that OBSID. If the file
exists, those channels will be flagged. If it doesn't, no channels are flagged
```

For the profiles in this directory, `BW_MHZ` is always set to 0.32 and `DT` to 0.5.
An example run:

```
python dedispersed_dynspec.py 1205001656_I_57.0_dedispersed_dynspec.dat 0.32 0.5
```

**NB:** This script runs the barycentric correction utility found in the `barycentric_correction_utility` folder. Ensure that the files `bc_corr` and `de430.bsp` are compiled and/or available, and that the corresponding path variables at the top of `dedispersed_profiles.py` are set correctly.
