# GLEAM-X_Periodic_Transient
Data and code that support the discovery paper of the long-period GLEAM-X transient

In alphabetical order, the scripts in this repository include:
- barycentric_correction_utility : Correct for the delays introduced to the pulses by the Earth's orbit around the Sun
- brightness_fluence_wrt_time : produce plots of the brightness and fluence of the pulses with respect to time (Fig 2)
- burst_periodicity : periodogram analysis to generate initial estimate of source periodicity
- calibration_solutions : radio interferometric calibration solutions that can be applied to MWA measurement sets downloaded from ASVO
- dedispersed_profiles : dispersion corrections applied to produce dynamic spectra (Extended Fig 1)
- dynamic_spectra : dynamic spectra from each observation (Fig 3)
- imaging : example WSClean scripts used to generate images of the transient
- iquv_image : plotting code for Extended Figure 7
- metafits : FITS format information about each observation
- ppdot : Period - Period derivative grid search and plotting code (Extended Fig 2)
- pulse_stack : pulse profiles w.r.t. time for each observation, aligned by P and Pdot (Fig 1 and Extended Fig 5)
- spectral_index : spectral index measurement over five contiguous observations (Extended Fig 3)
- xray_limit_estimates : X-ray luminosity and hardness limits from Swift observation (Extended Fig 4)

Code for Fig 4 of the paper (transient phase space) can be found in this repository: https://github.com/nhurleywalker/Transient_Phase_Space which is forked from https://github.com/FRBs/Transient_Phase_Space .
