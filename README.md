# GLEAM-X_Periodic_Transient
Data and code that support the discovery paper of the long-period GLEAM-X transient

In alphabetical order, the scripts in this repository include:
- barycentric_correction_utility : Correct for the delays introduced to the pulses by the Earth's orbit around the Sun
- brightness_fluence_wrt_time : produce plots of the brightness and fluence of the pulses with respect to time
- burst_periodicity : periodogram analysis to generate initial estimate of source periodicity
- calibration_solutions : radio interferometric calibration solutions that can be applied to MWA measurement sets downloaded from ASVO
- dedispersed_profiles : dispersion corrections applied to produce dynamic spectra
- dynamic_spectra : uncorrected dynamic spectra from each observation
- imaging : example WSClean scripts used to generate images of the transient
- iquv_image : plotting code for Extended Figure 7
- metafits : FITS format information about each observation
- ppdot : Period - Period derivative grid search and plotting code
- pulse_stack : pulse profiles w.r.t. time for each observation, aligned by P and Pdot
- spectral_index : spectral index measurement over five contiguous observations
- xray_limit_estimates : X-ray luminosity and hardness limits from Swift observation
