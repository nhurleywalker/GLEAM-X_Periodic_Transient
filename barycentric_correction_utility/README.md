# Barycentric correction utility

## Installation

1. Install CSPICE (https://naif.jpl.nasa.gov/naif/toolkit.html)
2. Edit the Makefile in this directory to ensure that the CSPICE library (`cspice.a`) and the header file (`SpiceUsr.h`) can be found by the compiler (alternatively, place those files in standard locations for your operating system, e.g. `/usr/lib` and `/usr/include`).
3. Run `make`

This will compile the executable `bc_corr` and also download the Solar System ephemeris `de430.bsp`.

## Usage

The following basic usage information can also be obtained by running `bc_corr` without arguments.

```
usage: bc_corr RA DEC MJD SPK
  RA in J2000 decimal hours
  DEC in J2000 decimal degrees
  MJD in decimal days
  SPK path to NASA planetary ephemeris file
    (e.g. from http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp)
```
