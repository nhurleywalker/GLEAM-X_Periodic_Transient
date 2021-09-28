/*****************************************************************************
 * bc_corr.c
 * ---------
 *
 * Basic barycentric time-of-arrival correction for a given epoch and look-
 * direction. Uses the CSPICE library for Solar System ephemerides.
 *
 * A lot of this is copied from the veff program:
 * https://github.com/robotopia/veff
 *
 * Sam McSweeney
 * 2021
 * sam.mcsweeney@curtin.edu.au
 ****************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SpiceUsr.h"
#include "vec.h"
#include "barycentre.h"

void usage( char *argv[], FILE *out )
{
    fprintf( out, "usage: %s RA DEC MJD SPK\n"
                  "  RA in J2000 decimal hours\n"
                  "  DEC in J2000 decimal degrees\n"
                  "  MJD in decimal days\n"
                  "  SPK path to NASA planetary ephemeris file\n"
                  "    (e.g. from http://naif.jpl.nasa.gov/pub/naif/"
                  "generic_kernels/spk/planets/de430.bsp)\n",
                  argv[0] );
}

int main( int argc, char *argv[] )
{
    // Parse the command line
    if (argc < 5)
    {
        usage( argv, stderr );
        exit(EXIT_FAILURE);
    }

    double  ra_hrs    = atof(argv[1]);   // Right Ascension in decimal hours
    double  dec_deg   = atof(argv[2]);   // Declination in decimal degrees
    double  epoch     = atof(argv[3]);   // The epoch in question (MJD)
    char   *ephemfile = argv[4];         // Path to the NASA planetary
                                         // ephemeris file

    // Get the correction and print to screen
    double correction = get_bc_correction( ra_hrs, dec_deg, epoch, ephemfile );

    // Print out the result
    printf( "%f\n", correction );

    return EXIT_SUCCESS;
}
