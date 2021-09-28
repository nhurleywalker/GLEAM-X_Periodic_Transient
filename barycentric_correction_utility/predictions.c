#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "SpiceUsr.h"
#include "vec.h"
#include "barycentre.h"

void usage( char *argv[], FILE *out )
{
    fprintf( out, "usage: %s RA DEC PERIOD PULSESTART PULSEEND DM FREQ SPK\n"
                  "  RA in J2000 decimal hours\n"
                  "  DEC in J2000 decimal degrees\n"
                  "  PERIOD in seconds\n"
                  "  Return GPS times from pulse PULSESTART to PULSEEND inclusive\n"
                  "  Correct for DM (pc/cm^{-3}) delay at frequency FREQ (MHz)\n"
                  "  SPK path to NASA planetary ephemeris file\n"
                  "    (e.g. from http://naif.jpl.nasa.gov/pub/naif/"
                  "generic_kernels/spk/planets/de430.bsp)\n",
                  argv[0] );
}

int main ( int argc, char *argv[] )
{
    // Parse the command line
    if (argc < 9)
    {
        usage( argv, stderr );
        exit(EXIT_FAILURE);
    }

    double  ra_hrs     = atof(argv[1]);   // Right Ascension in decimal hours
    double  dec_deg    = atof(argv[2]);   // Declination in decimal degrees
    double  period     = atof(argv[3]);   // Declination in decimal degrees
    int     pulsestart = atoi(argv[4]);   // Starting pulse number
    int     pulseend   = atoi(argv[5]);   // Ending pulse number
    double  DM         = atof(argv[6]);   // Dispersion measure in pc/cm^{-3}
    double  freq       = atof(argv[7]);   // Frequency in MHz
    char   *ephemfile  = argv[8];         // Path to the NASA planetary
                                          // ephemeris file

    // Set the barycentric arrival time of the 0th pulse without dispersion
    double mjdzero = 58120.951773;
    double gpszero = 1198968652.0;
    double mjd, gps;

    // Print output header
    printf( "# Created by\n#    " );
    int i;
    for (i = 0; i < argc; i++)
        printf( " %s", argv[i] );
    printf( "\n# pulse_number  uncorrected_mjd  uncorrected_gps  "
            "bc_correction  dm_correction  corrected_mjd  corrected_gps\n" );

    // Loop through the requested pulses
    int p; // pulse number
    double correction, bc_correction, dm_correction; // Corrections in seconds
    for (p = pulsestart; p <= pulseend; p++)
    {
        // Reset the correction
        correction = 0.0;

        // Add the appropriate number of rotations
        correction += p*period;

        // Add the DM delay
        dm_correction = 4.148008e3 * DM / (freq*freq);
        correction += dm_correction;

        // Add the barycentric delay
        mjd = mjdzero + SEC2MJD(correction);
        bc_correction = -get_bc_correction( ra_hrs, dec_deg, mjd, ephemfile );
        correction   += bc_correction;

        // Print out the result
        mjd = mjdzero + SEC2MJD(correction);
        gps = gpszero + correction;
        printf( "%d %f %d %f %f %f %d\n",
                p,                            // pulse number
                mjdzero + SEC2MJD(p*period),  // uncorrected MJD
                (int)(gpszero + p*period),    // uncorrected GPS
                bc_correction,                // barycentric correction
                dm_correction,                // DM correction
                mjd,                          // corrected MJD
                (int)gps                      // corrected GPS
              );
    }
}
