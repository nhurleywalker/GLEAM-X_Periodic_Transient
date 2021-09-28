#include <stdlib.h>
#include <stdio.h>
#include "SpiceUsr.h"
#include "vec.h"
#include "barycentre.h"

void get_earth_pos( double epoch, char *ephemfile, vec *r_earth )
/* Extracts the Earth's position (km) and velocity (km/s) in J2000 coordinates
 * from the supplied planetary ephemeris.
 */
{
    // Open ephemeris file
    double jd = epoch + 2400000.5;
    SpiceDouble et = (jd - j2000_c() ) * spd_c();
    SpiceDouble state[6];
    SpiceDouble lt;
    furnsh_c( ephemfile );
    spkezr_c( "earth", et, "j2000", "NONE", "solar system barycenter",
              state, &lt );
    unload_c( ephemfile );

    // Put answer into "output" variables
    r_earth->x = state[0];
    r_earth->y = state[1];
    r_earth->z = state[2];
}

double get_bc_correction( double ra_hrs, double dec_deg, double epoch, char *ephemfile )
{
    // Get Earth's position
    vec r_earth;
    get_earth_pos( epoch, ephemfile, &r_earth );

    // Get unit vector in direction of source
    double ra_rad  = HRS2RAD(ra_hrs);
    double dec_rad = DEG2RAD(dec_deg);

    vec rn_src;
    rn_src.x = cos(ra_rad)*cos(dec_rad);
    rn_src.y = sin(ra_rad)*cos(dec_rad);
    rn_src.z = sin(dec_rad);
    vec rn_earth;
    normalise( &r_earth, &rn_earth );

    // Get the barycentric correction
    double proj = proj_length( &r_earth, &rn_src ) * 1.0e3; // km * 1000 = metres
    //fprintf( stderr, "proj_length = %e m\n", proj );

    // Convert the projected length into seconds of light travel time
    return proj/LIGHTSPEED; // In seconds
}
