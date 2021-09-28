#ifndef __BARYCENTRE_H__
#define __BARYCENTRE_H__

#include <math.h>

#define LIGHTSPEED 299792458.0  /* (m/s) */
#define DEG2RAD(x)  ((x)*M_PI/180.0)
#define HRS2RAD(x)  ((x)*M_PI/12.0)
#define SEC2MJD(x)  ((x)/86400.0)

void get_earth_pos( double, char *, vec * );

double get_bc_correction( double, double, double, char * );



#endif
