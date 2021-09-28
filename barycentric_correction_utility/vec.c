#include <math.h>
#include "vec.h"

double dot( vec *v1, vec *v2 )
{
    return v1->x * v2->x +
           v1->y * v2->y +
           v1->z * v2->z;
}

void scale( vec *vin, double factor, vec *vout )
{
    // Multiply vector vin by a scalar factor
    vout->x = vin->x * factor;
    vout->y = vin->y * factor;
    vout->z = vin->z * factor;
}

double magnitude( vec *v )
{
    return sqrt(dot(v, v));
}

void normalise( vec *vin, vec *vout )
{
    // Normalise the vector length
    double mag = magnitude( vin );
    scale( vin, 1.0/mag, vout );
}

void cross( vec *v1, vec *v2, vec *vout )
{
    // Calculate the cross product of two vectors
    vout->x = (v1->y * v2->z) - (v1->z * v2->y);
    vout->y = (v1->z * v2->x) - (v1->x * v2->z);
    vout->z = (v1->x * v2->y) - (v1->y * v2->x);
}

void cross_norm( vec *v1, vec *v2, vec *vout )
{
    // Calculate the normalised cross product of two vectors
    cross( v1, v2, vout );
    normalise( vout, vout );
}

double proj_length( vec *v1, vec *v2 )
{
    // Calculate the projected length of v1 onto v2
    vec vtmp;
    normalise( v2, &vtmp );   // normalise v2 = "v2n"
    return dot( v1, &vtmp );  // calculate (v1 dot v2n)
}

void projection( vec *v1, vec *v2, vec *vout )
{
    // Calculate the projection of v1 onto v2.
    // The result, vout, will be parallel to v2.
    normalise( v2, vout );          // normalise v2 = "v2n"
    double len = dot( v1, vout );   // calculate (v1 dot v2n)
    scale( vout, len, vout );       // calculate (v1 dot v2n) times v2n
}


