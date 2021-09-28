#ifndef VEC_H
#define VEC_H

typedef struct vec_t
{
    double x;
    double y;
    double z;
} vec;

// Operations on vectors
double dot( vec *v1, vec *v2 );
void   scale( vec *vin, double factor, vec *vout );
double magnitude( vec *v );
void   normalise( vec *vin, vec *vout );
void   cross( vec *v1, vec *v2, vec *vout );
void   cross_norm( vec *v1, vec *v2, vec *vout );
double proj_length( vec *v1, vec *v2 );
void   projection( vec *v1, vec *v2, vec *vout );


#endif
